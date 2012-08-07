--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Search
-- Copyright : (c) [2009..2011] Trevor L. McDonell, Kevin Ying
-- License   : BSD
--
-- An implementation of the SEQUEST algorithm for fast cross-correlation based
-- identification of protein sequences.
--
--
-- References:
--
-- [1] J. K. Eng, B. Fischer, J. Grossmann, and M. J. MacCoss. "A fast sequest
--     cross correlation algorithm." Journal of Proteome Research,
--     7(10):4598-4602, 2008.
--
-- [2] J. K. Eng, A. L. McCormack, and I. John R. Yates. "An approach to
--     correlate tandem mass spectral data of peptides with amino acid sequences
--     in a protein database." Journal of the American Society for Mass
--     Spectrometry, 5(11):976-989, November 1994.
--
--------------------------------------------------------------------------------

module Sequence.Search (searchForMatches) where

import Mass
import Config
import Spectrum
import Spectrum.Data
import Sequence.Match
import Sequence.Fragment
import Sequence.Location
import Sequence.IonSeries
import Util.Misc                                hiding (sublist)
import Util.C2HS
import Util.Time
import Execution                                

import Data.Word
import Data.Maybe
import Data.List                                hiding (lookup)
import System.IO
import Control.Monad
import Prelude                                  hiding (lookup)

import Foreign (sizeOf)
import Foreign.Marshal.Array (peekArray) 
import Foreign.CUDA (DevicePtr, withHostPtr)
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Unboxed            as U
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Util              as CUDA
import qualified Foreign.CUDA.Algorithms        as CUDA
import Foreign.CUDA.Driver.Marshal (getMemInfo)
import qualified Foreign.CUDA.Runtime.Stream    as CUDA
import qualified Foreign.CUDA.BLAS                        as CUBLAS

import Debug.Trace

-- |Used only for search without mods. Candidates within mass limits
type CandidatesByMass = (Int, DevicePtr Word32)
type SpecCandidatesByMass = (Int, Int, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32)

-- |
type CandidatesByModMass = (Int, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32)

-- |Modable Peptides have enough of a acid for a modification to be applied
type CandidatesModable = (Int, {-DevicePtr Word32, -}DevicePtr Word32, DevicePtr Word32, DevicePtr Word32)

-- |Modified peptides from CandidatesModable
type ModCandidates = (Int, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32, Int)

-- |A device pointer to a theoretical spectrum
type IonSeries  = DevicePtr Float

-- |Theoretical spectrum information
type SpecData = ([Int], [Int], [Int], [Int])

-- |A section (after a split) of the database to be searched 
-- (start index to (r,c,n), num elements in section)
type SearchSection = (Int, Int)

--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- |Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> ExecutionPlan -> SequenceDB -> DeviceSeqDB -> HostModInfo -> DeviceModInfo -> [MS2Data] -> IO [(MS2Data, MatchCollection)]
searchForMatches cp ep sdb ddb hmi dmi ms2s = do
  --putTraceMsg $ " dbion length " ++ show (G.length $ dbIon sdb)
  --return []
  (t,matches) <- bracketTime $ searchWithoutMods cp ep sdb ddb ms2s
      --when (verbose cp) $ hPutStrLn stderr ("Search Without Modifications Elapsed time: " ++ showTime t)

      
      --when (verbose cp) $ hPutStrLn stderr ("Search With Modifications Elapsed time: " ++ showTime t2)

  {-results <- forM (zip3 ms2s matches matchesMods $ \(ms2, m, mm) -> do-}
  results <- forM (zip ms2s matches) $ \(ms2, m) -> do
      (t2,mm) <- bracketTime $ case dmi of
                        NoDMod -> return []
                        _      -> searchWithMods cp ep sdb ddb hmi dmi ms2 
      return $ (ms2,sortBy matchScoreOrder $ (mm ++ m))
  return results 

--
-- |Search without modifications
--
searchWithoutMods :: ConfigParams -> ExecutionPlan -> SequenceDB -> DeviceSeqDB -> [MS2Data] -> IO [MatchCollection]
searchWithoutMods cp ep sdb ddb ms2s = 
  filterCandidatesByMass cp ddb masses $ \specCandidatesByMass@(num_spec_cand_total,_,_,_,_) -> 
  if num_spec_cand_total == 0 then return [] else
  mkSpecXCorr ddb specCandidatesByMass ms2s lens $ \specThrys -> do 
    results <- sequestXC cp ep specs specThrys num_spec_cand_total 
    return $ map finish' $ zip3 ms2s specs results 
  
  where
    specs         = map (sequestXCorr cp) ms2s
    lens          = map G.length specs
    masses        = map ms2Mass ms2s
    finish' (ms2,spec,sxis) = 
      let finish (sx,i) = liftM (\f -> Match f sx (sp f) [] []) (lookup sdb i)
          sp f  = matchIonSequence cp (ms2charge ms2) (extractPeaks spec) f
      in  mapMaybe finish sxis

--
-- |Search with modifications
--
searchWithMods :: ConfigParams -> ExecutionPlan -> SequenceDB -> DeviceSeqDB -> HostModInfo -> DeviceModInfo -> MS2Data -> IO MatchCollection
searchWithMods cp ep sdb ddb hmi dmi ms2 = 
  let mass = ms2Mass ms2
  in
  filterCandidatesByModMass cp ddb dmi mass                $ \candsByModMass@(num_cand_mass, _, _, _) -> 
  if num_cand_mass == 0 then return [] else
  filterCandidatesByModability cp ddb dmi candsByModMass  $ \candsByMassAndMod@(num_cand_massmod,_,_,_) -> 
  if num_cand_massmod == 0 then return [] else
  genModCandidates cp ddb dmi candsByMassAndMod $ \modifiedCands ->
  mapMaybe finish `fmap` scoreModCandidates cp ep ddb hmi dmi modifiedCands (ms2charge ms2) (G.length spec) spec
  
  where
    spec          = sequestXCorr cp ms2
    peaks         = extractPeaks spec
    sp            = matchIonSequence cp (ms2charge ms2) peaks
    finish (sx,i,pmod,u) = liftM (\f -> Match (modifyFragment pmod u f) sx (sp f) pmod u) (lookup sdb i)
  
--
-- |Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found.
--
filterCandidatesByMass :: ConfigParams -> DeviceSeqDB -> [Float] -> (SpecCandidatesByMass -> IO b) -> IO b
filterCandidatesByMass cp ddb masses action =
  let d_pep_idx_r_sorted = devResIdxSort ddb
      num_spec           = length masses
  in 
  CUDA.allocaArray num_spec $ \d_spec_begin -> 
  CUDA.allocaArray num_spec $ \d_spec_end -> 
  CUDA.allocaArray num_spec $ \d_spec_num_pep -> 
  CUDA.withVector  (U.fromList masses)   $ \d_spec_masses -> 

  CUDA.findSpecCandsByMass d_spec_begin d_spec_end d_spec_num_pep (devResiduals ddb) d_pep_idx_r_sorted np d_spec_masses num_spec eps >>= \num_spec_cand_total -> 
    action(num_spec_cand_total, num_spec, d_spec_begin, d_spec_end, d_spec_num_pep)
  where
    np  = numFragments ddb
    eps = massTolerance cp

--
-- |For each modification and it's mass delta, find suitable peptides
-- Do this by finding begin and end indices to a list of peptides sorted by residual mass
--
filterCandidatesByModMass :: ConfigParams -> DeviceSeqDB -> DeviceModInfo -> Float -> (CandidatesByModMass -> IO b) -> IO b
filterCandidatesByModMass cp ddb dmi mass action =
  let (num_mod, _, _, d_mod_delta) = devModCombs dmi 
      d_pep_idx_r_sorted = devResIdxSort ddb
  in
  CUDA.allocaArray num_mod $ \d_begin ->      -- begin idx to d_pep_idx_r_sorted
  CUDA.allocaArray num_mod $ \d_end ->        -- end idx
  -- CUDA.allocaArray num_mod $ \d_num_pep ->    -- end[i] - begin[i]
  CUDA.allocaArray num_mod $ \d_num_pep_scan-> 
  CUDA.findBeginEnd d_begin d_end {-d_num_pep-} d_num_pep_scan (devResiduals ddb) d_pep_idx_r_sorted (numFragments ddb) d_mod_delta num_mod mass eps >>= \num_cand_mass -> action (num_cand_mass, d_begin, d_end, d_num_pep_scan)
  where
    eps = massTolerance cp

--
-- |Search a subset of peptides for peptides which a specific combination of modifications 
-- can be applied.
-- Peptides are deemed modifable if it has enough amino acids.
--
-- This generates a device array containing the indices of the relevant
-- peptides, together with the number of candidates found, and for each candidate the number of
-- modifiable acids to be generated from each peptide
--
filterCandidatesByModability :: ConfigParams -> DeviceSeqDB -> DeviceModInfo -> CandidatesByModMass -> (CandidatesModable -> IO b) -> IO b
filterCandidatesByModability cp ddb dmi (num_cand_mass, d_begin, d_end, d_num_pep_scan) action =
  let (num_ma, d_ma, _)                 = devModAcids dmi 
      (num_mod, d_mod_ma_count, _, _)   = devModCombs dmi 
      d_pep_idx_r_sorted                = devResIdxSort ddb
  in
  CUDA.allocaArray num_cand_mass $ \d_pep_idx ->        -- idx of peptide to tc tn
  CUDA.allocaArray num_cand_mass $ \d_pep_mod_idx ->    -- peptides corresponding modification
  CUDA.allocaArray (num_cand_mass*num_ma) $ \d_pep_ma_count -> -- modable acid count for each peptide
  CUDA.findModablePeptides d_pep_idx d_pep_mod_idx d_pep_ma_count num_cand_mass (devIons ddb) (devTerminals ddb) d_pep_idx_r_sorted d_begin d_end d_num_pep_scan d_mod_ma_count num_mod d_ma num_ma >>= \n -> do
    action (n, d_pep_idx, d_pep_mod_idx, d_pep_ma_count)

--
-- |Given some unmodified peptides, modifications and
-- count of each modifiable aa in each pep, generate modified peptides
-- See gen_mod_pep_cu for details
-- 
genModCandidates :: ConfigParams -> DeviceSeqDB -> DeviceModInfo -> CandidatesModable -> (ModCandidates -> IO b) -> IO b
genModCandidates cp ddb dmi (num_cand_massmod, d_pep_idx, d_pep_mod_idx, d_pep_ma_count) action =
  let (num_ma, _, _) = devModAcids dmi 
      (num_mod, d_mod_ma_count, d_mod_ma_count_sum, _) = devModCombs dmi 
  in

  -- calc number of modified peps generated from each pep
  CUDA.mallocArray num_cand_massmod          >>= \d_pep_num_mpep ->         -- number of modified peptides for each peptide
  CUDA.mallocArray (num_cand_massmod*num_ma) >>= \d_pep_ma_num_comb_scan -> -- above scaned for each peptide
  CUDA.calcTotalModCands d_pep_num_mpep d_pep_ma_num_comb_scan d_mod_ma_count d_pep_idx d_pep_mod_idx d_pep_ma_count num_cand_massmod num_mod num_ma >>= \num_mpep -> 

  CUDA.mallocArray num_mpep >>= \d_mpep_rank ->               
  CUDA.mallocArray num_mpep >>= \d_mpep_ith_cand ->           

  CUDA.allocaArray num_mpep $ \d_mpep_pep_idx -> --
  CUDA.allocaArray num_mpep $ \d_mpep_pep_mod_idx -> --
  CUDA.allocaArray num_mpep $ \d_mpep_mod_ma_count_sum_scan -> --

  CUDA.prepareGenMod d_mpep_pep_idx d_mpep_pep_mod_idx d_mpep_rank d_mpep_ith_cand d_mpep_mod_ma_count_sum_scan d_mod_ma_count_sum d_pep_idx d_pep_mod_idx d_pep_num_mpep num_cand_massmod num_mpep >>= \ len_unrank ->

  CUDA.allocaArray (len_unrank) $ \d_mpep_unrank -> do --
    CUDA.genModCands d_mpep_unrank d_mod_ma_count d_mpep_ith_cand d_mpep_rank d_mpep_mod_ma_count_sum_scan d_pep_mod_idx d_pep_ma_count d_pep_ma_num_comb_scan num_mpep num_ma 

    -- free mallocs
    CUDA.free d_pep_num_mpep
    CUDA.free d_pep_ma_num_comb_scan
    CUDA.free d_mpep_rank
    CUDA.free d_mpep_ith_cand

    action (num_mpep, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank)

--
-- |Score modified peptides by generating spectrums an then scoring them.
-- Will calculate remaining memory and split the work so device will not run out of memory
--
scoreModCandidates :: ConfigParams -> ExecutionPlan -> DeviceSeqDB -> HostModInfo -> DeviceModInfo -> ModCandidates -> Float -> Int -> Spectrum -> IO [(Float,Int,PepMod,[Int])]
scoreModCandidates cp ep ddb hmi dmi mcands chrg len expr = 
  let (num_mpep, _, _, _, _, _) = mcands

      -- each modified peptide will need the same amount of memory
      memPerPep = len * sizeOf (undefined :: Word32) -- d_mspec
                + sizeOf (undefined :: Word32) -- d_mpep_idx
                + sizeOf (undefined :: Float) -- d_score
      totalReqMem = num_mpep * memPerPep
  in
  --if num_mpep == 0 then return [] else 

  CUDA.withVector expr $ \d_expr -> 
  getMemInfo >>= \(freeMem,_) ->
  
  -- predict memory needed and split work accordingly
  let ratio    = (fromIntegral totalReqMem) / (fromIntegral freeMem) :: Double
      split    = (if ratio < 1 then 1 else 1 + ceiling ratio) -- +1 to stay well within limit
      sections = sectionsOfSplit num_mpep split
      most     = ceiling $ ((fromIntegral num_mpep) / (fromIntegral split) :: Double)-- most number in a section
  in

  CUDA.allocaArray (len*most) $ \d_mspecThry -> 
  CUDA.allocaArray most       $ \d_score -> 
  CUDA.allocaArray most       $ \d_mpep_idx -> 
    {-when (verbose cp) $ hPutStrLn stderr-}
                      {-$ "Free memory " ++ show freeMem ++ "\n"-}
                        {-++ "Memory required per peptide " ++ show memPerPep ++ "\n"-}
                        {-++ "Total memory required " ++ show totalReqMem ++ "\n"-}
                        {-++ "Total / Free " ++ show ratio ++ "\n"-}
                        {-++ "Number jobs to split into " ++ show split ++ "\n"-}
    liftM concat $ forM sections $ \s -> do
      mkModSpecXCorr ddb dmi mcands chrg len d_mspecThry s
      sequestXCMod cp ep hmi mcands expr d_expr d_mspecThry d_score d_mpep_idx s


--
-- |Generate theoretical spectrums
--
mkModSpecXCorr :: DeviceSeqDB -> DeviceModInfo -> ModCandidates -> Float -> Int -> IonSeries -> SearchSection -> IO ()
mkModSpecXCorr db dmi (num_mpep, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank) chrg len d_mspec (start, nsearch) =
  let (num_ma, d_ma, d_ma_mass) = devModAcids dmi 
      (_, d_mod_ma_count, _, d_mod_delta) = devModCombs dmi 
  in do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_mspec
    CUDA.memset  d_mspec bytes 0
    CUDA.addModIons d_mspec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_mpep_pep_idx d_mpep_pep_mod_idx d_mpep_unrank d_mpep_mod_ma_count_sum_scan len_unrank num_mpep d_mod_ma_count d_mod_delta d_ma d_ma_mass num_ma len ch start nsearch
  where
    n = len * nsearch
    ch = round $ max 1 (chrg - 1)

--
-- |Generate a theoretical spectral representation for each of the specified
-- candidates. This generates spectral peaks for all of the A-, B- and Y-ions of
-- the given sequences, retaining only the most intense peak in each bin.
--
-- On the device, this is stored as a dense matrix, each row corresponding to a
-- single sequence.
--
mkSpecXCorr :: DeviceSeqDB -> SpecCandidatesByMass -> [MS2Data] -> [Int] -> ((SpecData, IonSeries) -> IO b) -> IO b
mkSpecXCorr ddb specCands ms2s lens action =
  let (num_spec_cand_total, num_spec, 
       d_spec_begin, d_spec_end, d_spec_num_pep) = specCands
      chrgs = map ms2charge ms2s
      d_pep_idx_r_sorted = devResIdxSort ddb
  in do
    spec_num_pep' <- CUDA.peekListArray num_spec d_spec_num_pep
    spec_begin'   <- CUDA.peekListArray num_spec d_spec_begin

    let spec_num_pep = map fromIntegral spec_num_pep'
        spec_begin   = map fromIntegral spec_begin'
        spec_sum_len_scan = scanl (+) 0 $ map (\(a,b) -> a * b) (zip lens spec_num_pep) 
        total_len    = last spec_sum_len_scan 
        max_chs = map (\chrg -> round $ max 1 (chrg - 1)) chrgs

    CUDA.allocaArray total_len $ \d_specs -> do
      let bytes = fromIntegral $ total_len * CUDA.sizeOfPtr d_specs
      CUDA.memset  d_specs bytes 0
      forM_ (zip5 [0..num_spec-1] spec_num_pep spec_sum_len_scan max_chs lens) $ 
        \(spec_idx, num_idx, thry_start, chrg, len) -> do
          let d_idx = d_pep_idx_r_sorted `CUDA.advanceDevPtr` (spec_begin !! spec_idx) 
              d_specs' = d_specs `CUDA.advanceDevPtr` thry_start 
          
          CUDA.addIons d_specs' (devResiduals ddb) (devMassTable ddb) (devIons ddb) (devTerminals ddb) d_idx num_idx chrg len
      CUDA.sync
      action ((lens, spec_sum_len_scan, spec_num_pep, spec_begin), d_specs)


--
-- |Score each candidate sequence against the observed intensity spectra,
-- returning the most relevant results.
--
sequestXC :: ConfigParams -> ExecutionPlan -> [Spectrum] -> (SpecData, IonSeries) -> Int -> IO [[(Float,Int)]]
sequestXC cp ep exprs ((spec_lens, spec_sum_len_scan, spec_num_pep, spec_begin),d_thry) cand_total = 
  let n' = max (numMatches cp) (numMatchesDetail cp) 
      all_exprs = U.concat exprs
      spec_len_scan = scanl (+) 0 spec_lens
      spec_num_pep_scan = scanl (+) 0 spec_num_pep
      num_spec = length exprs
      spec_retrieve = map (\num_pep -> min n' num_pep) spec_num_pep
      spec_retrieve_scan = scanl (+) 0 spec_retrieve
      retrieve_total = last spec_retrieve_scan
  in if cand_total == 0 then return [] else
  CUDA.withVector  all_exprs $ \d_all_exprs ->
  CUDA.allocaArray cand_total $ \d_spec_cand_idx ->
  CUDA.allocaArray cand_total $ \d_scores -> 
  CUDA.withHostArray retrieve_total $ \h_spec_cand_idx_sorted ->
  CUDA.withHostArray retrieve_total $ \h_scores -> do
    -- For each spectrum score each candidate sequence
    -- Hopefully thec multiplication will run concurrently
    forM_ (zip7 [0..num_spec-1] (cudaStreams ep) spec_sum_len_scan spec_len_scan spec_num_pep_scan spec_num_pep spec_lens) $ 
      \(spec_idx, strm, thry_start, spec_start, score_start, num_pep, len) -> do
        let m        = num_pep 
            n        = len
            d_thry'  = d_thry `CUDA.advanceDevPtr` thry_start
            d_expr   = d_all_exprs `CUDA.advanceDevPtr` spec_start 
            d_score' = d_scores `CUDA.advanceDevPtr` score_start 

        -- Set the stream for CUBLAS, concurrent kernel execution works for compute 2.0 and up
        CUBLAS.setStream (cublasHandle ep) strm
        -- Because cublas uses col major storage (as opposed to row major) swap row and col values and use CUBLAS_OP_T (Transform)
        --CUDA.mvm   (cublasHandle ep) d_score d_thry d_expr m n
        CUBLAS.sgemv (cublasHandle ep) (CUBLAS.T) n m 1 d_thry' n d_expr 1 0 d_score' 1

    -- Sort the results
    -- and copy asynchronously, thrust sort will not run concurrently,
    -- so copying asynchronously will only reall have small speed up
    forM_ (zip6 spec_retrieve spec_retrieve_scan (cudaStreams ep) spec_num_pep_scan spec_num_pep spec_begin) $ 
      \(n, ret_start, strm, score_start, num_pep, sorted_begin) -> do
        CUDA.block strm
        let d_score = d_scores `CUDA.advanceDevPtr` score_start
            h_score = h_scores `CUDA.advanceHostPtr` ret_start
            d_idx    = d_spec_cand_idx `CUDA.advanceDevPtr` score_start
            h_idx    = h_spec_cand_idx_sorted `CUDA.advanceHostPtr` ret_start

        CUDA.rsort_idx d_score d_idx num_pep sorted_begin

        -- Retrieve the most relevant matches
        --
        CUDA.peekArrayAsync n d_score h_score (Just strm)
        CUDA.peekArrayAsync n d_idx h_idx (Just strm)

    -- retrieve results
    results <- forM (zip3 spec_retrieve spec_retrieve_scan (cudaStreams ep)) $ 
      \(n, re_start, strm) -> do
        CUDA.block strm
        let h_score = h_scores `CUDA.advanceHostPtr` re_start 
            h_idx    = h_spec_cand_idx_sorted `CUDA.advanceHostPtr` re_start

        sc <- withHostPtr h_score $ \ptr -> peekArray n ptr
        sorted_idx <- withHostPtr h_idx $ \ptr -> peekArray n ptr
        let ix = map (\i -> resIdxSort ep U.! i) $ map fromIntegral sorted_idx 

        return $ zipWith (\s i -> (s/10000, i)) sc ix
    return results 

--
-- |Modified candidates version
--
sequestXCMod :: ConfigParams -> ExecutionPlan -> HostModInfo -> ModCandidates -> Spectrum -> DevicePtr Float -> IonSeries -> DevicePtr Float -> DevicePtr Word32 -> SearchSection -> IO [(Float,Int,PepMod,[Int])]
sequestXCMod cp ep hmi (_, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank) expr d_expr d_thry d_score d_mpep_idx (start, nsearch) = 
  let num_mpep = nsearch
      (num_ma, ma, ma_mass) = modAcids hmi
      (_, mod_ma_count, mod_ma_count_sum, _) = modCombs hmi

      n' = max (numMatches cp) (numMatchesDetail cp) 
      h_mpep_idx = G.generate num_mpep (\i -> cIntConv i) :: (U.Vector Word32) 
  in
    -- There may be no candidates as a result of bad database search parameters,
    -- or if something unexpected happened (out of memory)
    --
  if num_mpep == 0 then return [] else do
      
    CUDA.copyVector h_mpep_idx d_mpep_idx

    -- Score and rank each candidate sequence
    --
    let m = num_mpep
        n = G.length expr

    --CUDA.mvm   (cublasHandle ep) d_score d_thry d_expr num_mpep (G.length expr)

    -- Because cublas uses col major storage (as opposed to row major) swap row and col values and use CUBLAS_OP_T 
    CUBLAS.sgemv (cublasHandle ep) (CUBLAS.T) n m 1 d_thry n d_expr 1 0 d_score 1
    CUDA.rsort d_score d_mpep_idx num_mpep

    -- Retrieve the most relevant matches
    --
    let n = min n' num_mpep
    score             <- CUDA.peekListArray n d_score
    mpep_idx'         <- CUDA.peekListArray n d_mpep_idx
    mpep_pep_idx'     <- CUDA.peekListArray num_mpep (d_mpep_pep_idx `CUDA.advanceDevPtr` start)
    mpep_pep_mod_idx' <- CUDA.peekListArray num_mpep (d_mpep_pep_mod_idx `CUDA.advanceDevPtr` start)
    mpep_unrank'      <- CUDA.peekListArray len_unrank d_mpep_unrank
    mpep_mod_ma_count_sum_scan' <- CUDA.peekListArray num_mpep (d_mpep_mod_ma_count_sum_scan `CUDA.advanceDevPtr` start)

    let mpep_idx         = map fromIntegral mpep_idx'
        mpep_pep_idx     = map fromIntegral mpep_pep_idx'
        mpep_unrank      = map fromIntegral mpep_unrank'
        mpep_pep_mod_idx = map fromIntegral mpep_pep_mod_idx'
        mpep_mod_ma_count_sum_scan = map fromIntegral mpep_mod_ma_count_sum_scan'

        pack s i =
            let mi = mpep_pep_mod_idx !! i
            in (s/10000, 
                (mpep_pep_idx !! i), 
                (getPepMod mi), 
                (getUnrank i mi))
        getPepMod mi = 
            let ma_count = sublist (mi*num_ma) num_ma mod_ma_count
            in  U.toList $ U.filter (\(_,x,_) -> x /= 0) $
                           U.zipWith3 (\a b c -> (a, b, c)) ma ma_count ma_mass
        
        getUnrank i mi = U.toList $ sublist (mpep_mod_ma_count_sum_scan !! i) 
                                    (mod_ma_count_sum U.! mi) $ U.fromList mpep_unrank

    return $ zipWith pack score mpep_idx
    where
    sublist beg num ls = U.drop beg $ U.take (beg + num) ls
