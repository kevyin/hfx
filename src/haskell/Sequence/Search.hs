--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Search
-- Copyright : (c) [2009..2010] Trevor L. McDonell
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
import Sequence.Match
import Sequence.Fragment
import Sequence.Location
import Sequence.IonSeries
import Util.Misc
import Util.C2HS
import Util.Time
import Util.Show

import Data.Word
import Data.Maybe
import Data.List                                hiding (lookup)
import Control.Monad
import System.IO
import Prelude                                  hiding (lookup)

import Foreign.CUDA (DevicePtr)
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Unboxed            as U
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Util              as CUDA
import qualified Foreign.CUDA.Algorithms        as CUDA
--import qualified Foreign.CUDA.Stream            as CUDA

import Debug.Trace

type CandidatesByMass = (Int, DevicePtr Word32)
type CandidatesByModMass = (DevicePtr Word32, DevicePtr Word32, Int, DevicePtr Word32, DevicePtr Word32)
type CandidatesModable = (Int, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32)
type ModCandidates = (Int, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32, DevicePtr Word32, Int)
type IonSeries  = DevicePtr Word32

type PepModDevice = (Int, DevicePtr Word8, DevicePtr Word8, Int, DevicePtr Float)
--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> SequenceDB -> DeviceSeqDB -> HostModInfo -> DeviceModInfo -> MS2Data -> IO MatchCollection
searchForMatches cp sdb ddb hmi dmi ms2 = do
  matches <- searchWithoutMods cp sdb ddb ms2 
  --(t,matches) <- bracketTime $ searchWithoutMods cp sdb ddb ms2 
  --when (verbose cp) $ hPutStrLn stderr ("searchWithoutMods Elapsed time: " ++ showTime t)

  matchesMods <- searchWithMods cp sdb ddb hmi dmi ms2 
  --(t,matchesMods) <- bracketTime $ searchWithMods cp sdb ddb hmi dmi ms2 
  --when (verbose cp) $ hPutStrLn stderr ("searchWithMods Elapsed time: " ++ showTime t)

  return $ sortBy matchScoreOrder $ matches ++ matchesMods
  --return $ sortBy matchScoreOrder $ matches

searchWithoutMods :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> IO MatchCollection
searchWithoutMods cp sdb ddb ms2 =
  filterCandidateByMass cp ddb mass                                   $ \candidatesByMass ->
  mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec)          $ \specThry   -> do
  mapMaybe finish `fmap` sequestXC cp candidatesByMass spec specThry
  where
    --modComb = (ma, ma_count, 1434.7175)
    --ma = map c2w ['A','C']
    --ma_count = [0,4]
    spec          = sequestXCorr cp ms2
    peaks         = extractPeaks spec

    finish (sx,i) = liftM (\f -> Match f sx (sp f) [] []) (lookup sdb i)
    sp            = matchIonSequence cp (ms2charge ms2) peaks
    mass          = (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 - 1) * massH) - (massH + massH2O)

--
--
searchWithMods :: ConfigParams -> SequenceDB -> DeviceSeqDB -> HostModInfo -> DeviceModInfo -> MS2Data -> IO MatchCollection
searchWithMods cp sdb ddb hmi dmi ms2 = 
  let mass = (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 - 1) * massH) - (massH + massH2O)
  in
  -- marshall over d_ma, d_ma_mass, num_ma,
  --               d_mod_ma_count
  --          calc d_mod_ma_count_sum
  --          calc d_mod_delta
  --          sort candidates by residual
  --
  -- findCandidates
  filterCandidateByModMass cp ddb dmi mass $ \candsByModMass -> 
  --    by mass
  --          find beg and end for each modcomb
  --          calc num_pep for each modcomb
  --          calc above scanned
  filterCandidatesByModability cp ddb dmi candsByModMass $ \candsByMassAndMod ->
  --    by modability
  --          alloc mem 
  --          calc pep_ma_count
  --               pep_idx
  --               pep_mod_idx
  --          filter pep_cand_idx
  genModCandidates cp ddb dmi candsByMassAndMod $ \modifiedCands ->
  -- generate modified candidates
  --          calc num_mpep for filtered cands
  --          calc ma_ncomb
  --          calc na_ncomb_scan
  --          calc mpep_unrank
  mkModSpecXCorr ddb dmi modifiedCands (ms2charge ms2) (G.length spec)  $ \mspecThry   -> do
  --mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec)              $ \specThry   -> -- @TODO delete later
  mapMaybe finish `fmap` sequestXCMod cp hmi dmi modifiedCands spec mspecThry -- @TODO
  where
    spec            = sequestXCorr cp ms2
    peaks           = extractPeaks spec
    finish (sx,i,pmod,u) = liftM (\f -> Match (modifyFragment pmod u f) sx (sp f) pmod u) (lookup sdb i)
    sp              = matchIonSequence cp (ms2charge ms2) peaks
    --pmod            = zipWith3 (\a b c -> (w2c a, fromIntegral b, c)) ma ma_count ma_mass
  
--searchAModComb :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> ([Word8], [Word8], Float) -> IO MatchCollection
--searchAModComb cp sdb ddb ms2 (ma, ma_count, mod_mass) = 
  ----traceShow ("Searching comb: ", pmod) $

  --filterCandidateByMass cp ddb mod_mass                                     $ \candidatesByMass ->
  --CUDA.withVector (U.fromList ma)       $ \d_mod_ma       ->   -- modifable acid
  --CUDA.withVector (U.fromList ma_count) $ \d_mod_ma_count ->   -- number of the acid to modify
  --CUDA.withVector (U.fromList ma_mass)  $ \d_mod_ma_mass  ->   -- 
  --let d_mods = (mod_num_ma, d_mod_ma, d_mod_ma_count, sum_ma_count, d_mod_ma_mass) in
  --filterCandidatesByModability cp ddb candidatesByMass d_mods                    $ \candidatesByMassAndMod ->
  --genModCandidates cp ddb candidatesByMassAndMod d_mods                         $ \modifiedCandidates ->
  --mkModSpecXCorr ddb d_mods modifiedCandidates (ms2charge ms2) (G.length spec)  $ \mspecThry   -> do
  ----mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec)              $ \specThry   -> -- @TODO delete later
  --mapMaybe finish `fmap` sequestXCMod cp modifiedCandidates sum_ma_count spec mspecThry -- @TODO
  ----(t,res) <- bracketTime $ mapMaybe finish `fmap` sequestXCMod cp modifiedCandidates sum_ma_count spec mspecThry -- @TODO
  ----when (verbose cp) $ hPutStrLn stderr ("sequestXCMod Elapsed time: " ++ showTime t)
  ----return res
  --where
    --ma_mass         = getMA_Mass cp 
    --sum_ma_count    = fromIntegral $ sum ma_count
    --mod_num_ma      = length ma
    --spec            = sequestXCorr cp ms2
    --peaks           = extractPeaks spec

    --finish (sx,i,u) = liftM (\f -> Match (modifyFragment pmod u f) sx (sp f) pmod u) (lookup sdb i)
    --sp              = matchIonSequence cp (ms2charge ms2) peaks
    --pmod            = zipWith3 (\a b c -> (w2c a, fromIntegral b, c)) ma ma_count ma_mass
    

--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found.
--
filterCandidateByMass :: ConfigParams -> DeviceSeqDB -> Float -> (CandidatesByMass -> IO b) -> IO b
filterCandidateByMass cp db mass action =
  CUDA.allocaArray np $ \d_idx -> -- do
  CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (n,d_idx)
  --(t,n) <- bracketTime $ CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) 
  --when (verbose cp) $ hPutStrLn stderr ("filterCandidatesByMass Elapsed time: " ++ showTime t)
  --action (n,d_idx)
  where
    np    = numFragments db
    delta = massTolerance cp

--
-- Searches for each mass made by modified candidates
filterCandidateByModMass :: ConfigParams -> DeviceSeqDB -> DeviceModInfo -> Float -> (CandidatesByModMass -> IO b) -> IO b
filterCandidateByModMass cp ddb dmi mass action =
  let -- (num_ma, d_ma, d_ma_mass) = devModAcids dmi 
      (num_mod, _, _, d_mod_delta) = devModCombs dmi 
      d_pep_idx_r_sorted = devResIdxSort dmi
  in
  CUDA.allocaArray num_mod $ \d_begin ->
  CUDA.allocaArray num_mod $ \d_end ->
  CUDA.allocaArray num_mod $ \d_num_pep -> do
  CUDA.allocaArray num_mod $ \d_num_pep_scan-> do
  CUDA.findBeginEnd d_begin d_end d_num_pep d_num_pep_scan (devResiduals ddb) d_pep_idx_r_sorted (numFragments ddb) d_mod_delta num_mod mass eps >>= \num_pep_total -> do
  action (d_begin, d_end, num_pep_total, d_num_pep, d_num_pep_scan)
  where
    --np    = numFragments db
    eps = massTolerance cp

--
-- Search a subset of peptides for peptides which a specific combination of modifications 
-- can be applied.
-- Peptides are deemed modifable if it has enough amino acids.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found, and for each candidate, a count of the
-- modifiable acids in the peptide
--
--
filterCandidatesByModability :: 
                       ConfigParams 
                    -> DeviceSeqDB 
                    -> DeviceModInfo
                    -> CandidatesByModMass     --- ^ subset of peptides, [idx to pep], number of pep
                    -> (CandidatesModable -> IO b)
                    -> IO b
filterCandidatesByModability cp ddb dmi (d_begin, d_end, num_pep_total, d_num_pep, d_num_pep_scan) action =
  let (num_ma, d_ma, _) = devModAcids dmi 
      (num_mod, d_mod_ma_count, _, _) = devModCombs dmi 
      d_pep_idx_r_sorted = devResIdxSort dmi
  in
  CUDA.allocaArray num_pep_total $ \d_pep_valid_idx-> -- filtered peptide indices to be returned here
  CUDA.allocaArray num_pep_total $ \d_pep_idx -> -- idx of peptide to tc tn
  CUDA.allocaArray num_pep_total $ \d_pep_mod_idx -> -- peptides corresponding modification
  CUDA.allocaArray (num_pep_total*num_ma) $ \d_pep_ma_count -> 
  CUDA.findModablePeptides d_pep_valid_idx d_pep_idx d_pep_mod_idx d_pep_ma_count num_pep_total (devIons ddb) (devTerminals ddb) d_pep_idx_r_sorted d_begin d_end d_num_pep_scan d_mod_ma_count num_mod d_ma num_ma >>= \n -> do
    action (n, d_pep_valid_idx, d_pep_idx, d_pep_mod_idx, d_pep_ma_count)


--
-- Given some unmodified peptides, a modification and
-- count of each modifiable aa in each pep, generate modified peptides
-- 
genModCandidates :: ConfigParams
                 -> DeviceSeqDB
                 -> DeviceModInfo
                 -> CandidatesModable
                 -> (ModCandidates -> IO b)
                 -> IO b
genModCandidates cp ddb dmi (num_pep, d_pep_valid_idx, d_pep_idx, d_pep_mod_idx, d_pep_ma_count) action =
  let (num_ma, _, _) = devModAcids dmi 
      (num_mod, d_mod_ma_count, d_mod_ma_count_sum, _) = devModCombs dmi 
  in

  -- calc number of modified peps generated from each pep
  --
  CUDA.allocaArray num_pep $ \d_pep_num_mpep -> 
  CUDA.allocaArray (num_pep*num_ma) $ \d_pep_ma_num_comb -> 
  CUDA.allocaArray (num_pep*num_ma) $ \d_pep_ma_num_comb_scan -> 
  CUDA.calcTotalModCands d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan d_mod_ma_count d_pep_idx d_pep_mod_idx d_pep_ma_count d_pep_valid_idx num_pep num_mod num_ma >>= \num_mpep -> 

  --CUDA.sum_Word32 d_mod_ma_count_sum num_mod >>= \ma_count_sum_total ->
  --CUDA.scan d_mod_ma_count_sum_scan d_mod_ma_count_sum num_mod >>= \ma_count_sum_total ->
  ----action num_mpep
  ----
  CUDA.allocaArray num_mpep $ \d_mpep_pep_idx -> 
  CUDA.allocaArray num_mpep $ \d_mpep_pep_mod_idx -> 
  CUDA.allocaArray num_mpep $ \d_mpep_rank -> 
  CUDA.allocaArray num_mpep $ \d_mpep_ith_valid-> 
  CUDA.allocaArray num_mpep $ \d_mpep_mod_ma_count_sum -> 
  CUDA.allocaArray num_mpep $ \d_mpep_mod_ma_count_sum_scan -> 

  CUDA.prepareGenMod d_mpep_pep_idx d_mpep_pep_mod_idx d_mpep_rank d_mpep_ith_valid d_mpep_mod_ma_count_sum d_mpep_mod_ma_count_sum_scan d_mod_ma_count_sum d_pep_idx d_pep_mod_idx d_pep_valid_idx d_pep_num_mpep num_pep num_mpep >>= \ len_unrank ->

  CUDA.allocaArray (len_unrank) $ \d_mpep_unrank -> do
    CUDA.genModCands d_mpep_unrank d_mod_ma_count d_mpep_ith_valid d_mpep_rank d_mpep_mod_ma_count_sum_scan d_pep_mod_idx d_pep_ma_count d_pep_valid_idx d_pep_ma_num_comb_scan num_mpep num_ma 

    action (num_mpep, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank)
      --(t,_) <- bracketTime $ CUDA.genModCands d_mpep_idx d_mpep_rank d_mpep_unrank num_mpep d_pep_idx d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan num_pep d_pep_ma_count d_mod_ma_count mod_num_ma
      --when (verbose cp) $ hPutStrLn stderr ("genModCands Elapsed time: " ++ showTime t)
      --action (num_mpep, d_mpep_idx, d_mpep_unrank)

mkModSpecXCorr :: DeviceSeqDB
               -> DeviceModInfo
               -> ModCandidates
               -> Float
               -> Int
               -> (IonSeries -> IO b)
               -> IO b
mkModSpecXCorr db dmi (num_mpep, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank) chrg len action =
  let (num_ma, d_ma, d_ma_mass) = devModAcids dmi 
      (_, d_mod_ma_count, _, d_mod_delta) = devModCombs dmi 
  in 
  CUDA.allocaArray n $ \d_mspec -> do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_mspec
    CUDA.memset  d_mspec bytes 0
    CUDA.addModIons d_mspec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_mpep_pep_idx d_mpep_pep_mod_idx d_mpep_unrank d_mpep_mod_ma_count_sum_scan len_unrank num_mpep d_mod_ma_count d_mod_delta d_ma d_ma_mass num_ma len ch
    --(t,_) <- bracketTime $ CUDA.addModIons d_mspec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_mpep_idx d_mpep_unrank num_mpep d_mod_ma d_mod_ma_count d_mod_ma_mass mod_num_ma ch len
    --hPutStrLn stderr ("addModIons Elapsed time: " ++ showTime t)

    action d_mspec
  where
    n = len * num_mpep
    ch = round $ max 1 (chrg - 1)

--
-- Generate a theoretical spectral representation for each of the specified
-- candidates. This generates spectral peaks for all of the A-, B- and Y-ions of
-- the given sequences, retaining only the most intense peak in each bin.
--
-- On the device, this is stored as a dense matrix, each row corresponding to a
-- single sequence.
--
mkSpecXCorr :: DeviceSeqDB -> CandidatesByMass -> Float -> Int -> (IonSeries -> IO b) -> IO b
mkSpecXCorr db (nIdx, d_idx) chrg len action =
  CUDA.allocaArray n $ \d_spec -> do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_spec

    CUDA.memset  d_spec bytes 0
    CUDA.addIons d_spec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_idx nIdx ch len

    action d_spec
  where
    n  = len * nIdx
    ch = round $ max 1 (chrg - 1)


--
-- Score each candidate sequence against the observed intensity spectra,
-- returning the most relevant results.
--
sequestXC :: ConfigParams -> CandidatesByMass -> Spectrum -> IonSeries -> IO [(Float,Int)]
sequestXC cp (nIdx,d_idx) expr d_thry = let n' = max (numMatches cp) (numMatchesDetail cp) in
  CUDA.withVector  expr $ \d_expr  ->
  CUDA.allocaArray nIdx $ \d_score -> do
    --when (verbose cp) $ hPutStrLn stderr ("Matched peptides: " ++ show nIdx)

    -- There may be no candidates as a result of bad database search parameters,
    -- or if something unexpected happened (out of memory)
    --
    if nIdx == 0 then return [] else do



    -- Score and rank each candidate sequence
    --
    CUDA.mvm   d_score d_thry d_expr nIdx (G.length expr)
    CUDA.rsort d_score d_idx nIdx

    -- Retrieve the most relevant matches
    --
    let n = min n' nIdx
    sc <- CUDA.peekListArray n d_score
    ix <- CUDA.peekListArray n d_idx

    return $ zipWith (\s i -> (s/10000,fromIntegral i)) sc ix

-- Modified candidates version
sequestXCMod :: ConfigParams -> HostModInfo -> DeviceModInfo -> ModCandidates -> Spectrum -> IonSeries -> IO [(Float,Int,PepMod,[Int])]
sequestXCMod cp hmi dmi (num_mpep, d_mpep_pep_idx, d_mpep_pep_mod_idx, d_mpep_unrank, d_mpep_mod_ma_count_sum_scan, len_unrank) expr d_thry = 
  let (num_ma, d_ma, d_ma_mass) = devModAcids dmi 
      (num_mod, d_mod_ma_count, d_mod_ma_count_sum, d_mod_delta) = devModCombs dmi 
      (_, ma, ma_mass) = modAcids hmi
      (_, mod_ma_count, mod_ma_count_sum, _) = modCombs hmi

      n' = max (numMatches cp) (numMatchesDetail cp) 
      h_mpep_idx = G.generate num_mpep (\i -> cIntConv i) :: (U.Vector Word32) 
  in
    -- There may be no candidates as a result of bad database search parameters,
    -- or if something unexpected happened (out of memory)
    --
  if num_mpep == 0 then return [] else 
  CUDA.withVector  expr             $ \d_expr -> 
  CUDA.withVector  h_mpep_idx       $ \d_mpep_idx -> 

  CUDA.allocaArray num_mpep         $ \d_score -> do
    --when (verbose cp) $ hPutStrLn stderr ("Candidate modified peptides: " ++ show num_mpep)


    -- Score and rank each candidate sequence
    --
    CUDA.mvm   d_score d_thry d_expr num_mpep (G.length expr)
    CUDA.rsort d_score d_mpep_idx num_mpep

    -- Retrieve the most relevant matches
    --
    let n = min n' num_mpep
    score             <- CUDA.peekListArray n d_score
    mpep_idx'         <- CUDA.peekListArray n d_mpep_idx
    mpep_pep_idx'     <- CUDA.peekListArray num_mpep d_mpep_pep_idx
    mpep_pep_mod_idx' <- CUDA.peekListArray num_mpep d_mpep_pep_mod_idx
    mpep_unrank'      <- CUDA.peekListArray len_unrank d_mpep_unrank
    mpep_mod_ma_count_sum_scan' <- CUDA.peekListArray num_mpep d_mpep_mod_ma_count_sum_scan

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
        
        --pmod            = zipWith3 (\a b c -> (w2c a, fromIntegral b, c)) ma ma_count ma_mass
        getUnrank i mi = U.toList $ sublist (mpep_mod_ma_count_sum_scan !! i) 
                                    (mod_ma_count_sum U.! mi) $ U.fromList mpep_unrank

    return $ zipWith pack score mpep_idx
    where
        
    sublist beg num ls = U.drop beg $ U.take (beg + num) ls
