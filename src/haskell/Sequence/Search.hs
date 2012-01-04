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

import Debug.Trace

type CandidatesByMass = (Int, DevicePtr Word32)
type CandidatesByMod  = (Int, DevicePtr Word32, DevicePtr Word32)
type ModCandidates = (Int, DevicePtr Word32, DevicePtr Word32)
type IonSeries  = DevicePtr Word32

type PepModDevice = (Int, DevicePtr Word8, DevicePtr Word8, Int, DevicePtr Float)
--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> IO MatchCollection
searchForMatches cp sdb ddb ms2 = do
  matches <- searchWithoutMods cp sdb ddb ms2 
  matchesMods <- searchWithMods cp sdb ddb ms2 
  return $ reverse $ sortBy matchScoreOrder $ matches ++ matchesMods

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
searchWithMods :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> IO MatchCollection
searchWithMods cp sdb ddb ms2 = traceShow modCombs $ do
  matches <- mapM (searchAModComb cp sdb ddb ms2) modCombs
  
  
  traceShow ("max", max_mass) $ traceShow ("min", min_mass) $ return $ concat matches
  where
    comparedbFrag (res1,_,_) (res2,_,_) = compare res1 res2
    (max_mass,_,_)  = U.maximumBy comparedbFrag $ dbFrag sdb 
    (min_mass,_,_)  = U.minimumBy comparedbFrag $ dbFrag sdb 
    max_ma          = 20  -- @TODO
    ma = map c2w ['A','C']
    ma_mass       = [15.15,-75.75] -- @TODO get this from config params
    --ma_count = [2,1]
    --modCombs = [(ma, ma_count, 1123.4)]
    --[(['A'],[2],343.1)]
    
    -- @TODO generate using liftM
    modCombs  = filter (\(_,_,m) -> if min_mass <= m && m <= max_mass then True else False) $ mods_raw 
    mods_raw  = map (\ma_cnt -> (ma, (map fromIntegral ma_cnt), mass - (calcDelta ma_cnt) )) ma_counts
        where 
        calcDelta ma_cnt = sum (zipWith (\c m -> m * fromIntegral c) ma_cnt ma_mass)
    ma_counts  = filter (\c -> if sum c < max_ma then True else False) ma_counts'
    ma_counts' = tail (addModDimensions ((length ma) - 1) [ replicate (length ma) 0 ])

    addModDimensions :: Int -> [[Int]] -> [[Int]]
    addModDimensions dim xs = 
        if dim >= 0 then expand $ addModDimensions (dim-1) xs else xs
        where
        expand (c:cs) = [c] ++ map (\cnt -> replace dim cnt c) [1..max_ma] ++ expand cs
        expand _      = []

    replace :: (Show a) => Int -> a -> [a] -> [a]
    replace idx val list = x ++ val : ys
        where
        (x,_:ys) = splitAt idx list

    mass          = (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 - 1) * massH) - (massH + massH2O)


    
searchAModComb :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> ([Word8], [Word8], Float) -> IO MatchCollection
searchAModComb cp sdb ddb ms2 (ma, ma_count, mod_mass) = 
  filterCandidateByMass cp ddb mod_mass                                     $ \candidatesByMass ->
  CUDA.withVector (U.fromList ma)       $ \d_mod_ma       ->   -- modifable acid
  CUDA.withVector (U.fromList ma_count) $ \d_mod_ma_count ->   -- number of the acid to modify
  CUDA.withVector (U.fromList ma_mass)  $ \d_mod_ma_mass  ->   -- 
  let d_mods = (mod_num_ma, d_mod_ma, d_mod_ma_count, sum_ma_count, d_mod_ma_mass) in
  filterCandidateByModability cp ddb candidatesByMass d_mods                    $ \candidatesByMassAndMod ->
  genModCandidates cp ddb candidatesByMassAndMod d_mods                         $ \modifiedCandidates ->
  mkModSpecXCorr ddb d_mods modifiedCandidates (ms2charge ms2) (G.length spec)  $ \mspecThry   -> -- do
  --mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec)              $ \specThry   -> -- @TODO delete later
  mapMaybe finish `fmap` sequestXCMod cp modifiedCandidates sum_ma_count spec mspecThry -- @TODO
  --return []
  where
    ma_mass         = [15.15,-75.75] -- @TODO get this from config params
    sum_ma_count    = fromIntegral $ sum ma_count
    mod_num_ma      = length ma
    spec            = sequestXCorr cp ms2
    peaks           = extractPeaks spec

    finish (sx,i,u) = liftM (\f -> Match (modifyFragment pmod u f) sx (sp f) pmod u) (lookup sdb i)
    sp              = matchIonSequence cp (ms2charge ms2) peaks
    pmod            = zipWith3 (\a b c -> (w2c a, fromIntegral b, c)) ma ma_count ma_mass
    

--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found.
--
filterCandidateByMass :: ConfigParams -> DeviceSeqDB -> Float -> (CandidatesByMass -> IO b) -> IO b
filterCandidateByMass cp db mass action =
  CUDA.allocaArray np $ \d_idx ->
  CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (n,d_idx)
  where
    np    = numFragments db
    delta = massTolerance cp

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
filterCandidateByModability :: 
                       ConfigParams 
                    -> DeviceSeqDB 
                    -> CandidatesByMass     --- ^ subset of peptides, [idx to pep], number of pep
                    -> PepModDevice               
                    -> (CandidatesByMod -> IO b)
                    -> IO b
filterCandidateByModability cp db (sub_nIdx, d_sub_idx) (mod_num_ma, d_mod_ma, d_mod_ma_count, _, _) action =
  CUDA.allocaArray sub_nIdx             $ \d_idx      ->     -- filtered results to be returned here
  CUDA.allocaArray (sub_nIdx*mod_num_ma)      $ \d_pep_ma_count -> -- peptide ma counts
  CUDA.findModablePeptides d_idx d_pep_ma_count (devIons db) (devTerminals db) d_sub_idx sub_nIdx d_mod_ma d_mod_ma_count mod_num_ma >>= \n -> do
    --when (verbose cp) $ hPutStrLn stderr ("filter by modability:\n" ++ 
                                          --"num searched: " ++ show sub_nIdx ++ "\n" ++
                                          --"num modable: " ++ show n ++ "\n" ++
                                          --"pep_ma_count not used: " ++ show (mod_num_ma*(sub_nIdx - n)) ++ "\n") 
    action (n, d_idx, d_pep_ma_count)
  --CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (d_idx,n)


--
-- Given some unmodified peptides, a modification and
-- count of each modifiable aa in each pep, generate modified peptides
-- 
genModCandidates :: ConfigParams
                 -> DeviceSeqDB
                 -> CandidatesByMod
                 -> PepModDevice
                 -> (ModCandidates -> IO b)
                 -> IO b
genModCandidates cp ddb (nPep, d_pep_idx, d_pep_ma_count) (mod_num_ma, d_mod_ma, d_mod_ma_count, sum_ma_count, _) action =
  -- calc number of modified peps generated from each pep
  CUDA.allocaArray nPep $ \d_pep_num_mpep -> 
  CUDA.allocaArray (nPep*mod_num_ma) $ \d_pep_ma_num_comb -> 
  CUDA.allocaArray (nPep*mod_num_ma) $ \d_pep_ma_num_comb_scan -> 
  CUDA.calcTotalModCands d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan nPep d_pep_ma_count d_mod_ma_count mod_num_ma >>= \total -> 
  --action total
  CUDA.allocaArray total                $ \d_mpep_rank -> 
  CUDA.allocaArray (total*sum_ma_count) $ \d_mpep_unrank -> 
  CUDA.allocaArray total                $ \d_mpep_idx -> do
      CUDA.genModCands d_mpep_idx d_mpep_rank d_mpep_unrank total d_pep_idx d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan nPep d_pep_ma_count d_mod_ma_count mod_num_ma
      action (total, d_mpep_idx, d_mpep_unrank)

mkModSpecXCorr :: DeviceSeqDB
               -> PepModDevice
               -> ModCandidates
               -> Float
               -> Int
               -> (IonSeries -> IO b)
               -> IO b
mkModSpecXCorr db (mod_num_ma, d_mod_ma, d_mod_ma_count,_ , d_mod_ma_mass) (total, d_mpep_idx, d_mpep_unrank) chrg len action =
  CUDA.allocaArray n $ \d_mspec -> do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_mspec
    CUDA.memset  d_mspec bytes 0
    CUDA.addModIons d_mspec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_mpep_idx d_mpep_unrank total d_mod_ma d_mod_ma_count d_mod_ma_mass mod_num_ma ch len

    action d_mspec
  where
    n = len * total
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
    when (verbose cp) $ hPutStrLn stderr ("Matched peptides: " ++ show nIdx)

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
sequestXCMod :: ConfigParams -> ModCandidates -> Int -> Spectrum -> IonSeries -> IO [(Float,Int,[Int])]
sequestXCMod cp (nMCands, d_mpep_idx, d_mpep_unrank) sum_ma_count expr d_thry = let 
  n' = max (numMatches cp) (numMatchesDetail cp) 
  h_mpep_idx_idx = G.generate nMCands (\i -> cIntConv i) :: (U.Vector Word32) in
    -- There may be no candidates as a result of bad database search parameters,
    -- or if something unexpected happened (out of memory)
    --
  if nMCands == 0 then return [] else 
  CUDA.withVector  expr             $ \d_expr  -> 
  CUDA.withVector  h_mpep_idx_idx   $ \d_mpep_idx_idx -> 

  CUDA.allocaArray nMCands          $ \d_score -> do
    when (verbose cp) $ hPutStrLn stderr ("Candidate modified peptides: " ++ show nMCands)


    -- Score and rank each candidate sequence
    --
    CUDA.mvm   d_score d_thry d_expr nMCands (G.length expr)
    CUDA.rsort d_score d_mpep_idx_idx nMCands

    -- Retrieve the most relevant matches
    --
    let n = min n' nMCands
    score         <- CUDA.peekListArray n d_score
    mpep_idx_idx' <- CUDA.peekListArray n d_mpep_idx_idx
    mpep_idx'     <- CUDA.peekListArray nMCands d_mpep_idx
    mpep_unrank'  <- CUDA.peekListArray (nMCands*sum_ma_count) d_mpep_unrank

    let mpep_idx_idx = map fromIntegral mpep_idx_idx'
        mpep_idx     = map fromIntegral mpep_idx'
        mpep_unrank  = map fromIntegral mpep_unrank'


    return $ zipWith (\s i -> (s/10000, (mpep_idx !! i), (sublist (i*sum_ma_count) (i*sum_ma_count + sum_ma_count) mpep_unrank))) score mpep_idx_idx
    where
    sublist beg end ls = drop beg $ take end ls
