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

import Data.Word
import Data.Maybe
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
type PepMod = (Int, DevicePtr Word8, DevicePtr Word8)
type IonSeries  = DevicePtr Word32

--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> IO MatchCollection
searchForMatches cp sdb ddb ms2 =
  filterCandidateByMass cp ddb mass                                   $ \candidatesByMass ->
  CUDA.withVector (U.fromList ma)       $ \d_mod_ma       ->     -- modifable acid
  CUDA.withVector (U.fromList ma_count) $ \d_mod_ma_count ->     -- number of the acid to modify
  filterCandidateByModability cp ddb candidatesByMass (mod_num_ma, d_mod_ma, d_mod_ma_count) $ \candidatesByMassMod ->
  genModCandidates cp ddb candidatesByMassMod (mod_num_ma, d_mod_ma, d_mod_ma_count) $ \modifiedCandidates ->
  mkModSpecXCorr ddb (mod_num_ma, d_mod_ma, d_mod_ma_count) modifiedCandidates (ms2charge ms2) (G.length spec) $ \mspecThry   ->
  mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec) $ \specThry   ->
  mapMaybe finish `fmap` sequestXC cp candidatesByMass spec specThry
  where
    ma            = map c2w ['S','A'] 
    ma_count      = [3,1]
    mod_num_ma = length ma
    spec          = sequestXCorr cp ms2
    peaks         = extractPeaks spec

    finish (sx,i) = liftM (\f -> Match f sx (sp f)) (lookup sdb i)
    sp            = matchIonSequence cp (ms2charge ms2) peaks
    mass          = (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 - 1) * massH) - (massH + massH2O)


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
                    -> PepMod               
                    -> (CandidatesByMod -> IO b)
                    -> IO b
filterCandidateByModability cp db (sub_nIdx, d_sub_idx) (mod_num_ma, d_mod_ma, d_mod_ma_count) action =
  CUDA.allocaArray sub_nIdx             $ \d_idx      ->     -- filtered results to be returned here
  CUDA.allocaArray (sub_nIdx*mod_num_ma)      $ \d_pep_ma_count -> -- peptide ma counts
  CUDA.findModablePeptides d_idx d_pep_ma_count (devIons db) (devTerminals db) d_sub_idx sub_nIdx d_mod_ma d_mod_ma_count mod_num_ma >>= \n -> 
    traceShow ("num peptides searched", sub_nIdx) $
    traceShow ("num modable", n) $
    traceShow ("pep_ma_count not used ", mod_num_ma*(sub_nIdx - n)) $ 
    action (n, d_idx, d_pep_ma_count)
  --CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (d_idx,n)


--
-- Given some unmodified peptides, a modification and
-- count of each modifiable aa in each pep, generate modified peptides
-- 
genModCandidates :: ConfigParams
                 -> DeviceSeqDB
                 -> CandidatesByMod
                 -> PepMod
                 -> (ModCandidates -> IO b)
                 -> IO b
genModCandidates cp ddb (nPep, d_pep_idx, d_pep_ma_count) (mod_num_ma, d_mod_ma, d_mod_ma_count) action =
  -- calc number of modified peps generated from each pep
  CUDA.allocaArray nPep $ \d_pep_num_mpep -> 
  CUDA.allocaArray (nPep*mod_num_ma) $ \d_pep_ma_num_comb -> 
  CUDA.allocaArray (nPep*mod_num_ma) $ \d_pep_ma_num_comb_scan -> 
  CUDA.calcTotalModCands d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan nPep d_pep_ma_count d_mod_ma_count mod_num_ma >>= \total -> 
  --action total
  CUDA.allocaArray total $ \d_mpep_mcomb -> do
  CUDA.allocaArray total $ \d_mpep_idx -> do
      CUDA.genModCands d_mpep_idx d_mpep_mcomb total d_pep_idx d_pep_num_mpep d_pep_ma_num_comb d_pep_ma_num_comb_scan nPep d_pep_ma_count d_mod_ma_count mod_num_ma
      action (total, d_mpep_idx, d_mpep_mcomb)

mkModSpecXCorr :: DeviceSeqDB
               -> PepMod
               -> ModCandidates
               -> Float
               -> Int
               -> (IonSeries -> IO b)
               -> IO b
mkModSpecXCorr db (mod_num_ma, d_mod_ma, d_mod_ma_count) (total, d_mpep_idx, d_mpep_mcomb) chrg len action =
  CUDA.allocaArray n $ \d_mspec -> do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_mspec
    CUDA.memset  d_mspec bytes 0
    CUDA.addModIons d_mspec (devResiduals db) (devMassTable db) (devIons db) (devTerminals db) d_mpep_idx d_mpep_mcomb total d_mod_ma d_mod_ma_count mod_num_ma ch len

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
    --CUDA.rsort d_score d_idx nIdx

    -- Retrieve the most relevant matches
    --
    let n = min n' nIdx
    sc <- CUDA.peekListArray n d_score
    ix <- CUDA.peekListArray n d_idx

    return $ zipWith (\s i -> (s/10000,fromIntegral i)) sc ix

