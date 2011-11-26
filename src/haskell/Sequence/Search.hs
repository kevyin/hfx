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


type Candidates = (DevicePtr Word32, Int)
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
  filterCandidateByModability cp ddb candidatesByMass pepMod             $ \candidatesByMassModability ->
  mkSpecXCorr ddb candidatesByMass (ms2charge ms2) (G.length spec) $ \specThry   ->
  mapMaybe finish `fmap` sequestXC cp candidatesByMass spec specThry
  where
    pepMod        = ((map c2w ['M','C']), [1,2])
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
filterCandidateByMass :: ConfigParams -> DeviceSeqDB -> Float -> (Candidates -> IO b) -> IO b
filterCandidateByMass cp db mass action =
  CUDA.allocaArray np $ \d_idx ->
  CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (d_idx,n)
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
-- @TODO also record aa counts for each peptide
--
filterCandidateByModability :: 
                       ConfigParams 
                    -> DeviceSeqDB 
                    -> (DevicePtr Word32, Int)  --- ^ subset of peptides, [idx to pep], number of pep
                    -> ([Word8],[Word8])        --- ^ peptide mod eg (['A','R'],[2,3])
                    -> ((DevicePtr Word32, Int) -> IO b)
                    -> IO b
filterCandidateByModability cp db (d_sub_idx, sub_nIdx) (ma, ma_count) action =
  CUDA.withVector (U.fromList ma)        $ \d_ma       ->    -- modifable acid
  CUDA.withVector (U.fromList ma_count)  $ \d_ma_count ->    -- number of the acid to modify
  CUDA.allocaArray sub_nIdx              $ \d_idx      ->    -- filtered results to be returned here
  CUDA.findModablePeptides d_idx (devIons db) (devTerminals db) d_sub_idx sub_nIdx d_ma d_ma_count >>= \n -> action (d_sub_idx, sub_nIdx)
  --CUDA.findIndicesInRange (devResiduals db) d_idx np (mass-delta) (mass+delta) >>= \n -> action (d_idx,n)


--
-- Generate a theoretical spectral representation for each of the specified
-- candidates. This generates spectral peaks for all of the A-, B- and Y-ions of
-- the given sequences, retaining only the most intense peak in each bin.
--
-- On the device, this is stored as a dense matrix, each row corresponding to a
-- single sequence.
--
mkSpecXCorr :: DeviceSeqDB -> Candidates -> Float -> Int -> (IonSeries -> IO b) -> IO b
mkSpecXCorr db (d_idx, nIdx) chrg len action =
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
sequestXC :: ConfigParams -> Candidates -> Spectrum -> IonSeries -> IO [(Float,Int)]
sequestXC cp (d_idx,nIdx) expr d_thry = let n = max (numMatches cp) (numMatchesDetail cp) in
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
    sc <- CUDA.peekListArray n d_score
    ix <- CUDA.peekListArray n d_idx

    return $ zipWith (\s i -> (s/10000,fromIntegral i)) sc ix

