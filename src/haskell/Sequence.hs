{-# LANGUAGE BangPatterns  #-}
{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence
-- Copyright : (c) 2010 Trevor L. McDonell
-- License   : BSD
--
-- Reading and processing sequence related stuff
--
--------------------------------------------------------------------------------

module Sequence
  (
    Protein,
    readFasta, countSeqs,
    numIons, ionMasses, peptides, digest,

    Fragment(..), fraglabel
  )
  where

import Mass
import Config

import Data.Word
import Data.List
import Control.Applicative

import qualified Bio.Sequence                      as F
import qualified Bio.Sequence.Fasta                as F
import qualified Data.ByteString.Lazy.Char8        as L
import qualified Codec.Compression.GZip            as GZip

import qualified Data.Vector.Unboxed               as U
import qualified Data.Vector.Generic               as G
import qualified Data.Vector.Fusion.Stream         as Stream
import qualified Data.Vector.Fusion.Stream.Monadic as S
import qualified Data.Vector.Fusion.Stream.Size    as S
import Data.Vector.Fusion.Util


--------------------------------------------------------------------------------
-- Reading FASTA Files
--------------------------------------------------------------------------------

type Protein = F.Sequence F.Amino

--
-- Lazily read sequences from a FASTA-formatted file. This is identical to the
-- code of the bio package, except the sequence type is cast on output and GZip
-- compressed files are deflated as necessary.
--
readFasta :: FilePath -> IO [Protein]
{-# INLINE readFasta #-}
readFasta fp = map F.castToAmino . F.mkSeqs . L.lines . prepare <$> L.readFile fp
  where
    prepare = if ".gz" `isSuffixOf` fp then GZip.decompress
                                       else id

--
-- Count the number of sequences in a file. Each sequence consist of a header
-- and a set of lines containing the sequence data. GZip compressed files are
-- inflated if necessary, which is the only addition over the bio package
-- function of the same name.
--
countSeqs :: FilePath -> IO Int
{-# INLINE countSeqs #-}
countSeqs fp = length . headers . prepare <$> L.readFile fp
  where
    headers = filter (('>' ==) . L.head) . filter (not . L.null) . L.lines
    prepare = if ".gz" `isSuffixOf` fp then GZip.decompress
                                       else id


--------------------------------------------------------------------------------
-- Fragments
--------------------------------------------------------------------------------

--
-- Sequence fragments (typically to store results for pretty printing)
--
data Fragment = Fragment
  {
    fragmass   :: Float,
    fragheader :: L.ByteString,
    fragdata   :: L.ByteString
  }
  deriving (Eq, Show)

fraglabel :: Fragment -> L.ByteString
fraglabel = head . L.words . fragheader


--------------------------------------------------------------------------------
-- Digestion
--------------------------------------------------------------------------------

fst3 :: (a,b,c) -> a
fst3 (a,_,_) = a

snd3 :: (a,b,c) -> b
snd3 (_,b,_) = b

thd3 :: (a,b,c) -> c
thd3 (_,_,c) = c


--
-- Determine the total number of ions represented by the protein sequences
--
numIons :: [Protein] -> Int
{-# INLINE numIons #-}
numIons = fromIntegral . foldl (+) 0 . map (L.length . F.seqdata)

--
-- Translate the amino acid character codes of the database into the
-- corresponding ion masses in a single, flattened vector.
--
ionMasses :: ConfigParams -> Int -> [Protein] -> U.Vector Float
{-# INLINE ionMasses #-}
ionMasses cp n db = G.unstream (ionMassesS cp n db)

ionMassesS :: ConfigParams -> Int -> [Protein] -> Stream.Stream Float
{-# INLINE [1] ionMassesS #-}
ionMassesS cp n db = S.Stream (return . step) (map F.seqdata db) (S.Exact (delay_inline max n 0))
  where
    {-# INLINE [0] step #-}
    step []     = S.Done
    step (s:ss) = case L.uncons s of
                    Nothing    -> S.Skip ss
                    Just (c,t) -> S.Yield (getAAMass cp c) (t:ss)


--
-- Scan the input sequences to generate a flattened vector of the peptide
-- fragments found using the current digestion rules. The output indices are
-- global offsets, not into the individual parent sequence. This is suitable for
-- indexing into the similarly globally concatenated 'ionMass' output.
--
peptides :: ConfigParams -> [Protein] -> U.Vector (Float,Word32,Word32)
peptides cp db
  = U.fromList . concat
  . zipWith global offset
  $ map (digest cp) db
  where
    offset   = scanl (+) 0 . map (fromIntegral . L.length . F.seqdata) $ db
    global o = map (\(r,c,n) -> (r,c+o,n+o))


--
-- Digest a single protein sequence with the given enzyme and cleavage rules,
-- removing any fragments outside the broad mass range of interest.
--
-- The mass of the peptide is the sum of the amino acid residue masses plus the
-- mass of the water molecule released in forming the peptide bond (plus one;
-- from Eq. 1 of Eng.[1])
--
digest :: ConfigParams -> Protein -> [(Float,Word32,Word32)]
digest cp = filter (inrange . fst3) . splice cp . fragment cp
  where
    inrange m = let m' = m + massH2O + massH
                in  minPeptideMass cp <= m' && m' <= maxPeptideMass cp


--
-- Split a protein sequence with the given digestion rule. This outputs the
-- total residual mass of the fragment, together with the indices of the C- and
-- N- terminals from parent sequence.
--
fragment :: ConfigParams -> Protein -> [(Float,Word32,Word32)]
fragment cp = unfoldr step . (0,) . F.seqdata
  where
    rule = fst . digestionRule $ cp

    -- add strictness here?
    --
    chunk c n s = let (h,t) = L.splitAt (n+1) s
                      mass  = L.foldl' (\m a -> m + getAAMass cp a) 0 h
                  in  Just ((mass, fromIntegral c, fromIntegral (c+n)), (c+n+1, t))

    step (!o,!s) = case L.findIndex rule s of
      Nothing -> if L.null s then Nothing
                             else chunk o (L.length s - 1) s
      Just !i -> chunk o i s


--
-- Generate additional sequences from missed cleavages of sequential, adjacent
-- peptide fragments.
--
splice :: ConfigParams -> [(Float,Word32,Word32)] -> [(Float,Word32,Word32)]
splice cp = loop
  where
    loop []     = []
    loop (p:ps) = scanl join p (take n ps) ++ loop ps

    n        = missedCleavages cp
    join a b = (fst3 a + fst3 b, snd3 a, thd3 b)

