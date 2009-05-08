{-
 - Utilities for protein sequences and peptide fragments
 -}


module Protein where

import Utils
import Config
import AminoAcid

import Data.Int
import Data.List
import qualified Bio.Sequence as S
import qualified Data.ByteString.Lazy.Char8 as L


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

type ProteinDatabase = [Protein]

--
-- A protein, represented by its amino acid character code sequence
--
data Protein = Protein
    {
        header    :: L.ByteString,      -- Description of the protein
        chain     :: L.ByteString,      -- Amino acid character sequence
        fragments :: [Peptide]          -- Peptide fragments digested from this protein
    }
    deriving (Eq, Show)

--
-- Extract the name and full description of a protein
--
name, description :: Protein -> String
name        = L.unpack . head . L.words . header
description = L.unpack . header

--
-- A subsequence of a protein
--
-- The mass of the peptide is the sum of the amino acid residue masses plus the
-- mass of the water molecule released in forming the peptide bond (plus one;
-- from Eq. 1 of Eng.[1])
--
-- The ladder represents the b-ion series, generated by successive breaks along
-- the spine of the peptide
--
data Peptide = Peptide
    {
--        parent    :: Protein,           -- Protein this fragment derives from
        residual  :: Float,             -- The sum of the residual masses of this peptide
        ladder    :: [Float],           -- Sequence ladder of b-ion series fragments
        terminals :: (Int64, Int64)     -- Location in the parent protein of this peptide
    }
    deriving (Eq, Show)

pmass   :: Peptide -> Float
pmass p =  residual p + (massH2O + 1.0)


--------------------------------------------------------------------------------
-- File Handlers
--------------------------------------------------------------------------------

--
-- Read a protein database from a Fasta formatted file
--
-- Each entry consists of a header (with a prefix of >) followed by a series of
-- lines containing the sequence data.
--
readFasta :: FilePath -> IO ProteinDatabase
readFasta fasta = do
    database <- S.readFasta fasta
    return   $  map (\(S.Seq h d _) -> Protein h d []) database


--------------------------------------------------------------------------------
-- Protein Fragments
--------------------------------------------------------------------------------

--
-- Record a new protein fragment
--
fragment :: ConfigParams -> Protein -> (Int64, Int64) -> Peptide
fragment cp protein indices = Peptide
    {
--        parent    = protein,
        residual  = subFoldBS' (\m a -> m + getAAMass cp a) 0.0 (chain protein) indices,
        ladder    = subScanBS' (\m a -> m + getAAMass cp a) 0.0 (chain protein) indices,
        terminals = indices
    }


--
-- Scan a protein sequence from the database looking for combinations of amino
-- acids, proceeding from the N to the C terminus.
--
-- Fragments are generated between each of the given list of amino acids where
-- cleavage occurs, and attached to the input protein.
--
-- This almost supports special digestion rules...
--
digestProtein :: ConfigParams -> Protein -> Protein
digestProtein cp protein = protein { fragments = seqs }
    where
        seqs      = filter inrange splices
        inrange p = minPeptideMass cp <= pmass p && pmass p <= maxPeptideMass cp

        indices   = L.findIndices (digestRule cp) (chain protein)
        frags     = simpleFragment cp protein (-1:indices)
        splices   = simpleSplice cp frags


--
-- Split a protein at the given amino acid locations. Be sure to include the
-- final fragment, with a dummy cleavage point at the end of the sequence.
--
simpleFragment :: ConfigParams -> Protein -> [Int64] -> [Peptide]
simpleFragment _  _ []     = []
simpleFragment cp p (i:is) = fragment cp p (i+1, n) : simpleFragment cp p is
    where n = if null is
                then L.length (chain p) - 1
                else head is


--
-- Include the possibility of a number of missed cleavages. All combinations of
-- [0..n] sequential peptides are joined. For example, with two missed
-- cleavages:
--      simpleSplice [a,b,c] -> [a,a:b,a:b:c, b,b:c, c]
--
-- The input list must consist of sorted and adjacent breaks of the original
-- protein sequence, but the output is unordered.
--
simpleSplice :: ConfigParams -> [Peptide] -> [Peptide]
simpleSplice _  []         = []
simpleSplice cp pep@(p:ps)
    | null ps              = [p]
    | otherwise            = scanl1 splice (take n pep) ++ simpleSplice cp ps
    where
        n          = missedCleavages cp + 1
        splice a b = Peptide
            {
                residual  = residual a + residual b,
                ladder    = ladder a ++ map (+residual a) (ladder b),
                terminals = (fst (terminals a), snd (terminals b))
            }


--------------------------------------------------------------------------------
-- Pretty printing
--------------------------------------------------------------------------------

--
-- Return a representation of the amino acid sequence of a peptide, including
-- the flanking residuals (if present)
--
slice :: Protein -> Peptide -> String
slice parent peptide = [ca,'.'] ++ body ++ ['.',na]
    where
        aseq  = chain parent
        (c,n) = terminals peptide
        body  = (L.unpack . L.take (n-c+1) . L.drop c) aseq

        l     = L.length aseq - 2
        ca    = if c > 0 then L.index aseq (c-1) else '-'
        na    = if n < l then L.index aseq (n+1) else '-'

