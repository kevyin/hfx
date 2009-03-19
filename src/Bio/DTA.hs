{-
 - Parse a DTA LC-MS/MS results file.
 -
 - The file format is very simple. The first line contains the singly protonated
 - peptide mass (MH+) and the peptide charge state as a pair of space separated
 - values. Subsequent lines contain space separated pairs of fragment ion m/z
 - ratio and intensity values. Note that the precursor peptide mass is
 - independent of the charge state.
 -
 - The filename usually used to identify the dataset, and each file contains
 - only a single MS/MS sample set.
 -}

module Bio.DTA
    (
      Spectrum(..),                     -- Data structures
      getPrecur, getParent, getSpec,    -- Accessor functions
      readDTA                           -- File formats
    ) where

import Numeric
import Control.Monad (liftM2)
import Text.ParserCombinators.Parsec


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- The mass spectroscopy data generated from a protein identification
-- experiment. It consists of:
--
data Spectrum = Spec
        Double                  -- The singly protonated peptide mass
        Double                  -- Peptide charge state
        [(Double, Double)]      -- The actual mass/charge ratio intensity measurements
    deriving (Eq, Show)


--------------------------------------------------------------------------------
-- DTA File Parser/Lexer
--------------------------------------------------------------------------------

--
-- The DTA file contains at least one line, each of which is terminated by an
-- end-of-line character (eol)
--
dtaFile :: GenParser Char st [(Double, Double)]
dtaFile =  endBy line eol

-- 
-- Each line contains exactly two data values, separated by white space. These
-- are returned as a pair of (mass/charge ratio, intensity) values. Detecting
-- signed values isn't really necessary, but is done for completeness.
--
line :: GenParser Char st (Double, Double)
line =  liftM2 (,) fval fval
    where fval = (fst . head . readSigned readFloat) `fmap` value

--
-- Each value is a floating point number. Discard any leading white space
-- encountered.
--
value :: GenParser Char st [Char]
value =  skipMany (oneOf " \t") >> getValue
    where getValue =  many1 (oneOf (['0'..'9']++"-."))
                  <?> "floating-point number"

--
-- The end of line character. Different operating systems use different
-- characters to mark the end-of-line, so just look for all combinations
--
eol :: GenParser Char st String
eol =  try (string "\n\r")
   <|> try (string "\r\n")
   <|> string "\r"
   <|> string "\n"
   <?> "end of line"


--------------------------------------------------------------------------------
-- Spectrum
--------------------------------------------------------------------------------

--
-- Encase the values read from the DTA file into a Spectrum data structure
--
mkSpec                  :: [(Double, Double)] -> Either String Spectrum
mkSpec []               =  Left "Error: empty spectrum"
mkSpec ((m,c):ss)
    | truncate' c /= c  =  Left "Error: invalid peptide charge state\nexpecting integer"
    | otherwise         =  Right (Spec m c ss)
    where
        truncate' = fromInteger . truncate

getPrecur (Spec p c _) = (p + c - 1) / c
getParent (Spec p _ _) = p
getSpec   (Spec _ _ s) = s

--------------------------------------------------------------------------------
-- File I/O
--------------------------------------------------------------------------------

--
-- Read the given file and return either an error of the Spectrum data.
--
readDTA      :: FilePath -> IO (Either String Spectrum)
readDTA name =  do
    dta <- parseFromFile dtaFile name
    case dta of
        Left  e -> return (Left ("Error parsing file: " ++ show e))
        Right s -> return (mkSpec s)
