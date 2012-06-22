{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Index
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Reading and writing digested sequence database
--
--------------------------------------------------------------------------------

module Sequence.Index where

import Mass
import Config
import Sequence.Fragment

import Data.Ix
import Data.Char
import Data.Maybe
import Data.Binary
import Data.Binary.Get
import Codec.Compression.GZip
import System.IO

import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Vector.Generic        as G

import Debug.Trace


--
-- Output the digested sequence database together with the relevant digestion
-- rules used to generate it, in a form suitable for the configuration parser.
--
-- Annoyingly, binary instances for generic vector match with String, so pack
-- into a ByteString first. Probably not too bad, as 'encode' would have done
-- this anyway.
--
writeIndex :: Handle -> ConfigParams -> SequenceDB -> IO ()
writeIndex hdl cp db = encodeFile (fromJust $ outputPath cp) db
{-
writeIndex :: Handle -> ConfigParams -> SequenceDB -> IO ()
writeIndex hdl cp db = L.hPut hdl header >> L.hPut hdl content
  where
    content = compress (encode db)
    header  = encode . L.pack . ('\n':) . unlines $
      [ "database         = " ++ (fromJust (databasePath cp))
      , "missed-cleavages = " ++ (show (missedCleavages cp))
      , "digestion-rule   = " ++ (takeWhile isDigit . dropWhile isSpace . snd $ digestionRule cp)
      , "min-peptide-mass = " ++ (show (minPeptideMass cp))
      , "max-peptide-mass = " ++ (show (maxPeptideMass cp))
      , "aa-mass-type     = " ++ (if aaMassTypeMono cp then "Mono" else "Average")
      ]
      ++ G.ifoldr' aaMod [] (aaMassTable cp)

    aaMod i m mods =
      let aa  = chr (ord 'A' + i)
          det = m - if aaMassTypeMono cp then getAAMassMono aa else getAAMassAvg aa
      in
      if det == 0 then mods
                  else ("add_" ++ [aa] ++ " = " ++ show det) : mods

-}


--
-- Read an indexed sequence database. Update the configuration options with
-- those used to generate the database.
--
readIndex :: ConfigParams -> FilePath -> IO (ConfigParams, SequenceDB)
readIndex cp fp = do
  f   <- L.readFile fp
  --let (opt,f',_) = runGetState get f 0
      --db         = decode (decompress f')
      --table      = G.replicate (rangeSize ('A','Z')) 0

  putTraceMsg "decoding file"
  db         <- decodeFile fp
  putTraceMsg "decoded file"
  putTraceMsg $ " dbion length " ++ show (G.length $ dbIon db)
  --(,db) `fmap` readConfig (L.unpack opt) fp (cp {aaMassTable = table, databasePath = Just fp})

  return (cp,db)

