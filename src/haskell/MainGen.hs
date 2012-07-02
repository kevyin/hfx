--------------------------------------------------------------------------------
-- |
-- Module    : MainGen
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Generate a list of peptides that meet the digestion criteria
--
--------------------------------------------------------------------------------

module Main where

import Config
import Sequence.Index
import Sequence.Fragment
import Util.Time
import Util.Show

import Data.Maybe
import Control.Monad
import System.IO
import System.Environment
import System.Exit
import qualified Data.Vector.Generic as G


main :: IO ()
main = do
  args   <- getArgs
  (cp,f) <- sequestConfig args
  let fp =  fromMaybe (error "Protein database not specified") (databasePath cp)
  case outputPath cp of
       Just path -> when (verbose cp) $ hPutStrLn stderr ("Output file: " ++ show path)
       _         -> do hPutStrLn stderr $ "ERROR: No output path was specified\n" 
                       exitFailure

  (t,dbs') <- bracketTime $ makeSeqDB cp fp (splitDB cp)
  when (verbose cp) $ hPutStrLn stderr ("Reading time: " ++ showTime t)

  let dbs = zip (if length dbs' > 1 then map show [1..(length dbs')] else [""]) dbs'

  forM_ dbs $ \(suf,db) -> do
      if null f
         then writeIndex stdout cp db
         else withFile ((head f) ++ suf) WriteMode (\h -> writeIndex h cp db)

      hPutStrLn stderr $ "Database: " ++ fp
      hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
      hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
      hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

