{-# LANGUAGE TupleSections, ScopedTypeVariables, PatternGuards #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Main
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- A simple test file for the peptide sequence matching algorithm
--
--------------------------------------------------------------------------------

module Main where

--
-- Custom libraries
--
import Config
import Sequence
import Spectrum
import Util.Time
import Util.Show
import Util.PrettyPrint

--
-- System libraries
--
import Data.List
import Data.Maybe
import Control.Monad
import Control.Exception
import System.Environment
import System.FilePath
import System.IO
import Prelude                          hiding ( lookup, catch )

import qualified Data.Vector.Generic    as G
import qualified Foreign.CUDA           as CUDA


--------------------------------------------------------------------------------
-- Program Defaults
--------------------------------------------------------------------------------

--defaultConfigFile :: FilePath
--defaultConfigFile =  "hfx.params"


--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------

main :: IO ()
main = do
  (cp,dta) <- sequestConfig =<< getArgs
  let fp   = fromMaybe (error "Protein database not specified") (databasePath cp)

  when (verbose cp && not (useCPU cp)) $ do
    dev   <- CUDA.get
    props <- CUDA.props dev
    hPutStrLn stderr $ "Device " ++ shows dev ": "
                                 ++ CUDA.deviceName props ++ ", "
                                 ++ "compute v" ++ shows (CUDA.computeCapability props) ", "
                                 ++ "global memory: " ++ showFFloatSIBase 1024 (fromIntegral $ CUDA.totalGlobalMem props :: Double) "B, "
                                 ++ "core clock: "    ++ showFFloatSI (fromIntegral $ 1000 * CUDA.clockRate props :: Double) "Hz"
    when (verbose cp) $ do
      hPutStrLn stderr $ "Database: " ++ fp

  --
  -- Load the proteins from file, marshal to the device, and then get to work!
  --
  when (verbose cp) $ hPutStrLn stderr ("Loading Database ...\n" )
  --(cp',dbs) <- loadDatabase cp fp (splitDB cp)
  (t,(cp',dbs)) <- bracketTime $ loadDatabase cp fp (splitDB cp)
  when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)

  when (verbose cp) $ hPutStrLn stderr ("Searching ...\n" )
  (t2,allMatches') <- bracketTime $ forM dbs $ \db -> do
    hmi <- makeModInfo cp' db
    withDevModInfo hmi $ \dmi -> do
      withDeviceDB cp' db $ \ddb -> forM dta $ \f -> do
        matches <- search cp' db hmi dmi ddb f
        when (showMatchesPerScan cp) $ printScanResults cp' f matches
        return matches
  when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t2)
          

  let sortXC = sortBy (\(_,_,a) (_,_,b) -> compare (scoreXC b) (scoreXC a))
      matchesRaw = concat $ concat $ concat allMatches'
      --matchesByScan = map sortXC $ groupBy (\(_,ms,_) (_,ms',_) -> ms == ms') matchesRaw
      --matchesByFile = map sortXC $ groupBy (\(f,_,_) (f',_,_) -> f == f') matchesRaw
      n = maximum $ [(numMatches cp'), (numMatchesDetail cp'), (numMatchesIon cp')]
      allMatchesN = take n $! sortXC $ matchesRaw

  printAllResults cp' allMatchesN
    


{-# INLINE loadDatabase #-}
loadDatabase :: ConfigParams -> FilePath -> Int -> IO (ConfigParams, [SequenceDB])
loadDatabase cp fp split = do
  (cp',dbs) <- case takeExtensions fp of
    ext | ".fasta" `isPrefixOf` ext -> (cp,) `liftM` makeSeqDB cp fp split
    ext | ".index" `isPrefixOf` ext -> readIndex cp fp >>= \(cp',sdb) -> return (cp',[sdb]) 
    _                               -> error ("Unsupported database type: " ++ show fp)
  return (cp',dbs)


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> HostModInfo -> DeviceModInfo ->  DeviceSeqDB -> FilePath -> IO [[(FilePath,MS2Data,Match)]]
search cp db hmi dmi dev fp =
  readMS2Data fp >>= \r -> case r of
    Left  s -> hPutStrLn stderr s >>= \_ -> return []
    Right d -> do 
      when (verbose cp) $ do
        --hPutStrLn stderr $ "Database: " ++ fp
        hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
        hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
        hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

      (t,allMatches) <- bracketTime $ forM d $ \ms2 -> do
        matches <- searchForMatches cp db dev hmi dmi ms2
        return $ map (\m -> (fp, ms2, m)) matches
        `catch`
        \(e :: SomeException) -> do 
                hPutStrLn stderr $ unlines [ "scan:   " ++ (show $ ms2info ms2), "reason: " ++ show e ]
                return []

      when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)
      return allMatches

      `catch`
      \(e :: SomeException) -> do 
          hPutStrLn stderr $ unlines
              [ "file:   " ++ fp
              , "reason: " ++ show e ]
          return []
