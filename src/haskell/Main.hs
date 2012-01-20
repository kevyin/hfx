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

  --
  -- Load the proteins from file, marshal to the device, and then get to work!
  --
  when (verbose cp) $ hPutStrLn stderr ("Loading Database ...\n" )
  (cp',db) <- loadDatabase cp fp
  hmi <- makeModInfo cp' db
  when (verbose cp) $ hPutStrLn stderr ("Searching ...\n" )
  withDevModInfo hmi $ \dmi ->
    withDeviceDB cp' db $ forM_ dta . search cp' db hmi dmi


{-# INLINE loadDatabase #-}
loadDatabase :: ConfigParams -> FilePath -> IO (ConfigParams, SequenceDB)
loadDatabase cp fp = do
  (cp',db) <- case takeExtensions fp of
    ext | ".fasta" `isPrefixOf` ext -> (cp,) `liftM` makeSeqDB cp fp
    ext | ".index" `isPrefixOf` ext -> readIndex cp fp
    _                               -> error ("Unsupported database type: " ++ show fp)

  when (verbose cp) $ do
    hPutStrLn stderr $ "Database: " ++ fp
    hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
    hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
    hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

  return (cp',db)


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> HostModInfo -> DeviceModInfo ->  DeviceSeqDB -> FilePath -> IO ()
search cp db hmi dmi dev fp =
  --readMS2Data fp >>= \r -> case r of
    --Left  s -> hPutStrLn stderr s
    --Right d -> forM_ d $ \ms2 -> do
      --(t,matches) <- bracketTime $ searchForMatches cp db dev hmi dmi ms2
      --when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)

      --printConfig cp fp ms2
      --printResults           $! take (numMatches cp)       matches
      --printResultsDetail     $! take (numMatchesDetail cp) matches
      ----printIonMatchDetail cp $! take (numMatchesIon cp)    matches
  --`catch`
 -- \(e :: SomeException) -> hPutStrLn stderr $ unlines
          --[ "file:   " ++ fp
          --, "reason: " ++ show e ]

--search cp db hmi dmi dev fp =
  readMS2Data fp >>= \r -> case r of
    Left  s -> hPutStrLn stderr s
    Right d -> do 
      (t,allMatches) <- bracketTime $ forM d $ \ms2 -> do
        matches <- searchForMatches cp db dev hmi dmi ms2
        when (showMatchesPerScan cp) $ do
            printConfig cp fp ms2
            printResults           $! take (numMatches cp)       matches
            printResultsDetail     $! take (numMatchesDetail cp) matches
            --printIonMatchDetail cp $! take (numMatchesIon cp)    matches
        return $ map (\m -> (ms2, m)) matches
        `catch`
        \(e :: SomeException) -> do 
                hPutStrLn stderr $ unlines [ "scan:   " ++ (show $ ms2info ms2), "reason: " ++ show e ]
                return []

      when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)

      let n = maximum $ [(numMatches cp), (numMatchesDetail cp), (numMatchesIon cp)]
          toPrint = take n $! sortBy (\(_,a) (_,b) -> compare (scoreXC b) (scoreXC a)) $ concat allMatches
      printAllResults cp toPrint
      `catch`
      \(e :: SomeException) -> hPutStrLn stderr $ unlines
              [ "file:   " ++ fp
              , "reason: " ++ show e ]

