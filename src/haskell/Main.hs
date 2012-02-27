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
import Util.Concurrent

--
-- System libraries
--
import Data.List
import Data.Maybe
import Control.Monad
import Control.Exception
import Control.Concurrent               (threadDelay)
import Control.Concurrent.MVar
import System.Environment
import System.FilePath
import System.IO
import Prelude                          hiding ( lookup, catch )
import Criterion.Measurement            as CM

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
  --(t,(cp',dbs)) <- bracketTime $ loadDatabase cp fp (splitDB cp)
  let gpus = 4
  (t,(cp',dbs)) <- bracketTime $ loadDatabase cp fp gpus 
        --when (verbose cp) $ do
          --hPutStrLn stderr $ "Database: " ++ fp
          --hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
          --hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
          --hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)
          --
  when (verbose cp) $ hPutStrLn stderr ("Load Database Elapsed time: " ++ showTime t)

  when (verbose cp) $ hPutStrLn stderr ("Loading Spectrums ...\n" )
  samples <- forM dta $ \f -> do
    r <- readMS2Data f 
    case r of
      Left  s -> hPutStrLn stderr s >>= \_ -> return (f,[])
      Right d -> do 
        return (f,d)

  when (verbose cp) $ hPutStrLn stderr ("Searching ...\n" )
  let hmi = makeModInfo cp'
  plans <- makePlans cp' gpus dbs samples

  t1 <- CM.getTime
  res <- flip forkJoin plans $ \plan -> do
  --res <- forM plans $ \plan -> do
    let db = database plan
    CUDA.set (device plan)
    replicateM 70 $ do
      withDeviceDB cp' db $ \ddb -> 
        threadDelay 1
      --withDevModInfo ddb hmi $ \dmi -> 
        --search plan cp' db hmi dmi ddb
    --threadDelay 5000000
    return (device plan)
  CUDA.sync
  putStrLn $ show res
  --numResPerPlan <- mapM (liftM sum . mapM readMVar . numRes ) plans

  --when (verbose cp) $ hPutStrLn stderr $ show numResPerPlan
  t2 <- CM.getTime
  when (verbose cp) $ hPutStrLn stderr (" GPU processing time: " ++ (show $ CM.secs $ t2 - t1))

  --when (verbose cp) $ hPutStrLn stderr ("Search Elapsed time: " ++ showTime t2)

  --matchesRaw <- retrieveMatches cp' plans
          

  --let sortXC = sortBy (\(_,_,a) (_,_,b) -> compare (scoreXC b) (scoreXC a))
      --matchesByScan = map sortXC $ groupBy (\(_,ms,_) (_,ms',_) -> ms == ms') $ sortBy (\(_,ms,_) (_,ms',_) -> compare (ms2info ms) (ms2info ms')) matchesRaw
      ----matchesByFile = map sortXC $ groupBy (\(f,_,_) (f',_,_) -> f == f') matchesRaw
      --n = maximum $ [(numMatches cp'), (numMatchesDetail cp'), (numMatchesIon cp')]
      --allMatchesN = take n $! sortXC $ matchesRaw

  --when (showMatchesPerScan cp) $ do
      --printScanResults cp' matchesByScan
  --printAllResults cp' allMatchesN
    


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
search :: SearchPlan -> ConfigParams -> SequenceDB -> HostModInfo -> DeviceModInfo ->  DeviceSeqDB -> IO ()
search plan cp db hmi dmi dev =
  forM_ [0..(nsamples plan - 1)] $ \i -> do
    searchForMatches (plan,i) cp db dev hmi dmi
  `catch`
  \(e :: SomeException) -> do 
      hPutStrLn stderr $ unlines
          [ "device: " ++ (show $ device plan)
          , "file:   " ++ fileName plan
          , "reason: " ++ show e ]
