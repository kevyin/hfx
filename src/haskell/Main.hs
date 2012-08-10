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
import Execution                        (withExecutionPlan)

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
import qualified Foreign.CUDA.Algorithms as CUDA    (prepareIons)  


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
  --(cp',dbs) <- loadDatabase cp fp
  (t,(cp',db)) <- bracketTime $ loadDatabase cp fp
  when (verbose cp) $ hPutStrLn stderr ("Load Database Elapsed time: " ++ showTime t)

  when (verbose cp) $ hPutStrLn stderr ("Searching ...\n" )
  let hmi = makeModInfo cp'
  --(t2,allMatches') <- bracketTime $ forM dbs $ \db -> do
  (t2, allMatches') <- bracketTime $ withDeviceDB cp' db $ \ddb -> 
      withDevModInfo ddb hmi $ \dmi -> do
        case dmi of
          NoDMod -> CUDA.prepareIons (devIons ddb) (numIons ddb) CUDA.nullDevPtr 0 
          _      -> let (num_ma, d_ma, _) = devModAcids dmi in do
            CUDA.prepareIons (devIons ddb) (numIons ddb) d_ma num_ma
        forM dta $ \f -> do
            matches <- search cp' db hmi dmi ddb f
            when (showMatchesPerScan cp) $ printScanResults cp' f matches
            return matches
  when (verbose cp) $ hPutStrLn stderr ("Search Elapsed time: " ++ showTime t2)
          

  let sortXC = sortBy (\(_,_,a) (_,_,b) -> compare (scoreXC b) (scoreXC a))
      matchesRaw = concat $ concat allMatches'
      --matchesByScan = map sortXC $ groupBy (\(_,ms,_) (_,ms',_) -> ms == ms') matchesRaw
      --matchesByFile = map sortXC $ groupBy (\(f,_,_) (f',_,_) -> f == f') matchesRaw
      n = maximum $ [(numMatches cp'), (numMatchesDetail cp'), (numMatchesIon cp')]
      allMatchesN = take n $! sortXC $ matchesRaw

  printAllResults cp' allMatchesN
    


{-# INLINE loadDatabase #-}
loadDatabase :: ConfigParams -> FilePath -> IO (ConfigParams, SequenceDB)
loadDatabase cp fp = do
  (cp',dbs) <- case takeExtensions fp of
    ext | ".fasta" `isPrefixOf` ext -> (cp,) `liftM` makeSeqDB cp fp
    ext | ".index" `isPrefixOf` ext -> readIndex cp fp >>= \(cp',sdb) -> return (cp',sdb) 
    _                               -> error ("Unsupported database type: " ++ show fp)
  return (cp',dbs)


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> HostModInfo -> DeviceModInfo ->  DeviceSeqDB -> FilePath -> IO [[(FilePath,MS2Data,Match)]]
search cp db hmi dmi dev fp =
  readMS2Data fp >>= \r -> case r of
    Left  s -> hPutStrLn stderr s >>= \_ -> return []
    Right d -> withExecutionPlan dev (maxConcurrentScans cp) $ \ep -> do
      when (verbose cp) $ do
        --hPutStrLn stderr $ "Database: " ++ fp
        hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
        hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
        hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

      let groupByN n xs = if n <= length xs 
                         then let (a, b) = splitAt n xs in [a] ++ groupByN n b
                         else if length xs > 0 then [xs] else []
          -- groupByN n _  = []
      -- allMatches <- forM (groupBy2 d) $ \ms2s -> do
      allMatches <- forM (groupByN (maxConcurrentScans cp) d) $ \ms2s -> do
        forM_ ms2s $ \ms2 -> do
            hPutStrLn stderr $ " searching scan: " ++ (show $ ms2info ms2)
        matches <- searchForMatches cp ep db dev hmi dmi ms2s
        return $ flip map matches (\(ms2,mc) -> map (\m -> (fp, ms2, m)) mc)
        `catch`
        \(e :: SomeException) -> do 
                hPutStrLn stderr $ unlines [ "First scan in group:   " ++ (show $ map ms2info $ ms2s ), "reason: " ++ show e ]
                return []

      return $ concat $ allMatches

      `catch`
      \(e :: SomeException) -> do 
          hPutStrLn stderr $ unlines
              [ "file:   " ++ fp
              , "reason: " ++ show e ]
          return []
