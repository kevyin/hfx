--------------------------------------------------------------------------------
-- |
-- Module    : Util.Concurrent
-- Copyright : (c) [2009..2012] Kevin Ying
-- License   : BSD
--
-- Easy forkOS threads
-- The function of the MVar is not only for returning values.
-- Because child forks terminate when their parents terminate, the MVar
-- ensures the forks are finished before the parent thread resumes
--
-- source:
--
--------------------------------------------------------------------------------

module Util.Concurrent
  (
    forkJoin
  )
  where

import Control.Concurrent             (forkOS)
import Control.Concurrent.MVar 

-- |forkJoin
-- perform each task in an OS fork and returing the results
-- Stores results in MVars
-- source: stackoverflow.com/questions/2233452/
forkJoin :: (a -> IO b) -> [a] -> IO [b]
forkJoin f xs = (fork f xs) >>= joinf

fork1 :: (a -> IO b) -> a -> IO (MVar b)
fork1 f x =
  do 
    cell <- newEmptyMVar
    forkOS (do {result <- f x; putMVar cell result})
    return cell

fork :: (a -> IO b) -> [a] -> IO [MVar b]
fork f = mapM (fork1 f)

joinf :: [MVar b] -> IO [b]
joinf = mapM takeMVar

