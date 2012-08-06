{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Execution
-- Copyright : (c) [2012] Kevin Ying
-- License   : BSD
--
-- Functions to organise GPU execution
--
--------------------------------------------------------------------------------

module Execution
  (
    ExecutionPlan(..), withExecutionPlan
  )
  where

import qualified Data.Vector.Unboxed            as U
import qualified Foreign.CUDA.BLAS              as CUBLAS
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Runtime.Stream    as CUDA
import Sequence.Fragment 
import Control.Exception
import Control.Monad


--------------------------------------------------------------------------------
-- The State Token
--------------------------------------------------------------------------------

--
-- A data structure to hold all of the parameters
-- threaded throughout the program
--
data ExecutionPlan = ExecutionPlan
  {
    cublasHandle :: CUBLAS.Handle,
    resIdxSort   :: U.Vector Int,
    cudaStreams  :: [CUDA.Stream]
  }
  --deriving (Show)

--------------------------------------------------------------------------------
-- Actions
--------------------------------------------------------------------------------

withExecutionPlan :: DeviceSeqDB -> Int -> (ExecutionPlan -> IO a) -> IO a
withExecutionPlan ddb numStreams action = do 
  bracket CUBLAS.create CUBLAS.destroy $ \cuHdl -> 
    bracket (replicateM numStreams CUDA.create) (flip forM_ CUDA.destroy) $ \streams -> do
      pep_idx_r_sorted' <- CUDA.peekListArray (numFragments ddb) (devResIdxSort ddb)
      let pep_idx_r_sorted = U.map fromIntegral $ U.fromList pep_idx_r_sorted'
      action $ ExecutionPlan cuHdl pep_idx_r_sorted streams

--------------------------------------------------------------------------------
-- Parse a configuration file
--------------------------------------------------------------------------------

