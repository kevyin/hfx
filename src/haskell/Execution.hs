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

import qualified Data.Vector.Unboxed  as U
import Foreign.CUDA.BLAS    as CUBLAS
import Foreign.CUDA         as CUDA
import Sequence.Fragment 
import Control.Exception
import Data.Word


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
    resIdxSort   :: U.Vector Int
  }
  --deriving (Show)

--------------------------------------------------------------------------------
-- Actions
--------------------------------------------------------------------------------

withExecutionPlan :: DeviceSeqDB -> (ExecutionPlan -> IO a) -> IO a
withExecutionPlan ddb action = do 
  bracket CUBLAS.create CUBLAS.destroy $ \handle -> do
    pep_idx_r_sorted' <- CUDA.peekListArray (numFragments ddb) (devResIdxSort ddb)
    let pep_idx_r_sorted = U.map fromIntegral $ U.fromList pep_idx_r_sorted'
    action $ ExecutionPlan handle pep_idx_r_sorted

--------------------------------------------------------------------------------
-- Parse a configuration file
--------------------------------------------------------------------------------

