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

import Foreign.CUDA.BLAS as CUBLAS


--------------------------------------------------------------------------------
-- The State Token
--------------------------------------------------------------------------------

--
-- A data structure to hold all of the parameters
-- threaded throughout the program
--
data ExecutionPlan = ExecutionPlan
  {
    cublasHandle :: CUBLAS.Handle
  }
  --deriving (Show)

--------------------------------------------------------------------------------
-- Actions
--------------------------------------------------------------------------------

withExecutionPlan :: (ExecutionPlan -> IO a) -> IO a
withExecutionPlan action = do 
    handle <- CUBLAS.create
    res <- action $ ExecutionPlan handle
    CUBLAS.destroy handle
    return res

--------------------------------------------------------------------------------
-- Parse a configuration file
--------------------------------------------------------------------------------

