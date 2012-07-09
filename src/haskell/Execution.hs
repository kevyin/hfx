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

--import Mass
--import Util.Misc
--import Util.Parsec

--import Control.Monad
--import Data.Char
--import Data.List
--import Data.Maybe
--import Data.Version
--import System.Console.GetOpt
--import System.Directory
--import System.IO
--import System.Exit
--import System.Environment
--import Text.Show.Functions ()
--import Text.ParserCombinators.Parsec

--import Data.Ix
--import Data.Vector.Unboxed (Vector)
--import qualified Data.Vector.Unboxed as U

--import Paths_hfx (version)
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

