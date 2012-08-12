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
    ExecutionPlan(..), withExecutionPlan,
    ScorePlan(..), groupScorePlan
  )
  where

import Sequence.Fragment 
import Control.Exception
import Control.Monad
import Data.List

import qualified Data.Vector.Unboxed            as U
import qualified Foreign.CUDA.BLAS              as CUBLAS
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Runtime.Stream    as CUDA


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
      action $ ExecutionPlan cuHdl streams


--------------------------------------------------------------------------------
-- Scoring spectrums
--------------------------------------------------------------------------------
-- Aids in planning the task of making spectrum and scoring
-- 
-- | A list of tuples 
--  (ith spectrum, 
--   starting pos of the peptides for this spec, 
--   number of peptides to score)
type ScorePlan = [(Int,Int,Int)]

--
-- |Groups plans to use as much of the memory on a card as possible
-- maxN is the maximum number of peptides able to be scored at one time
-- The input plan is expected to be a list containing the
--      spectrum,
--      starting postition of all its peptides (0)
--      number of peptides this spectrum has
-- This function will group spectrums together to fill maxN,
-- A Plan for one spectrum will only be separated when it has more peptides than
-- the maximum a card can handle
groupScorePlan :: Int -> ScorePlan -> [ScorePlan]
groupScorePlan maxN specs = groupScorePlan' specs 
  where
    groupScorePlan' ss = let (_,_,ns) = unzip3 ss
                             n_scan = scanl1 (+) ns in
                   case findIndex (maxN <) n_scan of
                     Nothing     -> [ss]
                     Just cutoff -> if cutoff /= 0 
                                    then let (grp, rest) = splitAt cutoff ss in
                                         [grp] ++ groupScorePlan' rest
                                    else groupScorePlan' ((splitS (head ss)) ++ (tail ss))
    groupScorePlan' _ = []

    splitS (spec,begin,n) = 
        let ns = ceiling ((fromIntegral maxN) / (fromIntegral n)) in
        flip map [0..ns-1] (\i -> let offs = maxN*i in 
                (spec,begin+offs,(if offs <= n then maxN else n - maxN))) 


