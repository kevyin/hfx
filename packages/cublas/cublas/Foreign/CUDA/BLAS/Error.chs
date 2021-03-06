{-# LANGUAGE CPP                      #-}
{-# LANGUAGE DeriveDataTypeable       #-}
{-# LANGUAGE ForeignFunctionInterface #-}

module Foreign.CUDA.BLAS.Error
  where

-- System
import Data.Typeable
import Control.Exception.Extensible

#include <cublas_v2.h>
{# context lib="cublas" #}


-- Error codes -----------------------------------------------------------------
--
{# enum cublasStatus_t as Status
  { underscoreToCase }
  with prefix="CUBLAS_STATUS" deriving (Eq, Show) #}

-- Describe each error code
--
describe :: Status -> String
describe s = case s of
               NotInitialized   -> "NotInitialized"
               AllocFailed      -> "AllocFailed"
               InvalidValue     -> "InvalidValue"
               ArchMismatch     -> "ArchMismatch"
               MappingError     -> "MappingError"
               ExecutionFailed  -> "ExecutionFailed"
               InternalError    -> "InternalError"
               _                -> "Haskell cublas library error" -- should never occur


-- Exceptions ------------------------------------------------------------------
--
data CUBLASException
  = ExitCode  Status
  | UserError String
  deriving Typeable

instance Exception CUBLASException

instance Show CUBLASException where
  showsPrec _ (ExitCode  s) = showString ("CUBLAS Exception: " ++ describe s)
  showsPrec _ (UserError s) = showString ("CUBLAS Exception: " ++ s)


-- | Raise a CUBLASException in the IO Monad
--
cublasError :: String -> IO a
cublasError s = throwIO (UserError s)


-- | Return the results of a function on successful execution, otherwise throw
-- an exception with an error string associated with the return code
--
resultIfOk :: (Status, a) -> IO a
resultIfOk (status,result) =
    case status of
        Success -> return  result
        _       -> throwIO (ExitCode status)


-- | Throw an exception with an error string associated with an unsuccessful
-- return code, otherwise return unit.
--
nothingIfOk :: Status -> IO ()
nothingIfOk status =
    case status of
        Success -> return  ()
        _       -> throwIO (ExitCode status)

