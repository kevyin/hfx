{-# LANGUAGE CPP, ForeignFunctionInterface, EmptyDataDecls #-}

module Foreign.CUDA.Cublas (
  Handle, create, destroy

)
where

-- #include "cublas_v2.h"

--{# context lib="cublas" #}

import Util.C2HS
-- System
import Foreign
import Foreign.C
import Control.Monad          (liftM)


--------------------------------------------------------------------------------
-- Data Types
--------------------------------------------------------------------------------
-- |
-- Cublas handle
newtype Handle = Handle { useHandle :: {# type cublasHandle_t #}}

--------------------------------------------------------------------------------
-- Event management
--------------------------------------------------------------------------------

-- |
-- Create a handle
create :: IO CublasHandle
create = resultIfOk =<< cublasCreate

{# fun unsafe cublasCreate
  { alloca-     `Handle'    peekHandle* } -> `Status' cToEnum #}
  where peekHandle = liftM Handle . peek

-- |
-- Destroy a handle
destroy :: Handle -> IO ()
destroy ha = nothingIfOk =<< cublasDestroy ha

{# fun unsafe cublasDestroy
  { useHandle `Handle' } -> `Status' cToEnum #}
