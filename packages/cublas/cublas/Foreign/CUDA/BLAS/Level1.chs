{-# LANGUAGE CPP                      #-}
{-# LANGUAGE ForeignFunctionInterface #-}

module Foreign.CUDA.BLAS.Level1 (

  -- * Single-precision float
  sasum, sscal,

) where

-- Friends
import Foreign.CUDA.BLAS.Error
import Foreign.CUDA.BLAS.Context
import Foreign.CUDA.BLAS.Internal.C2HS

import Foreign.CUDA.Ptr

-- System
import Foreign
import Foreign.C

#include <cublas_v2.h>
{# context lib="cublas" #}


-- Level-1 BLAS operations -----------------------------------------------------
--

-- | This function computes the sum of the absolute values of the elements of
-- vector x.
--
sasum :: Handle -> DevicePtr Float -> Int -> Int -> IO Float
sasum ctx dx n inc = resultIfOk =<< cublasSasum ctx n dx inc

{# fun unsafe cublasSasum_v2 as cublasSasum
  { useHandle   `Handle'
  ,             `Int'
  , useDev      `DevicePtr Float'
  ,             `Int'
  , alloca-     `Float' peek'*  } -> `Status' cToEnum #}
  where
    useDev      = useDevicePtr . castDevPtr
    peek'       = peek . castPtr


-- | This function scales the vector x by the scalar α and overwrites it with
-- the result.
--
sscal :: Handle -> Float -> DevicePtr Float -> Int -> Int -> IO ()
sscal hdl alpha xs nx incx = nothingIfOk =<< cublasSscal hdl nx alpha xs incx

{# fun unsafe cublasSscal_v2 as cublasSscal
  { useHandle   `Handle'
  ,             `Int'
  , with'*      `Float'
  , withDev     `DevicePtr Float'
  ,             `Int'                   } -> `Status' cToEnum #}
  where
    withDev     = useDevicePtr . castDevPtr
    with' x f   = with x (f . castPtr)

