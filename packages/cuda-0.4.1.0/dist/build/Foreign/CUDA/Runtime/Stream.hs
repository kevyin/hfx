-- GENERATED by C->Haskell Compiler, version 0.16.3 Crystal Seed, 24 Jan 2009 (Haskell)
-- Edit the ORIGNAL .chs file instead!


{-# LINE 1 "./Foreign/CUDA/Runtime/Stream.chs" #-}{-# LANGUAGE ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Foreign.CUDA.Runtime.Stream
-- Copyright : (c) [2009..2011] Trevor L. McDonell
-- License   : BSD
--
-- Stream management routines
--
--------------------------------------------------------------------------------

module Foreign.CUDA.Runtime.Stream (

  -- * Stream Management
  Stream(..),
  create, destroy, finished, block, defaultStream

) where


{-# LINE 21 "./Foreign/CUDA/Runtime/Stream.chs" #-}

-- Friends
import Foreign.CUDA.Runtime.Error
import Foreign.CUDA.Internal.C2HS

-- System
import Foreign
import Foreign.C
import Control.Monad                                    (liftM)


--------------------------------------------------------------------------------
-- Data Types
--------------------------------------------------------------------------------

-- |
-- A processing stream
--
newtype Stream = Stream { useStream :: ((Ptr ()))}
  deriving (Show)


--------------------------------------------------------------------------------
-- Functions
--------------------------------------------------------------------------------

-- |
-- Create a new asynchronous stream
--
create :: IO Stream
create = resultIfOk =<< cudaStreamCreate

cudaStreamCreate :: IO (Status, Stream)
cudaStreamCreate =
  alloca $ \a1' -> 
  cudaStreamCreate'_ a1' >>= \res ->
  peekStream  a1'>>= \a1'' -> 
  let {res' = cToEnum res} in
  return (res', a1'')
{-# LINE 55 "./Foreign/CUDA/Runtime/Stream.chs" #-}


-- |
-- Destroy and clean up an asynchronous stream
--
destroy :: Stream -> IO ()
destroy s = nothingIfOk =<< cudaStreamDestroy s

cudaStreamDestroy :: Stream -> IO (Status)
cudaStreamDestroy a1 =
  let {a1' = useStream a1} in 
  cudaStreamDestroy'_ a1' >>= \res ->
  let {res' = cToEnum res} in
  return (res')
{-# LINE 65 "./Foreign/CUDA/Runtime/Stream.chs" #-}


-- |
-- Determine if all operations in a stream have completed
--
finished   :: Stream -> IO Bool
finished s =
  cudaStreamQuery s >>= \rv -> do
  case rv of
      Success  -> return True
      NotReady -> return False
      _        -> resultIfOk (rv,undefined)

cudaStreamQuery :: Stream -> IO (Status)
cudaStreamQuery a1 =
  let {a1' = useStream a1} in 
  cudaStreamQuery'_ a1' >>= \res ->
  let {res' = cToEnum res} in
  return (res')
{-# LINE 80 "./Foreign/CUDA/Runtime/Stream.chs" #-}


-- |
-- Block until all operations in a Stream have been completed
--
block :: Stream -> IO ()
block s = nothingIfOk =<< cudaStreamSynchronize s

cudaStreamSynchronize :: Stream -> IO (Status)
cudaStreamSynchronize a1 =
  let {a1' = useStream a1} in 
  cudaStreamSynchronize'_ a1' >>= \res ->
  let {res' = cToEnum res} in
  return (res')
{-# LINE 90 "./Foreign/CUDA/Runtime/Stream.chs" #-}


-- |
-- The main execution stream (0)
--
defaultStream :: Stream
defaultStream = Stream nullPtr

--------------------------------------------------------------------------------
-- Internal
--------------------------------------------------------------------------------

peekStream :: Ptr ((Ptr ())) -> IO Stream
peekStream = liftM Stream . peek



foreign import ccall unsafe "Foreign/CUDA/Runtime/Stream.chs.h cudaStreamCreate"
  cudaStreamCreate'_ :: ((Ptr (Ptr ())) -> (IO CInt))

foreign import ccall unsafe "Foreign/CUDA/Runtime/Stream.chs.h cudaStreamDestroy"
  cudaStreamDestroy'_ :: ((Ptr ()) -> (IO CInt))

foreign import ccall unsafe "Foreign/CUDA/Runtime/Stream.chs.h cudaStreamQuery"
  cudaStreamQuery'_ :: ((Ptr ()) -> (IO CInt))

foreign import ccall unsafe "Foreign/CUDA/Runtime/Stream.chs.h cudaStreamSynchronize"
  cudaStreamSynchronize'_ :: ((Ptr ()) -> (IO CInt))