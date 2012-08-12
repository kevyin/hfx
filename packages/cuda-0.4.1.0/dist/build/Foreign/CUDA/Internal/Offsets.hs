{-# LINE 1 "Foreign/CUDA/Internal/Offsets.hsc" #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LINE 2 "Foreign/CUDA/Internal/Offsets.hsc" #-}
{-
 - Structure field offset constants.
 - Too difficult to extract using C->Haskell )=
 -}

module Foreign.CUDA.Internal.Offsets where


--------------------------------------------------------------------------------
-- Runtime API
--------------------------------------------------------------------------------


{-# LINE 15 "Foreign/CUDA/Internal/Offsets.hsc" #-}

devNameOffset, devMaxThreadDimOffset, devMaxGridSizeOffset :: Int

devNameOffset         = (0)
{-# LINE 19 "Foreign/CUDA/Internal/Offsets.hsc" #-}
devMaxThreadDimOffset = (292)
{-# LINE 20 "Foreign/CUDA/Internal/Offsets.hsc" #-}
devMaxGridSizeOffset  = (304)
{-# LINE 21 "Foreign/CUDA/Internal/Offsets.hsc" #-}


{-# LINE 23 "Foreign/CUDA/Internal/Offsets.hsc" #-}
devMaxTexture2DOffset, devMaxTexture3DOffset :: Int
devMaxTexture2DOffset = (384)
{-# LINE 25 "Foreign/CUDA/Internal/Offsets.hsc" #-}
devMaxTexture3DOffset = (412)
{-# LINE 26 "Foreign/CUDA/Internal/Offsets.hsc" #-}

{-# LINE 27 "Foreign/CUDA/Internal/Offsets.hsc" #-}

texRefAddressModeOffset, texRefChannelDescOffset :: Int
texRefAddressModeOffset = (8)
{-# LINE 30 "Foreign/CUDA/Internal/Offsets.hsc" #-}
texRefChannelDescOffset = (20)
{-# LINE 31 "Foreign/CUDA/Internal/Offsets.hsc" #-}

--------------------------------------------------------------------------------
-- Driver API
--------------------------------------------------------------------------------


{-# LINE 37 "Foreign/CUDA/Internal/Offsets.hsc" #-}

devMaxThreadDimOffset', devMaxGridSizeOffset' :: Int

devMaxThreadDimOffset' = (4)
{-# LINE 41 "Foreign/CUDA/Internal/Offsets.hsc" #-}
devMaxGridSizeOffset'  = (16)
{-# LINE 42 "Foreign/CUDA/Internal/Offsets.hsc" #-}

