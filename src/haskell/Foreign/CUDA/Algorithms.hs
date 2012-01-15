{-# LANGUAGE CPP, ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Foreign.CUDA.Algorithms
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- FFI Bindings to a collection of algorithms, implemented to run on the
-- graphics processor using CUDA.
--
--------------------------------------------------------------------------------

module Foreign.CUDA.Algorithms
  (
    findIndicesInRange, findModablePeptides, findBeginEnd,
    calcTotalModCands, prepareGenMod, genModCands,
    addIons, addModIons, 
    rsort,
    mvm
  )
  where

import Util.C2HS
import Data.Word
import Foreign
import Foreign.CUDA                             (DevicePtr, withDevicePtr)


findIndicesInRange :: DevicePtr Float -> DevicePtr Word32 -> Int -> Float -> Float -> IO Int
findIndicesInRange a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  cIntConv `fmap` findIndicesInRange'_ a1' a2' (cIntConv a3) a4 a5

foreign import ccall unsafe "algorithms.h findIndicesInRange_f"
  findIndicesInRange'_ :: Ptr Float -> Ptr Word32 -> Word32 -> Float -> Float -> IO Word32
 
findBeginEnd :: DevicePtr Word32 -> DevicePtr Word32 -> DevicePtr Word32 -> DevicePtr Word32 -> DevicePtr Float -> DevicePtr Word32 -> Int -> DevicePtr Float -> Int -> Float -> Float -> IO Int
findBeginEnd a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  --withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  cIntConv `fmap` findBeginEnd'_ a1' a2' a3' a4' a5' a6' (cIntConv a7) a8' (cIntConv a9) a10 a11

foreign import ccall unsafe "algorithms.h findBeginEnd_f"
  findBeginEnd'_ :: Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Float -> Ptr Word32 -> Word32 -> Ptr Float -> Word32 -> Float -> Float -> IO Int

  --CUDA.findModablePeptides d_pep_idx_valid d_pep_idx d_pep_mod_idx d_pep_ma_count (devIons ddb) (devTerminals ddb) d_pep_idx_r_sorted d_begin d_end d_num_pep_scan d_mod_ma_count num_mod d_ma num_ma >>= \n -> do
findModablePeptides :: DevicePtr Word32                     --- ^ result array, indices to modable peps
                    -> DevicePtr Word32                     --- ^ result array, indices to the corresponding peptide
                    -> DevicePtr Word32                     --- ^ result array, indices to corres mod comb
                    -> DevicePtr Word32                     --- ^ result array, pep ma counts
                    -> Int                                  --- ^ number of total peptide candidates by mass 
                    -> DevicePtr Word8                      --- ^ amino acid ion db
                    -> (DevicePtr Word32, DevicePtr Word32) --- ^ c and n terminals

                    -> DevicePtr Word32                     --- ^ d_pep_idx_r_sorted peptides sorted by residual

                    -> DevicePtr Word32                     --- ^ d_begin : start of cands
                    -> DevicePtr Word32                     --- ^ d_end : end of cands
                    -> DevicePtr Word32                     --- ^ d_num_pep_scan : num of cand for this mod comb
                    -> DevicePtr Word32                     --- ^ d_mod_ma_count : ma counts for this modcomb
                    -> Int                                  --- ^ total mod combs


                    -> DevicePtr Word8                      --- ^ ma - modable acids 
                    -> Int                                  --- ^ number of ma's
                    -> IO Int
findModablePeptides a1 a2 a3 a4 a5 a6 (a7, a8) a9 a10 a11 a12 a13 a14 a15 a16 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  --withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  withDevicePtr a9 $ \a9' ->
  withDevicePtr a10 $ \a10' ->
  withDevicePtr a11 $ \a11' ->
  withDevicePtr a12 $ \a12' ->
  withDevicePtr a13 $ \a13' ->
  --withDevicePtr a14 $ \a14' ->
  withDevicePtr a15 $ \a15' ->
  --withDevicePtr a16 $ \a16' ->
  cIntConv `fmap` findModablePeptides'_ a1' a2' a3' a4' (cIntConv a5) a6' a7' a8' a9' a10' a11' a12' a13' (cIntConv a14) a15' (cIntConv a16)

foreign import ccall unsafe "algorithms.h findModablePeptides"
  findModablePeptides'_ :: Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Ptr Word8 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Ptr Word8 -> Word32 -> IO Word32

calcTotalModCands :: DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> DevicePtr Word32
                  -> Int
                  -> Int
                  -> Int
                  -> IO Int
calcTotalModCands a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  cIntConv `fmap` calcTotalModCands'_ a1' a2' a3' a4' a5' a6' a7' a8' (cIntConv a9) (cIntConv a10) (cIntConv a11)

foreign import ccall unsafe "algorithms.h calcTotalModCands"
  calcTotalModCands'_ :: Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO Word32

prepareGenMod :: DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> DevicePtr Word32
              -> Int
              -> Int
              -> IO Int
prepareGenMod a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  withDevicePtr a9 $ \a9' ->
  withDevicePtr a10 $ \a10' ->
  cIntConv `fmap` prepareGenMod'_ a1' a2' a3' a4' a5' a6' a7' a8' a9' a10' (cIntConv a11) (cIntConv a12)

foreign import ccall unsafe "algorithms.h prepareGenMod"

  prepareGenMod'_ :: Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> IO Word32

genModCands :: DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> DevicePtr Word32
            -> Int
            -> Int
            -> IO ()
genModCands a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  withDevicePtr a9 $ \a9' ->
  withDevicePtr a10 $ \a10' ->
  --withDevicePtr a11 $ \a11' ->
  --withDevicePtr a12 $ \a12' ->
  genModCands'_ a1' a2' a3' a4' a5' a6' a7' a8' a9' a10' (cIntConv a11) (cIntConv a12)

foreign import ccall unsafe "algorithms.h genModCands"
  genModCands'_ :: Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> IO ()

addModIons :: DevicePtr Word32 -> DevicePtr Float -> DevicePtr Float -> DevicePtr Word8 -> (DevicePtr Word32, DevicePtr Word32) -> DevicePtr Word32 -> DevicePtr Word32 -> Int -> DevicePtr Word8 -> DevicePtr Word8 -> DevicePtr Float -> Int -> Int -> Int -> IO ()
addModIons a1 a2 a3 a4 (a5,a6) a7 a8 a9 a10 a11 a12 a13 a14 a15 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  withDevicePtr a8 $ \a8' ->
  --withDevicePtr a9 $ \a9' ->
  withDevicePtr a10 $ \a10' ->
  withDevicePtr a11 $ \a11' ->
  withDevicePtr a12 $ \a12' ->
  --withDevicePtr a13 $ \a13' ->
  --withDevicePtr a14 $ \a14' ->
  --withDevicePtr a15 $ \a15' ->
  addModIons'_ a1' a2' a3' a4' a5' a6' a7' a8' (cIntConv a9) a10' a11' a12' (cIntConv a13) (cIntConv a14) (cIntConv a15)

foreign import ccall unsafe "algorithms.h addModIons"
  addModIons'_ :: Ptr Word32 -> Ptr Float -> Ptr Float -> Ptr Word8 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Ptr Word8 -> Ptr Word8 -> Ptr Float-> Word32 -> Word32 -> Word32 -> IO ()

addIons :: DevicePtr Word32 -> DevicePtr Float -> DevicePtr Float -> DevicePtr Word8 -> (DevicePtr Word32, DevicePtr Word32) -> DevicePtr Word32 -> Int -> Int -> Int -> IO ()
addIons a1 a2 a3 a4 (a5,a6) a7 a8 a9 a10 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  withDevicePtr a7 $ \a7' ->
  addIons'_ a1' a2' a3' a4' a5' a6' a7' (cIntConv a8) (cIntConv a9) (cIntConv a10)

foreign import ccall unsafe "algorithms.h addIons"
  addIons'_ :: Ptr Word32 -> Ptr Float -> Ptr Float -> Ptr Word8 -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO ()


#if 0
addIonsIP
    :: DevicePtr Float                          -- [out] sequence scores
    -> DevicePtr Float                          -- experimental spectrum
    -> DevicePtr Float                          -- residual masses
    -> DevicePtr Float                          -- individual ion masses
    -> (DevicePtr Word32, DevicePtr Word32)     -- c- and n- terminal indices
    -> DevicePtr Word32                         -- indices of the sequences under consideration
    -> Int                                      -- number of sequences to consider
    -> Int                                      -- peptide charge state
    -> Int                                      -- length of input spectrum
    -> IO ()
addIonsIP d_score d_spec d_residual d_ions (d_tc, d_tn) d_idx num_idx max_charge len_spec =
  withDevicePtr d_score    $ \a1' ->
  withDevicePtr d_spec     $ \a2' ->
  withDevicePtr d_residual $ \a3' ->
  withDevicePtr d_ions     $ \a4' ->
  withDevicePtr d_tc       $ \a5' ->
  withDevicePtr d_tn       $ \a6' ->
  withDevicePtr d_idx      $ \a7' ->
  addIons_ip'_ a1' a2' a3' a4' a5' a6' a7' (cIntConv num_idx) (cIntConv max_charge) (cIntConv len_spec)

foreign import ccall unsafe "algorithms.h addIons_inplace"
  addIons_ip'_ :: Ptr Float -> Ptr Float -> Ptr Float -> Ptr Float -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO ()
#endif


rsort :: DevicePtr Float -> DevicePtr Word32 -> Int -> IO ()
rsort a1 a2 a3 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  rsort'_ a1' a2' (cIntConv a3)

foreign import ccall unsafe "algorithms.h sort_rf"
  rsort'_ :: Ptr Float -> Ptr Word32 -> Word32 -> IO ()


mvm :: DevicePtr Float -> DevicePtr Word32 -> DevicePtr Float -> Int -> Int -> IO ()
mvm a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  mvm'_ a1' a2' a3' (fromIntegral a4) (fromIntegral a5)

foreign import ccall unsafe "algorithms.h mvm_if"
  mvm'_ :: Ptr Float -> Ptr Word32 -> Ptr Float -> Word32 -> Word32 -> IO ()

--sum_Word32 :: DevicePtr Word32 -> Int -> IO Int
--sum_Word32 a1 a2 =
  --withDevicePtr a1 $ \a1' ->
  ----withDevicePtr a2 $ \a2' ->
  --cIntConv `fmap` sum_Word32'_ a1' (cIntConv a2)

--foreign import ccall unsafe "algorithms.h sum_Word32"
  --sum_Word32'_ :: Ptr Word32 -> Word32 -> IO Int

