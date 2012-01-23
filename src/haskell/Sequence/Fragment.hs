{-# LANGUAGE BangPatterns, CPP #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Fragment
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Reading and processing sequence related stuff
--
--------------------------------------------------------------------------------

module Sequence.Fragment
  (
    SequenceDB(..), DeviceSeqDB(..), HostModInfo(..), DeviceModInfo(..), Fragment(..),
    makeSeqDB, withDeviceDB, makeModInfo, withDevModInfo, fraglabel, PepMod, modifyFragment
  )
  where

import Mass
import Config
import Util.Misc
import Sequence.Fasta

import Prelude                                  hiding ( lookup )
import Data.List                                ( unfoldr, mapAccumL, findIndices )
import Data.Word
import Data.Binary
import Data.Vector.Binary                       ()
import Control.Monad
import Control.Applicative
import Data.Vector.Algorithms.Intro             as VA

import qualified Bio.Sequence                   as F
import qualified Data.ByteString.Lazy           as L
import qualified Data.ByteString.Lazy.Char8     as LC

import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Generic.Mutable    as GM
import qualified Data.Vector.Fusion.Stream      as S
import qualified Data.Vector.Fusion.Stream.Size as S

import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Util              as CUDA

import Debug.Trace

#define PHASE_STREAM [1]
#define PHASE_INNER  [0]

#define INLINE_STREAM INLINE PHASE_STREAM
#define INLINE_INNER  INLINE PHASE_INNER

type PepMod = [(Char, Int, Float)]

--------------------------------------------------------------------------------
-- Sequence Fragments
--------------------------------------------------------------------------------

data Fragment = Fragment
  {
    fragmass   :: Float,                -- Total mass of this fragment
    fragheader :: L.ByteString,         -- Full header describing this fragment
    fragdata   :: L.ByteString          -- Fragment sequence data, including flanking residuals
  }
  deriving (Eq, Show)

--
-- Fragment label (first word of full header)
--
fraglabel :: Fragment -> L.ByteString
fraglabel = head . LC.words . fragheader

--
-- fragdata modified 
--
--modFragData :: PepMod -> [Int] -> L.ByteString -> L.ByteString
--modFragData p u fd =

modifyFragment :: PepMod -> [Int] -> Fragment -> Fragment
modifyFragment pm unrank (Fragment fmass fheader fdata') = Fragment newMass fheader newData 
  where
    fdata = LC.unpack fdata'
    newMass = fmass  + (sum $ map (\(_,c,m) -> (fromIntegral c) * m) pm) 
    newData = LC.pack . modData 0  $ fdata

    modData :: Int -> [Char] -> [Char]
    modData i (x:xs) = if elem i modified
                       then x : '*' : modData (i+1) xs
                       else x : modData (i+1) xs
    modData _ _      = []

    modInfo = snd $ mapAccumL (\u_idx (ma,cnt,_) -> (u_idx + cnt, (ma, cnt, sublist u_idx (u_idx+cnt) unrank))) 0 pm

    modified = findModified modInfo :: [Int]

    findModified ((ma,cnt,u):ms) = if cnt > 0
                                   then (map (\ith -> occurences !! ith) u) ++ findModified ms
                                   else findModified ms
      where
        occurences = filter (\i -> 1 < i && i < ((length fdata) - 2)) $ findIndices ((==) ma) fdata -- filter out flanking acids
    findModified _ = []

--------------------------------------------------------------------------------
-- Sequence Database
--------------------------------------------------------------------------------

--
-- A collection of protein sequences
--
data SequenceDB = SeqDB
  {
    dbHeader  :: V.Vector L.ByteString,            -- sequence ion description headers
    dbIon     :: U.Vector Word8,                   -- flattened array of amino character codes
    dbIonSeg  :: U.Vector Word32,                  -- segmenting information for ions
    dbFrag    :: U.Vector (Float, Word32, Word32), -- (residual mass, c-idx, n-idx)
    dbFragSeg :: U.Vector Word32                   -- fragment segmenting information for deriving parent
  }
  deriving Show

--
-- Binary instance of the digested sequence database, for serialising to disk
-- and (hopefully) fast retrieval for later reuse.
--
instance Binary SequenceDB where
  {-# INLINE put #-}
  put (SeqDB h i is f fs) =
    let (r,c,n) = G.unzip3 f
    in  put h >> put i >> put is >> put r >> put c >> put n >> put fs

  {-# INLINE get #-}
  get = liftM5 SeqDB get get get (liftM3 G.zip3 get get get) get

data HostModInfo = 
    NoHMod
  | HostModInfo
  {
    -- num_ma, ma, ma_mass (mass change)
    modAcids         :: (Int, U.Vector Char, U.Vector Float), 
    -- num_mod, mod_ma_count (size: num_ma X num_mod), mod_ma_count_sum, mod_delta)
    modCombs         :: (Int, U.Vector Int, U.Vector Int, U.Vector Float), 
    resIdxSort       :: U.Vector Int
  }
  deriving Show

data DeviceModInfo = 
    NoDMod
  | DevModInfo
  {
    -- num_ma, d_ma, d_ma_mass (mass change)
    devModAcids         :: (Int, CUDA.DevicePtr Word8, CUDA.DevicePtr Float), 
    -- num_mod, d_mod_ma_count (size: num_ma X num_mod), d_mod_ma_count_sum, d_mod_delta)
    devModCombs         :: (Int, CUDA.DevicePtr Word32, CUDA.DevicePtr Word32, CUDA.DevicePtr Float), 
    devResIdxSort       :: CUDA.DevicePtr Word32
  }
  deriving Show

--
-- An extension of the above referencing memory actually stored on the graphics
-- device. Meant to be used in tandem.
--
data DeviceSeqDB = DevDB
  {
    numIons             :: Int,
    numFragments        :: Int,
    devIons             :: CUDA.DevicePtr Word8,
    devMassTable        :: CUDA.DevicePtr Float,
    devResiduals        :: CUDA.DevicePtr Float,
    devTerminals        :: (CUDA.DevicePtr Word32, CUDA.DevicePtr Word32)
  }
  deriving Show


--
-- Generate a sequence database from the given file. This is a flattened
-- representation of ion character codes and fragment information, together with
-- segmenting information describing the boundaries of the original Proteins
-- they derive from.
--
makeSeqDB :: ConfigParams -> FilePath -> Int -> IO [SequenceDB]
{-# INLINE makeSeqDB #-}
makeSeqDB cp fp split = do
  --ns <- countSeqs fp
  db' <- readFasta fp
  let ns'  = length db'
      sections = sectionsOfSplit ns' split

  forM sections $ \(beg, ns) -> do
      let db   = sublist beg ns db'
          hdr  = headers ns db
          iseg = ionSeg  ns db
          ions = ionSeq (fromIntegral (U.last iseg)) db
          nf   = countFrags cp ions

      f  <- GM.new nf               -- maximum number as estimated above, trim later
      fs <- GM.new (ns+1)           -- include the zero index

      --
      -- OPT: make versions that operate directly from (unboxed) streams?
      --
      let write !i !v = do GM.unsafeWrite f i v
                           return (i+1)

      let fill (!i,!n) (!x,!y) = do GM.unsafeWrite fs i (fromIntegral n)
                                    n' <- foldM write n (digest cp ions (x,y))
                                    return (i+1, n')

      --
      -- Now iterate over all of the protein sequences, keeping track of the segment
      -- number and number of fragments generated so far, so we can also fill in the
      -- fragment segmenting information.
      --
      nf' <- snd <$> G.foldM' fill (0,0) (G.zip iseg (G.tail iseg))
      GM.unsafeWrite fs ns (fromIntegral nf')

      f'  <- G.unsafeFreeze (GM.take nf' f)
      fs' <- G.unsafeFreeze fs

      return $ SeqDB hdr ions iseg f' fs'


--
-- Transfer sequence data to the graphics device and execute some computation
--
withDeviceDB :: ConfigParams -> SequenceDB -> (DeviceSeqDB -> IO a) -> IO a
{-# INLINE withDeviceDB #-}
withDeviceDB cp sdb action =
  let (r,c,n) = G.unzip3 (dbFrag sdb)
      mt      = aaMassTable cp
      ions    = dbIon sdb
      numIon  = G.length (dbIon  sdb)
      numFrag = G.length (dbFrag sdb)
  in
    CUDA.withVector r    $ \d_r    ->
    CUDA.withVector c    $ \d_c    ->
    CUDA.withVector n    $ \d_n    ->
    CUDA.withVector mt   $ \d_mt   ->
    CUDA.withVector ions $ \d_ions ->
      action (DevDB numIon numFrag d_ions d_mt d_r (d_c, d_n))

makeModInfo :: ConfigParams -> SequenceDB -> IO HostModInfo
{-# INLINE makeModInfo #-}
makeModInfo cp sdb = if (num_ma < 1) then return NoHMod
  else do
    let (r,_,_) = G.unzip3 (dbFrag sdb)
    pep_idx <- U.unsafeThaw $ U.enumFromN 0 (G.length r) 
    VA.sortBy (\i j -> compare (r G.! i) (r G.! j)) pep_idx
    pep_idx_r_sorted <- U.unsafeFreeze pep_idx

    return $ HostModInfo (num_ma, ma, ma_mass) (num_mod, mod_ma_count, mod_ma_count_sum, mod_delta) pep_idx_r_sorted
      where
        max_ma  = maxModableAcids cp
        ma      = U.fromList $ getMA cp
        ma_mass = U.fromList $ getMA_Mass cp
        num_ma  = U.length ma
        num_mod          = length ma_count
        mod_ma_count     = U.concat ma_count
        mod_ma_count_sum = U.fromList $ map U.sum ma_count
        mod_delta        = U.fromList $ map calcDelta ma_count

        calcDelta ma_cnt = U.sum (U.zipWith (\c m -> m * fromIntegral c) ma_cnt ma_mass)
        ma_count  = map U.fromList $ filter (\c -> if sum c < max_ma then True else False) ma_count'
        ma_count' = tail (addModDimensions ((U.length ma) - 1) [ replicate (U.length ma) 0 ])

        addModDimensions :: Int -> [[Int]] -> [[Int]]
        addModDimensions dim xs = 
            if dim >= 0 then expand $ addModDimensions (dim-1) xs else xs
            where
            expand (c:cs) = [c] ++ 
                            map (\cnt -> replace dim cnt c) [1..max_ma] ++ 
                            expand cs
            expand _      = []

        replace :: (Show a) => Int -> a -> [a] -> [a]
        replace idx val list = x ++ val : ys
            where
            (x,_:ys) = splitAt idx list

withDevModInfo :: HostModInfo -> (DeviceModInfo -> IO a) -> IO a
{-# INLINE withDevModInfo #-}
withDevModInfo hmi action = 
    let (num_ma, ma, ma_mass) = modAcids hmi
        (num_mod, mod_ma_count, mod_ma_count_sum, mod_delta) = modCombs hmi
        pep_idx_r_sorted = resIdxSort hmi
    in case hmi of
      NoHMod    -> action NoDMod
      _         -> 
        CUDA.withVector (U.map c2w ma)   $ \d_ma        ->
        CUDA.withVector ma_mass          $ \d_ma_mass   ->
        CUDA.withVector mod_delta        $ \d_mod_delta ->
        CUDA.withVector (U.map fromIntegral mod_ma_count)     $ \d_mod_ma_count ->
        CUDA.withVector (U.map fromIntegral mod_ma_count_sum) $ \d_mod_ma_count_sum -> 
        CUDA.withVector (U.map fromIntegral pep_idx_r_sorted) $ \d_r_sorted -> 

          action (DevModInfo (num_ma, d_ma, d_ma_mass) (num_mod, d_mod_ma_count, d_mod_ma_count_sum, d_mod_delta) d_r_sorted)
    
--------------------------------------------------------------------------------
-- Digestion
--------------------------------------------------------------------------------

--
-- Estimate the number of fragments that will be produced from a collection of
-- peptide sequences. This is an upper bound, as it does not take into account:
--   * fragments outside of the designated mass bounds
--   * split points at the end of a sequence, producing one fragment and not two
--
-- The latter accounts for a rather small discrepancy (<1%), while the former
-- may result in significant overestimation (>10%).
--
countFrags :: ConfigParams -> U.Vector Word8 -> Int
{-# INLINE countFrags #-}
countFrags cp aa = (missedCleavages cp + 1) * count aa
  where
    rule  = fst (digestionRule cp) . w2c
    count = U.length . U.filter rule


--
-- Extract sequence headers
--
headers :: Int -> [Protein] -> V.Vector L.ByteString
{-# INLINE headers #-}
headers n db
  = G.unstream
  $ S.fromList (map F.seqheader db) `S.sized` S.Exact n


--
-- Extract the amino acid character codes and segmenting information from the
-- given list of proteins
--
ionSeq :: Int -> [Protein] -> U.Vector Word8
{-# INLINE ionSeq #-}
ionSeq n db
  = G.unstream
  $ S.fromList (concatMap (L.unpack . F.seqdata) db) `S.sized` S.Exact n


ionSeg :: Int -> [Protein] -> U.Vector Word32
{-# INLINE ionSeg #-}
ionSeg n db
  = G.unstream
  $ S.fromList (scanl (+) 0 $ map (fromIntegral . L.length . F.seqdata) db) `S.sized` S.Exact (n+1)


--
-- Digest a single protein sequence with the given enzyme and cleavage rules,
-- removing any fragments outside the broad mass range of interest.
--
-- The mass of the peptide is the sum of the amino acid residue masses plus the
-- mass of the water molecule released in forming the peptide bond (plus one;
-- from Eq. 1 of Eng.[1])
--
digest :: ConfigParams -> U.Vector Word8 -> (Word32,Word32) -> [(Float,Word32,Word32)]
{-# INLINE digest #-}
digest cp ions (!x,!y) = filter (inrange . fst3) . splice cp . fragment cp x $ segment
  where
    segment   = G.unsafeSlice (fromIntegral x) (fromIntegral (y-x)) ions
    inrange m = let m' = m + massH2O + massH
                in  minPeptideMass cp <= m' && m' <= maxPeptideMass cp

--
-- Split a protein sequence with the given digestion rule. This outputs the
-- total residual mass of the fragment, together with the indices of the C- and
-- N- terminals in global index coordinates.
--
fragment :: ConfigParams -> Word32 -> U.Vector Word8 -> [(Float,Word32,Word32)]
{-# INLINE fragment #-}
fragment cp gid ions = unfoldr step (gid,ions)
  where
    rule = fst (digestionRule cp) . w2c
    mass = G.foldl' (+) 0 . G.map (getAAMass cp . w2c)

    step (!c,!v)
      | G.null v  = Nothing
      | otherwise = case G.findIndex rule v of
                      Nothing -> Just ((mass v, c, c + fromIntegral (G.length v)), (0, G.empty))
                      Just !i -> let n = fromIntegral i
                                     h = U.take (i+1) v
                                     t = U.drop (i+1) v
                                 in Just ((mass h, c, c+n), (c+n+1,t))


--
-- Generate additional sequences from missed cleavages of sequential, adjacent
-- peptide fragments.
--
splice :: ConfigParams -> [(Float,Word32,Word32)] -> [(Float,Word32,Word32)]
{-# INLINE splice #-}
splice cp = loop
  where
    loop []     = []
    loop (p:ps) = scanl fuse p (take n ps) ++ loop ps

    n        = missedCleavages cp
    fuse a b = (fst3 a + fst3 b, snd3 a, thd3 b)

