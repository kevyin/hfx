--------------------------------------------------------------------------------
-- 
-- Module    : Util.Statistics
--
-- Functions for calculating various statistics
--
--------------------------------------------------------------------------------

module Util.Statistics (evalues) where

import System.IO
import Config
import Data.List
import qualified Data.Vector.Unboxed as U
import Statistics.LinearRegression
import Control.Monad

type Bin  = Double
type Plot = ([Double], [Double])
--------------------------------------------------------------------------------
--
-- Calculating E-values (Expectation values)
--
-- An e-value represents the number of expected matches to occur in a search
-- by chance. The lower the E-value the more significant the score
-- E(x) = n * s(x)
-- where n is the total number of peptides scored and s is the survival function
-- of p(x) the probability distribution of x.
--
-- let binSize = evalueGroupLevel 
--
--------------------------------------------------------------------------------
-- Takes advantage that sc is sorted in ascending order

evalues :: ConfigParams -> [Double] -> IO [Double]
evalues cp sc = do
    (alpha, beta) <- linearFitLine cp sc
    --let ev = map (\x -> let x' = alpha + beta * (logBase 10 x) in (10**x') * fromIntegral(length sc)) sc 
    let ev = map (\x -> let x' = alpha + beta * (x/10000) in (10**x')) sc 

    when (verbose cp) $ do
        hPutStrLn stderr $ "Evalues"
        hPutStrLn stderr $ "Alpha = " ++ shows alpha " " ++
                           "Beta = "  ++ show  beta
        --outputScores "evalues.csv" ev
    return ev

linearFitLine :: ConfigParams -> [Double] -> IO (Double, Double)
linearFitLine cp sc = do
    when (verbose cp) $ do
        putStrLn $ "Number of peptides scored " ++ (show n)
        --putStrLn $ "cut " ++ (show cut)
        ----outputScores "scores.csv" sc
        outputCSV "hist.csv"    (x,h)
        --outputCSV "logxVh.csv"    (logx,h)
        outputCSV "xVlogh.csv"    (x,logh)
        --outputCSV "logxVlogh.csv"    (logx,logh)
        --outputCSV "probD.csv"    probD
        --outputCSV "probVlnx.csv"   (lnx, (snd probD))
        --outputCSV "survivalFn.csv"  survF
        --outputCSV "lnxVs.csv" (lnx, (snd survF))
        --outputCSV "lnxVlogs.csv" (lnx, logs)
        ----outputCSV "logx_Vlogs_.csv" (logx',logs')
        --outputCSV "logxVlogs.csv" (logx,logs)
        --outputCSV "logxVs.csv" (logx,(snd survF))
        --outputCSV "xVlogs.csv" ((fst survF),logs)
    return $ linearRegression (U.fromList x) (U.fromList logh)
    where 
          n      = length $ sc  -- total number of peptides scored
          cut    = round ((1 - (evalueCutOffPerc cp) / 100) * fromIntegral (n))
          --sc'    = filter (> 0) $ take cut sc
          sc''    = map (\x -> x + 1 + (-1) * (head sc)) $ take cut sc
          sc'    = filter (flip (>) 4000) sc''
          (x', h) = histogram cp sc'
          x      = map (\x -> x/10000) x'
          --logx   = map (logBase 10) x
          logh   = map (logBase 10) h
          
          --probD  = probHist cp n sc'
          --survF  = survivalFn probD
          ----cut2   = maybe 0 id $ findIndex (< -0.5) logs
          --lnx    = map (log) (fst survF)
          --logx   = map (logBase 10) (fst survF)
          --logs   = map (logBase 10) (snd survF)
          --logx'  = drop cut2 logx
          --logs'  = drop cut2 logs

histogram cp xs = (x, h)
    where
        binWidth  = fromIntegral $ evalueGroupLevel cp
        halfWidth = binWidth / 2
        (x, c)    = unzip $ uniqCount (map (valBin.ithBin) xs)
        h         = map (\x -> (x / (binWidth))) c
        -- given the value, which bin does it belong
        ithBin x = fromIntegral $ ceiling $ x / binWidth
        -- given the ith bin, what is it's middle value
        valBin i = i * binWidth - halfWidth

-- probability density estimation using histograms
-- bins are numbered from 1 and start at the origin
probHist :: (Integral a) => ConfigParams -> a -> [Double] -> Plot
probHist cp n xs = (x, p)
    where
        n'        = fromIntegral n
        binWidth  = fromIntegral $ evalueGroupLevel cp
        halfWidth = binWidth / 2
        (x, c)    = unzip $ uniqCount (map (valBin.ithBin) xs)
        p         = map (\x -> (x / (n'*binWidth))) c
        -- given the value, which bin does it belong
        ithBin x = fromIntegral $ ceiling $ x / binWidth
        -- given the ith bin, what is it's middle value
        valBin i = i * binWidth - halfWidth

-- Survival function s(x)
-- Assumes x is sorted in ascending order
-- Since x is likely to be sparse, fill the missing x values and corresponding s values
survivalFn :: Plot -> Plot
--survivalFn (x, p) = unzip $ fill $ zip x s
survivalFn (x, p) = (x, s)
    --unzip $ fill $ zip x s
    where
        -- s is effectively reverse cumulative sum (from right to left) of p 
        -- assuming list is in ascending order
        s = scanr1 (+) p 

-- Probability distribution p(x)
-- Assumes sc is sorted in ascending order
probDistn :: (Integral a) => a -> Plot -> Plot
probDistn n (x, freq) = let n' = fromIntegral n in (x, (map (\f -> f / n') freq))

-- Frequency histogram f(x)
freqHist :: [Double] -> Plot
freqHist sc = unzip $ uniqCount (map (fromIntegral . round) sc)

-- Assumes list is sorted
uniqCount :: (Eq a, Num b) => [a] -> [(a,b)]
uniqCount l = uniq' 1 l  
    where uniq' _ [] = []
          uniq' count (x:xs)  
              | x `elem` xs = uniq' (count + 1) xs  
              | otherwise   = (x, count) : uniq' 1 xs  

--decRound n d = fromIntegral (round (n * factor)) / factor
    --where factor = 10 ** fromIntegral d

--sigFig n s = fromIntegral (round (n * factor)) / factor
    --where shift  = s - (floor (logBase 10 n) + 1)
          --factor = 10 ** fromIntegral shift


outputCSV :: (Show a, Show b) => FilePath -> ([a],[b]) -> IO ()
outputCSV fp (x, y) = outputCSV' x y 
    where
        outputCSV' (a:as) (b:bs) = do
            appendFile fp $ (show a) ++ "," ++ (show b) ++ "\n"
            outputCSV' as bs
        outputCSV' _      _      = putStrLn $ "wrote to " ++ (show fp)

outputScores :: (Show a) => FilePath -> [a] -> IO ()
outputScores fp (x:xs) = do
    appendFile fp $ (show x) ++ "\n"
    outputScores fp xs
outputScores fp _    = putStrLn $ "wrote to " ++ (show fp)

-- produce a list of typles with first element enumerated from x to x'
fillWith x x' y = fg [(x)..(x')] y
    where
        fg (i:is) j = (i, j) : fg is j
        fg [] _ = []

-- returns: [(-3,40, (-2,30), (-1,30), (0,30), (1,30), (2,20), (3,20), (4,20), (5,20), (6,10))]
fill :: (Enum a, Num a, Ord a, Show a, Show b) => [(a, b)] -> [(a, b)]
fill [p] = [p]
fill a@(p@(x, y):ps) = let fs@((x', y'):_) = fill ps in fs `seq` 
    if (x < x' - 1) then
        --fill $ p : ((x' - 1), y') : fs  very recursive form
        p : (fillWith (x + 1) (x' - 1) y') ++ fs
    else
        p : fs
    --where
        --fs@((x', y'):_) = fill ps
fill [] = error ("Statistics.hs: fill encountered empty list")

