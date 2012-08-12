module Paths_cuda (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
catchIO = Exception.catch


version :: Version
version = Version {versionBranch = [0,4,1,0], versionTags = []}
bindir, libdir, datadir, libexecdir :: FilePath

bindir     = "/home/kevin/.cabal/bin"
libdir     = "/home/kevin/.cabal/lib/cuda-0.4.1.0/ghc-7.4.1"
datadir    = "/home/kevin/.cabal/share/cuda-0.4.1.0"
libexecdir = "/home/kevin/.cabal/libexec"

getBinDir, getLibDir, getDataDir, getLibexecDir :: IO FilePath
getBinDir = catchIO (getEnv "cuda_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "cuda_libdir") (\_ -> return libdir)
getDataDir = catchIO (getEnv "cuda_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "cuda_libexecdir") (\_ -> return libexecdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
