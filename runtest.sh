#!/bin/bash

dbFile=$1

filename="${dbFile##*/}"                      # Strip longest match of */ from start
dir="${dbFile:0:${#dbFile} - ${#filename}}" # Substring from 0 thru pos of filename
base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
ext="${filename:${#base} + 1}"                  # Substring from len of base thru end
if [[ -z "$base" && -n "$ext" ]]; then          # If we have an extension and no base, it's really the base
  base=".$ext"
  ext=""
fi

indexFile="$base.index"
tmpdbFile="tmpdb.out"
tmpinFile="tmpin.out"

echo "Generating $indexFile"
time ./dist/build/generate-peptides/generate-peptides -v -p hfx.params -d $dbFile -o $indexFile

echo "running normal database to $tmpdbFile"
time ./dist/build/hfx/hfx -v -p hfx.params -d $dbFile ./tests/data/test.dta > $tmpdbFile

echo "running normal database to $tmpinFile"
time ./dist/build/hfx/hfx -v -p hfx.params -d $indexFile ./tests/data/test.dta > $tmpinFile

diff $tmpdbFile $tmpinFile
