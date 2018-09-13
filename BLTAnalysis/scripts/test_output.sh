#!/bin/bash 

set -e

targName=$1
targSuff=$2
nStart=$3
inputdir="../test/batch/""$4"

nFiles=$(ls -1 ${inputdir}/${targName}_*.root | wc -l)

root.exe -q -b "../macros/sumEvents.cc($nFiles, $nStart, \"$inputdir/$targName\", \"$targSuff\")"
