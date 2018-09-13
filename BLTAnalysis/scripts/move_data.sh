#!/bin/bash 

set -e

srcName=$1
targName=$2
startIdx=$3
targDir=$4

nFiles=$(ls -1 ${srcName}*.root | wc -l)

echo "Moving ${nFiles} \"${srcName}\" files to ${targDir} as \"${targName}\""

i=$startIdx
for filename in `ls -1 ${srcName}*.root`
do
    echo "Copying ${filename} to ${targDir}/${targName}_${i}.root"
    if xrdcp ${filename} root://cmseos.fnal.gov/${targDir}/${targName}_${i}.root
    then
        echo "Deleting local copy of ${filename}"
        rm ${filename}
    fi
    echo ""

    i=$((i + 1))
done

nFiles=$((i - startIdx))
echo "Moved ${nFiles} files to ${targDir}"
