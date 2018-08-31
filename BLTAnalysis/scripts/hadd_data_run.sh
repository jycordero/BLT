#! /bin/sh

set -e


nTarget=$1
srcPref=$2
dataName=$3
dataTag=$4
srcDir=$5

# Don't try to do anything if the wrong version of ROOT is loaded
myROOTSYS="/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt"
if [ "$ROOTSYS" != "$myROOTSYS" ]
then
    echo "ROOT 6.10 not loaded"
else

    srcName=${srcPref}_${dataName}${dataTag}
    targName=${dataName}${dataTag}

    files=($(ls -1 ${srcDir}/${srcName}_*.root))
    nSource=${#files[@]}
    nHadd=$((nSource / nTarget))
    remainder=$((nSource - (nTarget * nHadd)))

    echo "Source directory contains ${nSource} \"${srcName}\" files"
    if [ ${remainder} -ne 0 ]
    then
        nHadd1=$nHadd
        nHadd2=$((nHadd + 1))
        echo "${nTarget} target \"${targName}\" files => ${nHadd1} or ${nHadd2} source \"${srcName}\" files each"
    else
        nHadd1=$nHadd
        echo "${nTarget} target \"${targName}\" files => ${nHadd1} source \"${srcName}\" files each"
    fi


    startIndex=0

    for ((i=0; i<nTarget; i+=1))
    do
        echo ""

        # Pick nHadd to ensure most even file splitting
        if [ "$remainder" -ne 0 ]
        then
            nHadd=$nHadd2
            remainder=$((remainder - 1))
        else
            nHadd=$nHadd1
        fi

        rootFile=${targName}_${i}.root

        hadd ${rootFile} ${files[@]:startIndex:nHadd}
        root.exe -q -b "../macros/renameTree.cc(\"${rootFile}\", \"tree_${targName}\", \"tree_${dataName}\")"

        if rootrm ${rootFile}:tree_${targName}
        then
            echo "Deleted tree_${targName}"
        fi
        if rootmv ${rootFile}:TotalEvents_${targName} ${rootFile}:TotalEvents_${dataName}
        then
            echo "Moved TotalEvents_${targName} to TotalEvents_${dataName}"
        fi
        if rootmv ${rootFile}:AcceptedEvents_${targName} ${rootFile}:AcceptedEvents_${dataName}
        then
            echo "Moved AcceptedEvents_${targName} to AcceptedEvents_${dataName}"
        fi

        startIndex=$((startIndex + nHadd))
    done
fi
