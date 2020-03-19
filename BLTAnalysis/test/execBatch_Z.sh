#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2
SUFFIX=$3
SELECTION=$4
PERIOD=$5

### Transfer files, prepare directory ###
TOPDIR=$PWD

# lpc
export SCRAM_ARCH=slc6_amd64_gcc530
export CMSSW_VERSION=CMSSW_10_2_13
source /cvmfs/cms.cern.ch/cmsset_default.sh

# nut3
#export SCRAM_ARCH=slc6_amd64_gcc530
#export CMSSW_VERSION=CMSSW_8_0_24_patch1
#source /software/tier3/osg/cmsset_default.sh 

# Setup CMSSW environment
#scram project CMSSW $CMSSW_VERSION
#cd $CMSSW_VERSION/src
#eval `scram runtime -sh`

#cp $TOPDIR/source.tar.gz .
tar -xzf source.tar.gz
cd $CMSSW_VERSION/src/
scramv1 b ProjectRename
cmsenv
cd BLT/BLTAnalysis/test
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt

echo $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
echo $PATH
pwd
cat input.txt
### Run the analyzer
#if [[ $DATANAME == DoubleMuon* ]]; then
#	echo $DATANAME ALL EVENTS 
#	zgAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#else 
#	echo $DATANAME REDUCED EVENTS
#	zgAnalyzer input.txt 300000 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#fi

echo "zgAnalyzer input.txt -1"
echo $DATANAME 
echo $SUFFIX 
echo $SELECTION 
echo $PERIOD 
echo $COUNT


#zgAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
zAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#zAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
### Copy output and cleanup ###
cp output_${DATANAME}_${COUNT}.root ${_CONDOR_SCRATCH_DIR}