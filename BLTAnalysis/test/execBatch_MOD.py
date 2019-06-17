#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2

### Specify addtional arguments here ####
SUFFIX=$3
SELECTION=$4
PERIOD=$5

### Transfer files, prepare directory ###
TOPDIR=$PWD

#export SCRAM_ARCH=slc6_amd64_gcc630
#export SCRAM_ARCH=slc6_amd64_gcc491
#export CMSSW_VERSION=CMSSW_7_4_14
export SCRAM_ARCH=slc6_amd64_gcc530
#export CMSSW_VERSION=CMSSW_8_0_26_patch1
export CMSSW_VERSION=CMSSW_8_0_20
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh

########################################
echo "--------- 1 ---------"
tar -xzf source.tar.gz
echo "~~~ PATH"
echo $PWD
echo "~~~ LS"
ls

echo 
echo

echo "-------- 2 ---------"
# temporary fix
mv $CMSSW_VERSION tmp
echo "~~~ PATH"
echo $PWD
echo "~~~ LS"
ls
echo "~~~ tmp"
ls tmp
echo 
echo


echo "-------- 3 ---------"
scram project CMSSW $CMSSW_VERSION
#cp -r tmp/src/* $CMSSW_VERSION/src
mv tmp/$CMSSW_VERSION/src/* $CMSSW_VERSION/src
echo "~~~ PATH"
echo $PWD
echo "~~~ tmp/cmssw/src"
ls tmp/$CMSSW_VERSION/src
echo "~~~ cmssw/src"
ls $CMSSW_VERSION/src

echo 
echo

echo "------------ 4 -------------"
cd $CMSSW_VERSION/src
echo "~~~ PATH"
echo $PWD
echo "~~~ LS"
ls
echo 
echo


echo "------------ 5 -------------"
cmsenv
#scram b -j8
cd BLT/BLTAnalysis/test
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt
echo "~~~ PATH"
echo $PWD
echo "~~~ LS"
ls
echo 
echo


echo "----------- 6 ------------"
cd $CMSSW_VERSION/src
echo "~~~ PATH"
echo $PWD
echo "~~~ LS"
ls
echo
echo "~~~ BEGGING RUNNING cmsRUN ~~~~~~"
echo 
echo
### Run the analyzer
#DimuonAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
cmsRun makingBacon_MC_25ns_MINIAOD.py inputFiles_load=input.txt outputFile=output_${SUFFIX}_${COUNT}.root

### Copy output and cleanup ###
echo "------------  6--------------"
cp output_${SUFFIX}_${COUNT}.root ${_CONDOR_SCRATCH_DIR}
echo "~~~ PATH"
echo $PWD
echo
echo "~~~ LS"
ls
echo 
echo
