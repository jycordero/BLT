BLT
===

BLT is a framework for analyzing bacon ntuples originally developed by @jiafulow.

This branch is designed for the BF(Z -> 4l) measurement over the 8 TeV data.

Setup
=====

This branch uses CMSSW_7_4_12 on nu tier3:

### Build and source environment

```
PATH=$PATH:$HOME/bin

export PATH
export OSG_APP=/software/tier3/osg
export SCRAM_ARCH=slc6_amd64_gcc481
source /software/tier3/osg/cmsset_default.sh
cmsrel CMSSW_7_4_12
cd CMSSW_7_4_12/src
cmsenv
```

### Checkout dependencies

Depends on NWUHEP/BaconAna tag 04:

```
git clone -b 04 git@github.com:NWUHEP/BaconAna
```

### Checkout and compile BLT code

```
git clone -b lacey_8TeV git@github.com:NWUHEP/BLT.git
scram b -j 12
```

## Running the analyzer

7 input arguments are mandatory: [input file] [no of events] [dataset] [datasetgroup] [selection] [period] [jobid]
[input file]

Currently [selection] = "single" ("double") for single-lepton (double-lepton) triggers

```
cd BLT/BLTAnalysis/test
MultileptonAnalyzer input/DYJetsToLL_M-50.txt 1000 DYJetsToLL_M-50 DYJetsToLL single 2012 0
```

## Running a BLT analyzer with condor

Must be on tier3:

```
cd BLT/BLTAnalysis/test
./batch_cfg.py
```
