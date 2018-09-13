BLT
===

BLT is a framework for analyzing bacon ntuples originally developed by @jiafulow.

This branch is designed for the BF(Z -> 4l) measurement over the full 2016 dataset (Moriond 2017).

Setup
=====

This branch uses CMSSW_8_0_26_patch1 on cmslpc-sl6:

### Build and source environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
```

### Checkout dependencies

Depends on NWUHEP/BaconAna topic_jbueghly branch:

```
git clone -b topic_jbueghly git@github.com:NWUHEP/BaconAna
```

### Checkout and compile BLT code

```
git clone -b lacey_13TeV git@github.com:NWUHEP/BLT.git
scram b -j 12
```

## Running the analyzer

7 input arguments are mandatory: [input file] [no of events] [dataset] [datasetgroup] [selection] [period] [jobid]
[input file]

Currently [selection] = "emu"

```
cd BLT/BLTAnalysis/test
MultileptonAnalyzer input/DYJetsToLL_M-50.txt 1000 DYJetsToLL_M-50 DYJetsToLL emu 2016 0
```

## Running a BLT analyzer with condor

Must be on cmslpc

```
cd BLT/BLTAnalysis/test
./batch_lpc_cfg.py
```

