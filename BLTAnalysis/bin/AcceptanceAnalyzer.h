// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef ACCEPTANCEANALYZER_HH
#define ACCEPTANCEANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <regex>


class AcceptanceAnalyzer: public BLTSelector {
public:
    AcceptanceAnalyzer();
    ~AcceptanceAnalyzer();

    void   Begin(TTree *tree);
    void   Init(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TH1D *hGenEvents, *hAccEvents;

    // Params and cuts
    std::unique_ptr<Parameters> params;




    /////////////////////////
    //                     //
    //   FIDUCIAL REGION   //
    //                     //
    /////////////////////////


//  // ATLAS 7 TeV (high mass ee)
//  const Int_t NBINS = 13;
//  Float_t xbins[NBINS+1] = {116, 130, 150, 170, 190, 210, 230, 250, 300, 400, 500, 700, 1000, 1500};
//  const double PT_MIN     = 25;
//  const double ETA_MAX    = 2.5;


//  // ATLAS 8 TeV (high mass ee & mumu)
//  Int_t nbins = 12;
//  Float_t xbins[12+1] = {116, 130, 150, 175, 200, 230, 260, 300, 380, 500, 700, 1000, 1500};
//  const double PT1_MIN    = 40,       PT2_MIN =   30;
//  const double ETA_MAX    = 2.5;


    // ATLAS 8 TeV (Z mass ee & mumu)
    Int_t nbins = 7;
    Float_t xbins[7+1] = {46, 66, 80, 91, 102, 116, 150, 200};
    const double PT1_MIN    = 20,       PT2_MIN =   20;
    const double ETA_MAX    = 2.4;
};


#endif  // ACCEPTANCEANALYZER_HH
