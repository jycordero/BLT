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
#include <utility>
#include <iterator>
#include <regex>


class AcceptanceAnalyzer: public BLTSelector {
public:
    AcceptanceAnalyzer();
    ~AcceptanceAnalyzer();

    void    Begin(TTree *tree);
    void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();
    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile *outFile;

    // Params and cuts
    std::unique_ptr<Parameters>         params;



    //--- SELECTION CUTS ---//

    // Phase space requirements
    Float_t M_MIN = 80,         M_MAX = 100,        MLL_MIN = 4;

    // Fiducial requirements
    // (decided they are the same for both flavors...)
    Float_t PT1_MIN = 20,    PT2_MIN = 10,      PT_MIN = 5,     ETA_MAX = 2.5;



    //--- HELPER FUNCTIONS ---//

    TLorentzVector GetP4Sum(vector<TLorentzVector);

    bool SortDecPt(const pair<TLorentzVector, Int_t> &i_, const pair<TLorentzVector, Int_t> &j_);



    //ClassDef(AcceptanceAnalyzer,0);
};


#endif  // ACCEPTANCEANALYZER_HH
