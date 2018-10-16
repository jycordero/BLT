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

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

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
    TTree *outTree;

    // Params and cuts
    std::unique_ptr<Parameters>         params;

    // Histograms
    TH1D *hPhaseSpaceEvents, *hFiducialEvents;



    //--- SELECTION CUTS ---//

    // Phase space requirements
    Float_t M_MIN = 80,         M_MAX = 100,        MLL_MIN = 4;

    // Fiducial requirements
    // (decided they are the same for both flavors...)
    Float_t PT1_MIN = 20,    PT2_MIN = 10,      PT_MIN = 5,     ETA_MAX = 2.5;



    //--- BRANCHES ---//
    
    // Event
    Int_t       runNumber,  lumiSection;
    Long64_t    evtNumber;
    Float_t     genWeight;
    Bool_t      isFiducial;

    // Counters
    UShort_t    nGenMuons,  nGenElectrons,  nGenLeptons;

    // Gen particles
    TClonesArray *genMuonsP4 = new TClonesArray("TLorentzVector"), &genMuonsP4ptr = *genMuonsP4;
    std::vector<Short_t>    genMuonsQ,      genMuonStatus;

    TClonesArray *genElectronsP4 = new TClonesArray("TLorentzVector"), &genElectronsP4ptr = *genElectronsP4;
    std::vector<Short_t>    genElectronsQ,  genElectronStatus;



    //--- HELPER FUNCTIONS ---//

//  TLorentzVector GetP4Sum(const std::vector<TLorentzVector>&);

//  bool SortDecPt(const std::pair<TLorentzVector, Int_t>&, const std::pair<TLorentzVector, Int_t>&);



    //ClassDef(AcceptanceAnalyzer,0);
};


#endif  // ACCEPTANCEANALYZER_HH
