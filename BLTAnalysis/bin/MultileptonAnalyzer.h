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


#ifndef MULTILEPTONANALYZER_HH
#define MULTILEPTONANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"
#include "BLT/BLTAnalysis/interface/ElectronCorrector.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

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
#include <iterator>
#include <regex>


class MultileptonAnalyzer: public BLTSelector {
public:
    MultileptonAnalyzer();
    ~MultileptonAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    RoccoR *muonCorr;
    TRandom3 *rng;

    // electron scale and smear corrector (trash)
    // (from topic_wbranch)
    EnergyScaleCorrection *electronScaler;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string>                 muonTriggerNames, electronTriggerNames;



    //--- BRANCHES ---//

    // Event
    Int_t runNumber, lumiSection;
    Long64_t evtNumber;
    UShort_t nPV;
    Float_t genWeight, PUWeight, nPU;
    UShort_t nPartons;
    Float_t met, metPhi;
    Bool_t evtMuonTriggered, evtElectronTriggered;


    // Counters
    UShort_t nMuons, nElectrons, nLeptons;
    UShort_t nLooseMuons, nMVAElectrons;
    UShort_t nHZZMuons, nHZZElectrons, nHZZLeptons;
    UShort_t nGenMuons, nGenElectrons, nGenLeptons;

    
    // Muons
    TClonesArray *muonsP4 = new TClonesArray("TLorentzVector"), &muonsP4ptr = *muonsP4;
    std::vector<Short_t>    muonsQ;
    std::vector<Bool_t>     muonFiredLeg1, muonFiredLeg2;

    std::vector<Bool_t>     muonIsLoose, muonIsHZZ; 
    std::vector<Float_t>    muonEnergySF, muonHZZIDSF;
    std::vector<Float_t>    muonTrigEffLeg1Data, muonTrigEffLeg1MC, muonTrigEffLeg2Data, muonTrigEffLeg2MC;

    std::vector<Float_t>    muonCombIso, muonsTrkIso, muonD0, muonDz, muonSIP3d, muonMuNChi2, muonPtErr;
    std::vector<UShort_t>   muonNMatchStn, muonNPixHits, muonNTkLayers;
    std::vector<Bool_t>     muonIsGlobal, muonIsTracker;
    std::vector<Short_t>    muonBestTrackType, muonNValidHits;


    // Electrons
    TClonesArray *electronsP4 = new TClonesArray("TLorentzVector"), &electronsP4ptr = *electronsP4;
    std::vector<Short_t>    electronsQ;
    std::vector<Bool_t>     electronFiredLeg1, electronFiredLeg2;

    std::vector<Bool_t>     electronIsMVA, electronIsHZZ;
    std::vector<Float_t>    electronEnergySF, electronHZZIDSF;
    std::vector<Float_t>    electronTrigEffLeg1Data, electronTrigEffLeg1MC, electronTrigEffLeg2Data, electronTrigEffLeg2MC;

    std::vector<Float_t>    electronCombIso, electronsTrkIso, electronD0, electronDz, electronSIP3d, electronScEta;
    std::vector<Float_t>    electronMVA, electronSieie, electronEnergyInv, electronHOverE, electronDEtaIn, electronDPhiIn;
    std::vector<UShort_t>   electronNMissHits;
    std::vector<Bool_t>     electronIsConv;


    // Gen particles
    TClonesArray *genMuonsP4 = new TClonesArray("TLorentzVector"), &genMuonsP4ptr = *genMuonsP4;
    std::vector<Short_t> genMuonsQ, genMuonStatus;

    TClonesArray *genElectronsP4 = new TClonesArray("TLorentzVector"), &genElectronsP4ptr = *genElectronsP4;
    std::vector<Short_t> genElectronsQ, genElectronStatus;



    //--- HELPER FUNCTIONS ---//

    float MetKluge(float);

    float GetMuonIsolation(const baconhep::TMuon*);
    float GetRochesterCorrection(const baconhep::TMuon*, RoccoR*, TRandom3*, bool);
    bool PassMuonTightID(const baconhep::TMuon*);
    bool PassMuonHZZTightID(const baconhep::TMuon*);

    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetElectronPtSF(baconhep::TElectron*, EnergyScaleCorrection*, TRandom3*, int);
    bool PassElectronGoodID(const baconhep::TElectron*, std::unique_ptr<ParticleSelector>&, std::unique_ptr<Cuts>&);
    bool PassElectronHZZTightID(const baconhep::TElectron*, std::unique_ptr<ParticleSelector>&, std::unique_ptr<Cuts>&, float);

    //ClassDef(MultileptonAnalyzer,0);
};


#endif  // MULTILEPTONANALYZER_HH
