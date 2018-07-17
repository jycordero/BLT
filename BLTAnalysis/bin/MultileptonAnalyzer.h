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

//  std::vector<string>                 triggerNames;


    // Branches in the output file
    // event data
    Int_t runNumber, lumiSection;
    Long64_t evtNumber;
    UShort_t nPV;
    Float_t eventWeight, PUWeight, PUVar, nPU;
    UShort_t nPartons;
    Float_t met, metPhi;
    Bool_t evtMuonTriggered, evtElectronTriggered;

    // counters
    UShort_t nMuons, nElectrons, nLeptons;
    UShort_t nLooseMuons, nLooseElectrons, nLooseLeptons;
    UShort_t nVetoElectrons, nMediumElectrons;
    UShort_t nTightMuons, nTightElectrons, nTightLeptons;
    UShort_t nHZZMuons, nHZZElectrons, nHZZLeptons;
    UShort_t nGenMuons, nGenElectrons, nGenLeptons;

    // muon info
    TClonesArray *muonsP4 = new TClonesArray("TLorentzVector"), &muonsP4ptr = *muonsP4;
    std::vector<Short_t> muonsQ;
    std::vector<Float_t> muonCombIso, muonsTrkIso;

    std::vector<Bool_t> muonIsPF, muonIsGlobal, muonIsTracker;
    std::vector<Bool_t> muonIsLoose, muonIsTight, muonIsHZZ, muonTriggered; 

    std::vector<Float_t> muonSF;
    std::vector<Float_t> muonLooseIDWeight, muonLooseIDVar, muonTightIDWeight, muonTightIDVar, muonHZZIDWeight, muonHZZIDVar;
    std::vector<Float_t> muonLooseIsoWeight, muonLooseIsoVar, muonTightIsoWeight, muonTightIsoVar;
    std::vector<Float_t> muonTriggerEffData, muonTriggerEffMC, muonTriggerErrData, muonTriggerErrMC;

    std::vector<Float_t> muonD0, muonDz, muonSIP3d;
    std::vector<Float_t> muonMuNChi2, muonPtErr;
    std::vector<UShort_t> muonNMatchStn, muonNPixHits, muonNTkLayers;
    std::vector<Short_t> muonBestTrackType, muonNValidHits;

    // electron info
    TClonesArray *electronsP4 = new TClonesArray("TLorentzVector"), &electronsP4ptr = *electronsP4;
    std::vector<Short_t> electronsQ;
    std::vector<Float_t> electronCombIso, electronsTrkIso;

    std::vector<Bool_t> electronPassVetoIso, electronPassLooseIso, electronPassMediumIso, electronPassTightIso;
    std::vector<Bool_t> electronIsVeto, electronIsLoose, electronIsMedium, electronIsTight, electronIsHZZ, electronTriggered;

    std::vector<Float_t> electronSF;
    std::vector<Float_t> electronRecoWeight, electronRecoVar, electronHZZRecoWeight, electronHZZRecoVar;
    std::vector<Float_t> electronTriggerEffData, electronTriggerEffMC, electronTriggerErrData, electronTriggerErrMC;

    std::vector<Float_t> electronD0, electronDz, electronSIP3d;
    std::vector<Float_t> electronScEta, electronSieie, electronEnergyInv;
    std::vector<Float_t> electronHOverE, electronDEtaIn, electronDPhiIn;
    std::vector<UShort_t> electronNMissHits;
    std::vector<Bool_t> electronIsConv;

    // gen-level particles
    TClonesArray *genMuonsP4 = new TClonesArray("TLorentzVector"), &genMuonsP4ptr = *genMuonsP4;
    std::vector<Short_t> genMuonsQ, genMuonStatus;

    TClonesArray *genElectronsP4 = new TClonesArray("TLorentzVector"), &genElectronsP4ptr = *genElectronsP4;
    std::vector<Short_t> genElectronsQ, genElectronStatus;

//  std::vector<Short_t> genIntermID;
//  std::vector<Float_t> genIntermMass;


    // Helper functions 
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
