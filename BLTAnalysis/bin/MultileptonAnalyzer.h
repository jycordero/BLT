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

    std::vector<string> triggerNames;


    // Branches in the output file

    // event data
    UInt_t runNumber, lumiSection;
    UShort_t nPV, nPartons, nJets;
    ULong64_t evtNumber;
    Float_t eventWeight, nPU, PUWeight;
    Float_t met, metPhi;
    Bool_t passTrigger, passHLT_IsoMu24, passHLT_IsoTkMu24, passHLT_Ele27_WPTight_Gsf;

    // leptons
    TClonesArray *muonsP4 = new TClonesArray("TLorentzVector");
    TClonesArray &muonsP4ptr = *muonsP4;
    TClonesArray *electronsP4 = new TClonesArray("TLorentzVector");
    TClonesArray &electronsP4ptr = *electronsP4;
    std::vector<Short_t> muonsQ, electronsQ;

    // muon info
    std::vector<Bool_t> muonIsPF, muonIsGLB, muonPassStdCuts, muonPassTrigger; 
    std::vector<Float_t> muonIDEff, muonTightIsoEff, muonLooseIsoEff;
    std::vector<Float_t> muonTriggerEffData, muonTriggerEffMC;
    std::vector<Float_t> muonCombIso, muonsTrkIso;
    std::vector<Float_t> muonSF, muonMuNChi2, muonD0, muonDz;
    std::vector<UShort_t> muonNMatchStn, muonNPixHits, muonNTkLayers;
    std::vector<UShort_t> muonNValidHits;

    // electron info
    std::vector<Bool_t> electronIsConv, electronPassID, electronPassIso;
    std::vector<Bool_t> electronPassStdCuts, electronPassTrigger;
    std::vector<Float_t> electronRecoEff, electronTriggerEffData, electronTriggerEffMC;
    std::vector<Float_t> electronCombIso, electronsTrkIso, electronEnergyInv;
    std::vector<Float_t> electronSF, electronScEta, electronD0, electronDz, electronSieie;
    std::vector<Float_t> electronHOverE, electronDEtaIn, electronDPhiIn;
    std::vector<UShort_t> electronNMissHits;

    // gen-level particles
    TClonesArray *genMuonsP4 = new TClonesArray("TLorentzVector");
    TClonesArray &genMuonsP4ptr = *genMuonsP4;
    TClonesArray *genElectronsP4 = new TClonesArray("TLorentzVector");
    TClonesArray &genElectronsP4ptr = *genElectronsP4;
    std::vector<Short_t> genMuonsQ, genElectronsQ;
    std::vector<Short_t> genIntermID;
    std::vector<Float_t> genIntermMass;

    // counters
    UShort_t nMuons, nElectrons, nLeptons;
    UShort_t nStdMuons, nStdElectrons, nStdLeptons;
    UShort_t nGenMuons, nGenElectrons, nGenLeptons;


    // Helper functions 
    float MetKluge(float);
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);

    //ClassDef(MultileptonAnalyzer,0);
};


#endif  // MULTILEPTONANALYZER_HH
