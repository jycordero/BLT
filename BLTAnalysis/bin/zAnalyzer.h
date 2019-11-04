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


#ifndef ZGANALYZER_HH
#define ZGANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"
//#include "BLT/BLTAnalysis/interface/RoccoR_2016ReReco.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

// C++ headers
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>
#define MZ 91.1876
#define WZ 2.4952

class zAnalyzer: public BLTSelector {
public:
    zAnalyzer();
    ~zAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    TH1D *hTotalEventsGen;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    RoccoR *muonCorr;
    //RoccoR_2016ReReco *muonCorr;
    TRandom3 *rng;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
   
    Float_t genWeight,eventWeight,puWeight;
 
    // event data
    UInt_t runNumber, lumiSection, nPV, nPartons;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t nPU,Rho;
    Float_t xPV, yPV, zPV;
    UInt_t nJets, nCentralJets, nFwdJets, nBJets, nMuons, nElectrons, nTaus, nPhotons;

    // generator level data
    Float_t genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi;
    Float_t genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi;
    Float_t genPhotonPt, genPhotonEta, genPhotonPhi;
    Int_t genLeptonOneId, genLeptonTwoId;
    Bool_t genPhotonFHPFS, genPhotonIPFS;
   
    // physics object Lorentz vectors
    Float_t leptonOnePt, leptonOneEta, leptonOnePhi;
    Float_t leptonOnePassPt, leptonOnePassEta, leptonOnePassPhi;

    Float_t leptonOnePtKin, leptonOnePassPtKin;
    Float_t leptonOnePtKinJames, leptonOnePassPtKinJames;

    // Additional lepton data
    Float_t leptonOneIso, leptonOnePassIso;
    Int_t leptonOneMother, leptonOnePassMother;
    Int_t leptonOneFlavor, leptonOnePassFlavor;
    Float_t leptonOneD0, leptonOnePassD0;
    Float_t leptonOneDZ, leptonOnePassDZ;
    Float_t leptonOneCharge, leptonOnePassCharge;
    Float_t leptonOneTag, leptonOnePassTag;
    Float_t leptonOneRecoWeight, leptonOnePassRecoWeight;

    Bool_t leptonOneECALDriven, leptonOnePassECALDriven;

    // physics object Lorentz vectors
    Float_t leptonTwoPt, leptonTwoEta, leptonTwoPhi;
    Float_t leptonTwoPassPt, leptonTwoPassEta, leptonTwoPassPhi;

    Float_t leptonTwoPtKin, leptonTwoPassPtKin;
    Float_t leptonTwoPtKinJames, leptonTwoPassPtKinJames;

    // Additional lepton data
    Float_t leptonTwoIso, leptonTwoPassIso;
    Int_t leptonTwoMother, leptonTwoPassMother;
    Int_t leptonTwoFlavor, leptonTwoPassFlavor;
    Float_t leptonTwoD0, leptonTwoPassD0;
    Float_t leptonTwoDZ, leptonTwoPassDZ;
    Float_t leptonTwoCharge, leptonTwoPassCharge;
    Float_t leptonTwoTag, leptonTwoPassTag;
    Float_t leptonTwoRecoWeight, leptonTwoPassRecoWeight;

    Bool_t leptonTwoECALDriven, leptonTwoPassECALDriven;
    // physics object Lorentz vectors
    Float_t leptonProbePt, leptonProbeEta, leptonProbePhi;
    Float_t leptonProbeFailPt, leptonProbeFailEta, leptonProbeFailPhi;

    Float_t leptonProbePtKin, leptonProbeFailPtKin;
    Float_t leptonProbePtKinJames, leptonProbeFailPtKinJames;

    // Additional lepton data
    Float_t leptonProbeFailIso, leptonProbeFailFailIso;
    Int_t leptonProbeFailMother, leptonProbeFailFailMother;
    Int_t leptonProbeFailFlavor, leptonProbeFailFailFlavor;
    Float_t leptonProbeFailD0, leptonProbeFailFailD0;
    Float_t leptonProbeFailDZ, leptonProbeFailFailDZ;
    Float_t leptonProbeFailCharge, leptonProbeFailFailCharge;
    Float_t leptonProbeFailTag, leptonProbeFailFailTag;
    Float_t leptonProbeFailRecoWeight, leptonProbeFailFailRecoWeight;

    Bool_t leptonProbeFailECALDriven, leptonProbeFailFailECALDriven;
    
    // physics object Lorentz vectors
    Float_t leptonProbePassPt, leptonProbePassEta, leptonProbePassPhi;
    Float_t leptonProbePassPassPt, leptonProbePassPassEta, leptonProbePassPassPhi;

    Float_t leptonProbePassPtKin, leptonProbePassPassPtKin;
    Float_t leptonProbePassPtKinJames, leptonProbePassPassPtKinJames;

    // Additional lepton data
    Float_t leptonProbePassIso, leptonProbePassPassIso;
    Int_t leptonProbePassMother, leptonProbePassPassMother;
    Int_t leptonProbePassFlavor, leptonProbePassPassFlavor;
    Float_t leptonProbePassD0, leptonProbePassPassD0;
    Float_t leptonProbePassDZ, leptonProbePassPassDZ;
    Float_t leptonProbePassCharge, leptonProbePassPassCharge;
    Float_t leptonProbePassTag, leptonProbePassPassTag;
    Float_t leptonProbePassRecoWeight, leptonProbePassPassRecoWeight;

    Bool_t leptonProbePassECALDriven, leptonProbePassPassECALDriven;

    // tau data
    Int_t tauDecayMode;
    Float_t tauMVA; 
    //UInt_t tauPhotonMult, tauChHadMult;

    // photon data
    Float_t photonOnePt, photonOneEta, photonOnePhi;
    Float_t photonOneR9;
    Float_t photonOneMVA;
    Float_t photonOneERes;
    Float_t photonOneSieie, photonOneHoverE, photonOneIneu, photonOneIph, photonOneIch;
    Bool_t passElectronVeto;

    Float_t photonOneSieip;      
    Float_t photonOneSipip;     
    Float_t photonOneSrr;       
    Float_t photonOneE2x2;      
    Float_t photonOneE5x5;      
    Float_t photonOneScEtaWidth; 
    Float_t photonOneScPhiWidth; 
    Float_t photonOneScRawE; 
    Float_t photonOnePreShowerE; 
    Float_t photonOneScBrem; 



    //Int_t genOneId, genTwoId, genOneMother, genTwoMother, genCategory;
    //TLorentzVector genOneP4, genTwoP4;
    //Bool_t fromHardProcessFinalState, isPromptFinalState, hasPhotonMatch;
    Bool_t vetoDY, genIsoPass;
    Bool_t TagFromZ, ProbeFromZ;
    Bool_t ProbePass;

    // dilepton data
    Float_t dileptonPt, dileptonEta, dileptonPhi, dileptonM;
    Float_t dileptonDEta, dileptonDPhi, dileptonDR;
    Float_t dileptonMKin;
    Float_t dileptonMKinJames;
    
    // dileptonProbe data
    Float_t dileptonProbePt, dileptonProbeEta, dileptonProbePhi, dileptonProbeM;
    Float_t dileptonProbeDEta, dileptonProbeDPhi, dileptonProbeDR;
    Float_t dileptonProbeMKin;
    Float_t dileptonProbeMKinJames;

    // dileptonProbePass data
    Float_t dileptonProbePassPt, dileptonProbePassEta, dileptonProbePassPhi, dileptonProbePassM;
    Float_t dileptonProbePassDEta, dileptonProbePassDPhi, dileptonProbePassDR;
    Float_t dileptonProbePassMKin;
    Float_t dileptonProbePassMKinJames;

    // dileptonProbeFail data
    Float_t dileptonProbeFailPt, dileptonProbeFailEta, dileptonProbeFailPhi, dileptonProbeFailM;
    Float_t dileptonProbeFailDEta, dileptonProbeFailDPhi, dileptonProbeFailDR;
    Float_t dileptonProbeFailMKin;
    Float_t dileptonProbeFailMKinJames;
    // dilepton vertex data
    //Float_t dileptonVertexOneX, dileptonVertexOneY, dileptonVertexOneZ;
    //Float_t dileptonVertexTwoX, dileptonVertexTwoY, dileptonVertexTwoZ;
    //Float_t dileptonVertexOneXErr, dileptonVertexOneYErr, dileptonVertexOneZErr;
    //Float_t dileptonVertexTwoXErr, dileptonVertexTwoYErr, dileptonVertexTwoZErr;
    //Float_t dileptonVertexChi2One, dileptonVertexDOFOne;
    //Float_t dileptonVertexChi2Two, dileptonVertexDOFTwo;
    
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetPhotonIsolation(const baconhep::TPhoton*, float);
    float GetGenIsolation(const TGenParticle*);

    void EvalMuonEnergyResolution(std::map<string, float>, std::map<string, int>, float&, float&, float&, float&, float&, float&);
    void EvalElectronEnergyResolution(std::map<string, float>, float&, float&, float&, float&, float&, float&);
    void find_optimized(double*, double&, double&);

    bool SignalRegionPass(const baconhep::TPhoton*);

    TH2F *PhotonProbe;
    TH2F *PhotonProbePass;
    TH2F *PhotonProbeFail;

    //ClassDef(zAnalyzer,0);
};

double GetDoubleSidedCB(double x, double mean, double sigma, double alphaL,
                                          double powerL, double alphaR, double powerR)
{
   // Returns value of double-sided Crystal Ball function for given parameters.
   //
   // Negative powerL and/or powerR means an infinite value.
    double a = (x - mean) / sigma;
    // left power-law or exponential tail
   if (a < -alphaL) {
      if (powerL < 0) // infinite powerL
         return exp(0.5 * alphaL*alphaL + alphaL*a);
       double b = powerL/alphaL;
      return exp(-0.5 * alphaL*alphaL) * pow(b/(b - alphaL - a), powerL);
   }
    // Gaussian core
   if (a <= alphaR)
      return exp(-0.5*a*a);
    // right exponential tail
   if (powerR < 0) // infinite powerR
      return exp(0.5 * alphaR*alphaR - alphaR*a);
    // right power-law tail
   double b = powerR/alphaR;
   return exp(-0.5 * alphaR*alphaR) * pow(b/(b - alphaR + a), powerR);
}
 Double_t NegativeProbability(Double_t* x, Double_t* p)
{
   // Negative probability density function to minimize.
    // lepton energies to optimize
   double e1 = x[0];
   double e2 = x[1];
    double m02 = p[0];       // lepton mass squared
   double cosAlpha = p[1];  // cosine of the angle between lepton directions
    // two-lepton invariant mass squared
   double mz2 = 2 * (e1*e2 + m02
                     - sqrt(e1*e1 - m02) * sqrt(e2*e2 - m02) * cosAlpha);
    // relativistic Breit-Wigner
   double a = mz2 - MZ*MZ;
   double probZ = 1/(a*a + mz2*mz2 * WZ*WZ/(MZ*MZ));
    // p[2] = first reconstructed energy
   // p[3] = second reconstructed energy
   // p[4], ... = parameters of energy resolution functions
    double prob1 = GetDoubleSidedCB(e1/p[2], p[4], p[5], p[6], p[7], p[8], p[9]);
   double prob2 = GetDoubleSidedCB(e2/p[3], p[10], p[11], p[12], p[13], p[14], p[15]);
    return -probZ * prob1 * prob2;
}

#endif  // HZGANALYZER_HH
