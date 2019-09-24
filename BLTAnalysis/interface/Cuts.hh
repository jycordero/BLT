#ifndef CUTS_HH
#define CUTS_HH

#include <string>
#include <vector>

class Cuts {
public:
    Cuts();
    ~Cuts() {}

    float leadJetPt, trailJetPt, leadMuPt, trailMuPt, leadElPt, trailElPt, gPtOverMass, gPt,
          zMassLow, zMassHigh, metLow, metHigh, zgMassLow, zgMassHigh;
    float EAMu[6];
    float EAEl[7];
    float EAPho[7][3];




    //////////////////////////////////////////////////////
    // Muon Cuts
    struct muIDCuts {
        //float pt;
	//float eta;
	float IsPF;
        float IsGLB;
        float IsTRK;
        float NormalizedChi2;
        float NumberOfValidMuonHits;
        float NumberOfMatchedStations;
        float NumberOfValidPixelHits;
        float TrackLayersWithMeasurement;
        float dxy;
        float dz;
        std::string cutName;
    } looseMuID, tightMuID, vetoMuID,TightMuID, muHZZID;

    struct muIsoCuts {
        float chIso04;
        float nhIso04;
        float phIso04;
        float relCombIso04;
        std::string cutName;
    } looseMuIso, tightMuIso;

    struct muDetIsoCuts {
        float hcalIso;
        float ecalIso;
        float trkIso;
        float relCombIso;
        std::string cutName;
    } looseMuDetIso, tightMuDetIso, amumuMuDetIso;
    /////////////////////////////////////////////////////






    ////////////////////////////////////////////////////
    // Electron Cuts
    struct elIDCuts {
        //broken into [0] barrel and [1] endcap
        float dEtaIn[2];
        float dPhiIn[2];
        float sigmaIetaIeta[2];
        float HadOverEm[2];
        float dxy[2];
        float dz[2];
        float fabsEPDiff[2];
        float ConversionMissHits[2];
        float dr03TkSumPt[2];
        float dr03EcalRecHitSumEt[2];
        float dr03HcalTowerSumEt[2];
        unsigned int   numberOfLostHits[2];
	float energyInverse[2];
        std::string cutName;
    } vetoElID, looseElID, mediumElID, tightElID, mvaPreElID;

    struct elIsoCuts {
        float chIso04;
        float nhIso04;
        float phIso04;
        float relCombIso04;
        float chIso;
        float nhIso;
        float phIso;
        float relCombIso;
        std::string cutName;
    } looseElIso, mediumElIso, tightElIso;

    struct elMVACuts {
        float mvaVal[6];
        float pt[2];
        float eta[3];
        int missHits[6];
        float sip[6];
        int conversionVeto[6];
        std::string cutName;
    } hzgMVAID, hzzMVAID, hwwMVAID, looseMVAElID, tightMVAElID, looseFallMVAElID;

    ////////////////////////////////////////////////////
    // Photon Cuts
    struct phIDCuts {
        //broken into [0] barrel and [1] endcap
        float PassedEleSafeVeto[2];
        float HadOverEm[2];
        float sigmaIetaIeta[2];
        float HcalIso[2];
        float TrkIso[2];
        float ChPfIso[2];
        std::string cutName;
    } loosePhID, mediumPhID, preSelPhID,preSelPhIDV2;

    struct phIsoCuts {
        float chIso[2];
        float nhIso[2];
        float phIso[2];
        float relCombIso[2];
        std::string cutName;
    } loosePhIso, mediumPhIso, preSelPhIso;

    struct phMVACuts {
        float mvaValCat1;
        float mvaValCat2;
        float mvaValCat3;
        float mvaValCat4;
        std::string cutName;
    } catPhMVAID, looseMVAPhID, tightMVAPhID;


    ////////////////////////////////////////////////////

    struct vbfJetIDCuts {
        float betaStarC[2];
        float dR2Mean[4];
        std::string cutName;
    } vbfJetID;




    ////////////////////////////////////////////////////
    // Jet Cuts
    struct jetIDCuts {
        unsigned int NumConst[4];
	unsigned int NumNeuPart[4];
        unsigned int CHM[4];
        float NHF[4];
        float NEMF[4];
        float CHF[4];
        float CEMF[4];
        float MUF[4];
        float CSV[4];
        std::string cutName[4];
    } looseJetID, tightJetID, bJetID;
};

#endif  // CUTS_HH
