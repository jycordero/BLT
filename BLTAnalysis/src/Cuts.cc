#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <iostream>
#include <stdexcept>


Cuts::Cuts() {
    leadJetPt   = 30;
    trailJetPt  = 30;
    leadMuPt    = 20;
    trailMuPt   = 10;
    leadElPt    = 20;
    trailElPt   = 10;
    gPtOverMass = 15./100.;
    gPt         = 15.;
    zMassLow    = 50;
    zMassHigh   = 999999;
    metLow      = -999999;
    metHigh     = 999999;
    zgMassLow   = 100;
    zgMassHigh  = 190;

    //combined rel ISO, 2012 Data, 0.5 GeV
    EAMu[0] =   0.674; //         eta < 1.0
    EAMu[1] =   0.565; // 1.0   < eta < 1.5
    EAMu[2] =   0.442; // 1.5   < eta < 2.0
    EAMu[3] =   0.515; // 2.0   < eta < 2.2
    EAMu[4] =   0.821; // 2.2   < eta < 2.3
    EAMu[5] =   0.660; // 2.3   < eta < 2.4

    EAEl[0] =   0.13; //         eta < 1.0
    EAEl[1] =   0.14; // 1.0   < eta < 1.5
    EAEl[2] =   0.07; // 1.5   < eta < 2.0
    EAEl[3] =   0.09; // 2.0   < eta < 2.2
    EAEl[4] =   0.11; // 2.2   < eta < 2.3
    EAEl[5] =   0.11; // 2.3   < eta < 2.4
    EAEl[6] =   0.14; // 2.4   < eta

    //   ch      nh       ph
    float EAPhoTemp[7][3] = {
        {0.012,  0.030,   0.148}, //         eta < 1.0
        {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
        {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
        {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
        {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
        {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
        {0.012,  0.072,   0.266}  // 2.4   < eta
    };

    for (unsigned int i =0; i<7; i++)
        for (unsigned int j =0; j<3; j++)
            EAPho[i][j] = EAPhoTemp[i][j];

    vetoElID.cutName                      = "vetoElID";
    vetoElID.dEtaIn[0]                    = 0.00749;
    vetoElID.dPhiIn[0]                    = 0.228;
    vetoElID.sigmaIetaIeta[0]             = 0.0115;
    vetoElID.HadOverEm[0]                 = 0.356;
    vetoElID.dxy[0]                       = 1;
    vetoElID.dz[0]                        = 1;
    vetoElID.fabsEPDiff[0]                = 99999;
    vetoElID.ConversionMissHits[0]        = 0;
    vetoElID.numberOfLostHits[0]	  = 2;
    vetoElID.energyInverse[0]		  = 0.299;

    vetoElID.dEtaIn[1]                    = 0.00895;
    vetoElID.dPhiIn[1]                    = 0.213;
    vetoElID.sigmaIetaIeta[1]             = 0.037;
    vetoElID.HadOverEm[1]                 = 0.211;
    vetoElID.dxy[1]                       = 1;
    vetoElID.dz[1]                        = 1;
    vetoElID.fabsEPDiff[1]                = 999999;
    vetoElID.ConversionMissHits[1]        = 0;
    vetoElID.numberOfLostHits[1]	  = 3;
    vetoElID.energyInverse[1]		  = 0.15;

    // tight electorn ID
    tightElID.cutName                    = "tightElID";
    tightElID.dEtaIn[0]                  = 0.00308;
    tightElID.dPhiIn[0]                  = 0.0816;
    tightElID.sigmaIetaIeta[0]           = 0.00998;
    tightElID.HadOverEm[0]               = 0.0414;
    tightElID.dxy[0]                     = 0.05;
    tightElID.dz[0]                      = 0.1;
    tightElID.fabsEPDiff[0]              = 0.0588;
    tightElID.ConversionMissHits[0]      = 1;
    tightElID.numberOfLostHits[0]	 = 1;
    tightElID.energyInverse[0]		 = 0.0129;

    tightElID.dEtaIn[1]                  = 0.00605;
    tightElID.dPhiIn[1]                  = 0.0394;
    tightElID.sigmaIetaIeta[1]           = 0.0292;
    tightElID.HadOverEm[1]               = 0.0641;
    tightElID.dxy[1]                     = 0.1;
    tightElID.dz[1]                      = 0.2;
    tightElID.fabsEPDiff[1]              = 0.0571;
    tightElID.ConversionMissHits[1]      = 1;
    tightElID.numberOfLostHits[1]	 = 1;
    tightElID.energyInverse[1]		 = 0.0129;

    // medium electron ID
    mediumElID.cutName                    = "mediumElID";
    mediumElID.dEtaIn[0]                  = 0.00311;
    mediumElID.dPhiIn[0]                  = 0.103;
    mediumElID.sigmaIetaIeta[0]           = 0.00998;
    mediumElID.HadOverEm[0]               = 0.253;
    mediumElID.dxy[0]                     = 1;
    mediumElID.dz[0]                      = 1;
    mediumElID.fabsEPDiff[0]              = 0.05;
    mediumElID.ConversionMissHits[0]      = 0;
    mediumElID.numberOfLostHits[0]	  = 1;
    mediumElID.energyInverse[0]		  = 0.134;

    mediumElID.dEtaIn[1]                  = 0.00609;
    mediumElID.dPhiIn[1]                  = 0.045;
    mediumElID.sigmaIetaIeta[1]           = 0.0298;
    mediumElID.HadOverEm[1]               = 0.0878;
    mediumElID.dxy[1]                     = 1;
    mediumElID.dz[1]                      = 1;
    mediumElID.fabsEPDiff[1]              = 0.05;
    mediumElID.ConversionMissHits[1]      = 0;
    mediumElID.numberOfLostHits[1]	  = 1;
    mediumElID.energyInverse[1]		  = 0.13;
    // loose electron ID
    looseElID.cutName                     = "looseElID";
    looseElID.dEtaIn[0]                   = 0.00477;
    looseElID.dPhiIn[0]                   = 0.222;
    looseElID.sigmaIetaIeta[0]            = 0.011;
    looseElID.HadOverEm[0]                = 0.298;
    looseElID.dxy[0]                      = 1;
    looseElID.dz[0]                       = 1;
    looseElID.fabsEPDiff[0]               = 0.05;
    looseElID.ConversionMissHits[0]       = 0;
    looseElID.numberOfLostHits[0]	  = 1;
    looseElID.energyInverse[0]		  = 0.241;

    looseElID.dEtaIn[1]                   = 0.00868;
    looseElID.dPhiIn[1]                   = 0.213;
    looseElID.sigmaIetaIeta[1]            = 0.0314;
    looseElID.HadOverEm[1]                = 0.101;
    looseElID.dxy[1]                      = 1;
    looseElID.dz[1]                       = 1;
    looseElID.fabsEPDiff[1]               = 0.05;
    looseElID.ConversionMissHits[1]       = 0;
    looseElID.numberOfLostHits[1]	  = 1;
    looseElID.energyInverse[1]		  = 0.14;

    mvaPreElID.cutName                    = "mvaPreElID";
    mvaPreElID.dEtaIn[0]                  = 99999;
    mvaPreElID.dPhiIn[0]                  = 9999;
    mvaPreElID.sigmaIetaIeta[0]           = 0.014;
    mvaPreElID.HadOverEm[0]               = 0.15;
    mvaPreElID.dxy[0]                     = 99999;
    mvaPreElID.dz[0]                      = 99999;
    mvaPreElID.fabsEPDiff[0]              = 99999;
    mvaPreElID.ConversionMissHits[0]      = 99999;
    mvaPreElID.dr03TkSumPt[0]             = 0.2;
    mvaPreElID.dr03EcalRecHitSumEt[0]     = 0.2;
    mvaPreElID.dr03HcalTowerSumEt[0]      = 0.2;
    mvaPreElID.numberOfLostHits[0]        = 0;

    mvaPreElID.dEtaIn[1]                  = 99999;
    mvaPreElID.dPhiIn[1]                  = 99999;
    mvaPreElID.sigmaIetaIeta[1]           = 0.035;
    mvaPreElID.HadOverEm[1]               = 0.10;
    mvaPreElID.dxy[1]                     = 99999;
    mvaPreElID.dz[1]                      = 99999;
    mvaPreElID.fabsEPDiff[1]              = 99999;
    mvaPreElID.ConversionMissHits[1]      = 99999;
    mvaPreElID.dr03TkSumPt[1]             = 0.2;
    mvaPreElID.dr03EcalRecHitSumEt[1]     = 0.2;
    mvaPreElID.dr03HcalTowerSumEt[1]      = 0.2;
    mvaPreElID.numberOfLostHits[1]        = 0;
    
    looseFallMVAElID.cutName              = "looseFall17V2";
    looseMVAElID.cutName                  = "looseMVAElID";
    tightMVAElID.cutName                  = "tightMVAElID";

    looseElIso.cutName                    = "looseElIso";
    looseElIso.chIso04                    = 99999;
    looseElIso.nhIso04                    = 99999;
    looseElIso.phIso04                    = 99999;
    looseElIso.relCombIso04               = 0.4;
    looseElIso.chIso                    = 99999;
    looseElIso.nhIso                    = 99999;
    looseElIso.phIso                    = 99999;
    looseElIso.relCombIso               = 0.15;

    mediumElIso.cutName                   = "mediumElIso";
    mediumElIso.chIso04                   = 99999;
    mediumElIso.nhIso04                   = 99999;
    mediumElIso.phIso04                   = 99999;
    mediumElIso.relCombIso04              = 0.15;
    mediumElIso.chIso                   = 99999;
    mediumElIso.nhIso                   = 99999;
    mediumElIso.phIso                   = 99999;
    mediumElIso.relCombIso              = 0.15;

    tightElIso.cutName                   = "tightElIso";
    tightElIso.chIso04                   = 99999;
    tightElIso.nhIso04                   = 99999;
    tightElIso.phIso04                   = 99999;
    tightElIso.relCombIso04              = 0.1;
    tightElIso.chIso                   = 99999;
    tightElIso.nhIso                   = 99999;
    tightElIso.phIso                   = 99999;
    tightElIso.relCombIso              = 0.1;

    /* Muon ID */
	vetoMuID.cutName 			= "vetoMuID";
	//vetoMuID.pt				= 10;
	//vetoMuID.eta				= 2.4;
	vetoMuID.IsPF				= 1;
	vetoMuID.IsGLB				= 1;
	vetoMuID.NormalizedChi2			= 10;
	vetoMuID.NumberOfValidMuonHits		= 0;
	vetoMuID.NumberOfValidPixelHits		= 0;
	vetoMuID.NumberOfMatchedStations	= 1;
	vetoMuID.TrackLayersWithMeasurement 	= 5; 
	vetoMuID.dz				= 0.2;
	vetoMuID.dxy				= 0.5;
/*
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15*/        
 
    muHZZID.cutName                  = "HZZ";
    TightMuID.cutName                  = "tight";

    looseMuID.cutName			= "looseMuID";
    looseMuID.IsPF			= 1;
    looseMuID.IsGLB			= 1;
    looseMuID.IsTRK			= 1;
      
    tightMuID.cutName                     = "tightMuID";
    tightMuID.IsPF                        = 1;
    tightMuID.IsGLB                       = 1;
    tightMuID.IsTRK                       = 1;
    tightMuID.NormalizedChi2              = 10;
    tightMuID.NumberOfValidMuonHits       = 0;
    tightMuID.NumberOfMatchedStations     = 1;
    tightMuID.NumberOfValidPixelHits      = 0;
    tightMuID.TrackLayersWithMeasurement  = 5;
    tightMuID.dxy                         = 0.2;
    tightMuID.dz                          = 0.5;


    /* Muon ISO */
    amumuMuDetIso.cutName                  = "amumuMuDetIso";
    amumuMuDetIso.hcalIso                  = 99999;
    amumuMuDetIso.ecalIso                  = 99999;
    amumuMuDetIso.trkIso                   = 0.1;
    amumuMuDetIso.relCombIso               = 0.2;

    looseMuDetIso.cutName                  = "looseMuDetIso";
    looseMuDetIso.hcalIso                  = 99999;
    looseMuDetIso.ecalIso                  = 99999;
    looseMuDetIso.trkIso                   = 99999;
    looseMuDetIso.relCombIso               = 0.2;

    tightMuDetIso.cutName                  = "tightMuDetIso";
    tightMuDetIso.hcalIso                  = 99999;
    tightMuDetIso.ecalIso                  = 99999;
    tightMuDetIso.trkIso                   = 99999;
    tightMuDetIso.relCombIso               = 0.2;

    looseMuIso.cutName                    = "looseMuIso";
    looseMuIso.chIso04                    = 99999;
    looseMuIso.nhIso04                    = 99999;
    looseMuIso.phIso04                    = 99999;
    looseMuIso.relCombIso04               = 0.2;

    tightMuIso.cutName                    = "tightMuIso";
    tightMuIso.chIso04                    = 99999;
    tightMuIso.nhIso04                    = 99999;
    tightMuIso.phIso04                    = 99999;
    tightMuIso.relCombIso04               = 0.12;
    
    /* Photon */
    loosePhID.cutName                     = "loosePhID";
    loosePhID.PassedEleSafeVeto[0]        = 1;
    loosePhID.HadOverEm[0]                = 0.0597;
    loosePhID.sigmaIetaIeta[0]            = 0.01031;

    loosePhID.PassedEleSafeVeto[1]        = 1;
    loosePhID.HadOverEm[1]                = 0.0481;
    loosePhID.sigmaIetaIeta[1]            = 0.03013;

    loosePhIso.cutName                    = "loosePhIso";
    loosePhIso.chIso[0]                 = 1.295;
    loosePhIso.nhIso[0]                 = 10.910;
    loosePhIso.phIso[0]                 = 3.630;

    loosePhIso.chIso[1]                 = 1.011;
    loosePhIso.nhIso[1]                 = 5.931;
    loosePhIso.phIso[1]                 = 6.641;

    /*
    mediumPhID.cutName                    = "mediumPhID"; 
    mediumPhID.PassedEleSafeVeto[0]       = 1;
    mediumPhID.HadOverEm[0]               = 0.08;
    mediumPhID.sigmaIetaIeta[0]           = 0.015;

    mediumPhID.PassedEleSafeVeto[1]       = 1;
    mediumPhID.HadOverEm[1]               = 0.05;
    mediumPhID.sigmaIetaIeta[1]           = 0.045;
    */
    /////////////// OLD MEDIUM ID  2017
    mediumPhID.cutName                    = "mediumPhID"; 
    mediumPhID.PassedEleSafeVeto[0]       = 1;
    mediumPhID.HadOverEm[0]               = 0.02197;
    mediumPhID.sigmaIetaIeta[0]           = 0.01015;
    
    mediumPhID.PassedEleSafeVeto[1]       = 1;
    mediumPhID.HadOverEm[1]               = 0.0326;
    mediumPhID.sigmaIetaIeta[1]           = 0.0272;

    /* /////////////// OLD MEDIUM ID
    mediumPhID.PassedEleSafeVeto[0]       = 1;
    mediumPhID.HadOverEm[0]               = 0.0396;
    mediumPhID.sigmaIetaIeta[0]           = 0.01022;

    mediumPhID.PassedEleSafeVeto[1]       = 1;
    mediumPhID.HadOverEm[1]               = 0.0219;
    mediumPhID.sigmaIetaIeta[1]           = 0.03001;
    */

    preSelPhIso.cutName                  = "preSelPhIso";
    preSelPhIso.chIso[0]                 = 0;
    preSelPhIso.nhIso[0]                 = 0;
    preSelPhIso.phIso[0]                 = 15;

    preSelPhIso.chIso[1]                 = 0;
    preSelPhIso.nhIso[1]                 = 0;
    preSelPhIso.phIso[1]                 = 15;

    ////////////// MEDIYM ISO 2017
    mediumPhIso.cutName                   = "mediumPhIso";
    mediumPhIso.chIso[0]                = 1.141;
    mediumPhIso.nhIso[0]                = 1.189;
    mediumPhIso.phIso[0]                = 2.08;

    mediumPhIso.chIso[1]                = 1.015;
    mediumPhIso.nhIso[1]                = 2.718;
    mediumPhIso.phIso[1]                = 3.867;

    /*
    mediumPhIso.cutName                   = "mediumPhIso";
    mediumPhIso.chIso[0]                = 0.441;
    mediumPhIso.nhIso[0]                = 2.725;
    mediumPhIso.phIso[0]                = 2.571;

    mediumPhIso.chIso[1]                = 0.442;
    mediumPhIso.nhIso[1]                = 1.715;
    mediumPhIso.phIso[1]                = 3.863;
    */

    preSelPhIDV2.cutName                    = "preSelPhIDV2";
    preSelPhIDV2.PassedEleSafeVeto[0]       = 1;
    preSelPhIDV2.HadOverEm[0]               = 0.082;
    preSelPhIDV2.sigmaIetaIeta[0]           = 0.014;
    preSelPhIDV2.HcalIso[0]                 = 50;
    preSelPhIDV2.TrkIso[0]                  = 50;
    preSelPhIDV2.ChPfIso[0]                 = 4;

    preSelPhIDV2.PassedEleSafeVeto[1]       = 1;
    preSelPhIDV2.HadOverEm[1]               = 0.075;
    preSelPhIDV2.sigmaIetaIeta[1]           = 0.034;
    preSelPhIDV2.HcalIso[1]                 = 4;
    preSelPhIDV2.TrkIso[1]                  = 4;
    preSelPhIDV2.ChPfIso[1]                 = 4;

    preSelPhID.cutName                    = "preSelPhID";
    preSelPhID.PassedEleSafeVeto[0]       = 1;
    preSelPhID.HadOverEm[0]               = 0.08;
    preSelPhID.sigmaIetaIeta[0]           = 0.015;
    preSelPhID.HcalIso[0]                 = 15;
    preSelPhID.TrkIso[0]                  = 15;
    preSelPhID.ChPfIso[0]                 = 15;

    preSelPhID.PassedEleSafeVeto[1]       = 1;
    preSelPhID.HadOverEm[1]               = 0.05;
    preSelPhID.sigmaIetaIeta[1]           = 0.045;
    preSelPhID.HcalIso[1]                 = 15;
    preSelPhID.TrkIso[1]                  = 15;
    preSelPhID.ChPfIso[1]                 = 15;

    catPhMVAID.cutName                    = "catPhMVAID";
    catPhMVAID.mvaValCat1                 = 0.126;
    catPhMVAID.mvaValCat2                 = 0.107;
    catPhMVAID.mvaValCat3                 = 0.126;
    catPhMVAID.mvaValCat4                 = 0.135;

    //looseMVAPhID.cutName                  = "looseMVAPhID";
    looseMVAPhID.cutName                  = "loose";
    tightMVAPhID.cutName                  = "tightMVAPhID";

    hzgMVAID.cutName                        = "hzgMVAID";
    hzgMVAID.mvaVal[0]                      = -0.9;
    hzgMVAID.mvaVal[1]                      = -0.5;
    hzgMVAID.mvaVal[2]                      = -99;
    hzgMVAID.mvaVal[3]                      = -99;
    hzgMVAID.mvaVal[4]                      = -99;
    hzgMVAID.mvaVal[5]                      = -99;
    hzgMVAID.pt[0]                          =  10;
    hzgMVAID.pt[1]                          =  20;
    hzgMVAID.eta[0]                         = -99;
    hzgMVAID.eta[1]                         = -99;
    hzgMVAID.eta[2]                         = -99;
    hzgMVAID.missHits[0]                    =  0;
    hzgMVAID.missHits[1]                    =  0;
    hzgMVAID.missHits[2]                    =  0;
    hzgMVAID.missHits[3]                    =  0;
    hzgMVAID.missHits[4]                    =  0;
    hzgMVAID.missHits[5]                    =  0;
    hzgMVAID.conversionVeto[0]              =  0;
    hzgMVAID.conversionVeto[1]              =  0;
    hzgMVAID.conversionVeto[2]              =  0;
    hzgMVAID.conversionVeto[3]              =  0;
    hzgMVAID.conversionVeto[4]              =  0;
    hzgMVAID.conversionVeto[5]              =  0;
    hzgMVAID.sip[0]                         =  0;
    hzgMVAID.sip[1]                         =  0;
    hzgMVAID.sip[2]                         =  0;
    hzgMVAID.sip[3]                         =  0;
    hzgMVAID.sip[4]                         =  0;
    hzgMVAID.sip[5]                         =  0;

    hzzMVAID.cutName                        = "hzzMVAID";
    hzzMVAID.mvaVal[0]                      =  -0.211;
    hzzMVAID.mvaVal[1]                      =  -0.396;
    hzzMVAID.mvaVal[2]                      =  -0.215;
    hzzMVAID.mvaVal[3]                      = -0.870;
    hzzMVAID.mvaVal[4]                      = -0.838;
    hzzMVAID.mvaVal[5]                      =  -0.763;
    hzzMVAID.pt[0]                          =  5.0;
    hzzMVAID.pt[1]                          =  10.0;
    hzzMVAID.eta[0]                         =  0.8;
    hzzMVAID.eta[1]                         =  1.479;
    hzzMVAID.eta[2]                         =  2.5;

    hwwMVAID.cutName                        = "hwwMVAID";
    hwwMVAID.mvaVal[0]                      =  0.00;
    hwwMVAID.mvaVal[1]                      =  0.10;
    hwwMVAID.mvaVal[2]                      =  0.62;
    hwwMVAID.mvaVal[3]                      =  0.94;
    hwwMVAID.mvaVal[4]                      =  0.85;
    hwwMVAID.mvaVal[5]                      =  0.92;
    hwwMVAID.pt[0]                          =  10;
    hwwMVAID.pt[1]                          =  20;
    hwwMVAID.eta[0]                         =  0.8;
    hwwMVAID.eta[1]                         =  1.479;
    hwwMVAID.eta[2]                         =  2.5;
    hwwMVAID.missHits[0]                    =  0;
    hwwMVAID.missHits[1]                    =  0;
    hwwMVAID.missHits[2]                    =  0;
    hwwMVAID.missHits[3]                    =  0;
    hwwMVAID.missHits[4]                    =  0;
    hwwMVAID.missHits[5]                    =  0;
    hwwMVAID.conversionVeto[0]              =  1;
    hwwMVAID.conversionVeto[1]              =  1;
    hwwMVAID.conversionVeto[2]              =  1;
    hwwMVAID.conversionVeto[3]              =  1;
    hwwMVAID.conversionVeto[4]              =  1;
    hwwMVAID.conversionVeto[5]              =  1;
    hwwMVAID.sip[0]                         = -99;
    hwwMVAID.sip[1]                         = -99;
    hwwMVAID.sip[2]                         = -99;
    hwwMVAID.sip[3]                         = -99;
    hwwMVAID.sip[4]                         = -99;
    hwwMVAID.sip[5]                         = -99;
    
    //Jet ID
    vbfJetID.cutName                        = "vbfJetID";
    vbfJetID.betaStarC[0]                   = 0.2;
    vbfJetID.dR2Mean[0]                     = 0.06;
    vbfJetID.betaStarC[1]                   = 0.3;
    vbfJetID.dR2Mean[1]                     = 0.05;
    vbfJetID.dR2Mean[2]                     = 0.05;
    vbfJetID.dR2Mean[3]                     = 0.055;

    // Loose Jet ID
    looseJetID.cutName[0]                      = "looseJetID";
    looseJetID.NHF[0]                          = 0.99;
    looseJetID.NEMF[0]                         = 0.99;
    looseJetID.NumConst[0]                     = 1;

    looseJetID.MUF[1]                          = 9999;
    looseJetID.CHF[1]                          = 0;
    looseJetID.CHM[1]                          = 0;
    looseJetID.CEMF[1]                         = 0.99;
    looseJetID.CSV[1]                          = -99;

    looseJetID.cutName[2]                      = "looseJetID";
    looseJetID.NHF[2]                          = 0.98;
    looseJetID.NEMF[2]                         = 0.01;
    looseJetID.NumNeuPart[2]                   = 2;

    looseJetID.cutName[3]                      = "looseJetID";
    looseJetID.NEMF[3]                         = 0.90;
    looseJetID.NumNeuPart[3]                     = 10;

    // tight Jet ID    
    tightJetID.cutName[0]                      = "tightJetID";
    tightJetID.NHF[0]                          = 0.90;
    tightJetID.NEMF[0]                         = 0.90;
    tightJetID.NumConst[0]                     = 1;

    tightJetID.MUF[1]                          = 9999;
    tightJetID.CHF[1]                          = 0;
    tightJetID.CHM[1]                          = 0;
    tightJetID.CEMF[1]                         = 0.99;
    tightJetID.CSV[1]                          = -99;

    tightJetID.cutName[2]                      = "tightJetID";
    tightJetID.NHF[2]                          = 0.98;
    tightJetID.NEMF[2]                         = 0.01;
    tightJetID.NumNeuPart[2]                   = 2;

    tightJetID.cutName[3]                      = "tightJetID";
    tightJetID.NEMF[3]                         = 0.90;
    tightJetID.NumNeuPart[3]                   = 10;
    
/*
    bJetID.cutName                          = "bJetID";
    bJetID.NHF                              = 0.99;
    bJetID.NEMF                             = 0.99;
    bJetID.NumConst                         = 1;
    bJetID.MUF                              = 0.8;
    bJetID.CHF                              = 0;
    bJetID.CHM                              = 0;
    bJetID.CEMF                             = 0.99;
    bJetID.CSV                              = 0.898;  // medium WP
*/
}
