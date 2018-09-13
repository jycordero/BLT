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




    //--- ELECTRON ID ---//
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for

    // Tight
    tightElID.cutName                    = "tightElID";

    tightElID.sigmaIetaIeta[0]           = 0.00998;     // barrel (|eta| <= 1.479)
    tightElID.dEtaIn[0]                  = 0.00308;
    tightElID.dPhiIn[0]                  = 0.0816;
    tightElID.HadOverEm[0]               = 0.0414;
    tightElID.fabsEPDiff[0]              = 0.0129;
    tightElID.ConversionMissHits[0]      = 1;

    tightElID.sigmaIetaIeta[1]           = 0.0292;      // endcap (|eta| > 1.479)
    tightElID.dEtaIn[1]                  = 0.00605;
    tightElID.dPhiIn[1]                  = 0.0394;
    tightElID.HadOverEm[1]               = 0.0641;
    tightElID.fabsEPDiff[1]              = 0.0129;
    tightElID.ConversionMissHits[1]      = 1;


    // Medium
    mediumElID.cutName                    = "mediumElID";

    mediumElID.sigmaIetaIeta[0]           = 0.00998;    // barrel (|eta| <= 1.479)
    mediumElID.dEtaIn[0]                  = 0.00311;
    mediumElID.dPhiIn[0]                  = 0.103;
    mediumElID.HadOverEm[0]               = 0.253;
    mediumElID.fabsEPDiff[0]              = 0.134;
    mediumElID.ConversionMissHits[0]      = 1;

    mediumElID.sigmaIetaIeta[1]           = 0.0298;     // endcap (|eta| > 1.479)
    mediumElID.dEtaIn[1]                  = 0.00609;
    mediumElID.dPhiIn[1]                  = 0.045;
    mediumElID.HadOverEm[1]               = 0.0878;
    mediumElID.fabsEPDiff[1]              = 0.13;
    mediumElID.ConversionMissHits[1]      = 1;


    // Loose
    looseElID.cutName                     = "looseElID";

    looseElID.sigmaIetaIeta[0]            = 0.011;      // barrel (|eta| <= 1.479)
    looseElID.dEtaIn[0]                   = 0.00477;
    looseElID.dPhiIn[0]                   = 0.222;
    looseElID.HadOverEm[0]                = 0.298;
    looseElID.fabsEPDiff[0]               = 0.241;
    looseElID.ConversionMissHits[0]       = 1;

    looseElID.sigmaIetaIeta[1]            = 0.0314;     // endcap (|eta| > 1.479)
    looseElID.dEtaIn[1]                   = 0.00868;
    looseElID.dPhiIn[1]                   = 0.213;
    looseElID.HadOverEm[1]                = 0.101;
    looseElID.fabsEPDiff[1]               = 0.14;
    looseElID.ConversionMissHits[1]       = 1;


    // Veto
    vetoElID.cutName                      = "vetoElID";

    vetoElID.sigmaIetaIeta[0]             = 0.0115;     // barrel (|eta| <= 1.479)
    vetoElID.dEtaIn[0]                    = 0.00749;
    vetoElID.dPhiIn[0]                    = 0.228;
    vetoElID.HadOverEm[0]                 = 0.356;
    vetoElID.fabsEPDiff[0]                = 0.299;
    vetoElID.ConversionMissHits[0]        = 2;

    vetoElID.sigmaIetaIeta[1]             = 0.037;      // endcap (|eta| > 1.479)
    vetoElID.dEtaIn[1]                    = 0.00895;
    vetoElID.dPhiIn[1]                    = 0.213;
    vetoElID.HadOverEm[1]                 = 0.211;
    vetoElID.fabsEPDiff[1]                = 0.15;
    vetoElID.ConversionMissHits[1]        = 3;


    // HZZ electron MVA
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2017#Electrons
    hzzMVA.cutName                        = "hzzMVA";
    hzzMVA.bdt[0]                         = -0.211;
    hzzMVA.bdt[1]                         = -0.396;
    hzzMVA.bdt[2]                         = -0.215;
    hzzMVA.bdt[3]                         = -0.870;
    hzzMVA.bdt[4]                         = -0.838;
    hzzMVA.bdt[5]                         = -0.763;
    hzzMVA.pt[0]                          = 5;
    hzzMVA.pt[1]                          = 10;
    hzzMVA.eta[0]                         = 0.8;
    hzzMVA.eta[1]                         = 1.479;
    hzzMVA.eta[2]                         = 2.5;



    //--- ELECTRON ISO ---//

    // Tight
    tightElIso.cutName                  = "tightElIso";
    tightElIso.relCombIso[0]            = 0.0588;       // barrel (|eta| <= 1.479)
    tightElIso.relCombIso[1]            = 0.0571;       // endcap (|eta| > 1.479)

    // Medium
    mediumElIso.cutName                 = "mediumElIso";
    mediumElIso.relCombIso[0]           = 0.0695;       // barrel (|eta| <= 1.479)
    mediumElIso.relCombIso[1]           = 0.0821;       // endcap (|eta| > 1.479)

    // Loose
    looseElIso.cutName                  = "looseElIso";
    looseElIso.relCombIso[0]            = 0.0994;       // barrel (|eta| <= 1.479)
    looseElIso.relCombIso[1]            = 0.107;        // endcap (|eta| > 1.479)
    
    // Veto
    vetoElIso.cutName                   = "vetoElIso";
    vetoElIso.relCombIso[0]             = 0.175;        // barrel (|eta| <= 1.479)
    vetoElIso.relCombIso[1]             = 0.159;        // endcap (|eta| > 1.479)



    //--- MUON ID ---//

    // Tight
    tightMuID.cutName                     = "tightMuID";
    tightMuID.IsPF                        = 1;
    tightMuID.IsGLB                       = 1;
    tightMuID.IsTRK                       = 1;
    tightMuID.NormalizedChi2              = 10.;
    tightMuID.NumberOfValidMuonHits       = 0;
    tightMuID.NumberOfMatchedStations     = 1;
    tightMuID.NumberOfValidPixelHits      = 0;
    tightMuID.TrackLayersWithMeasurement  = 5;
    tightMuID.dxy                         = 0.2;
    tightMuID.dz                          = 0.5;



    //--- MUON ISO ---//
    amumuMuDetIso.cutName                    = "amumuMuDetIso";
    amumuMuDetIso.hcalIso                  = 99999;
    amumuMuDetIso.ecalIso                  = 99999;
    amumuMuDetIso.trkIso                   = 0.1;
    amumuMuDetIso.relCombIso               = 0.2;

    looseMuDetIso.cutName                    = "looseMuDetIso";
    looseMuDetIso.hcalIso                  = 99999;
    looseMuDetIso.ecalIso                  = 99999;
    looseMuDetIso.trkIso                   = 99999;
    looseMuDetIso.relCombIso               = 0.2;

    tightMuDetIso.cutName                    = "tightMuDetIso";
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



    //--- PHOTONS ----//

    loosePhID.cutName                     = "loosePhID";
    loosePhID.PassedEleSafeVeto[0]        = 1;
    loosePhID.HadOverEm[0]                = 0.05;
    loosePhID.sigmaIetaIeta[0]            = 0.012;

    loosePhID.PassedEleSafeVeto[1]        = 1;
    loosePhID.HadOverEm[1]                = 0.05;
    loosePhID.sigmaIetaIeta[1]            = 0.034;

    loosePhIso.cutName                    = "loosePhIso";
    loosePhIso.chIso[0]                 = 2.6;
    loosePhIso.nhIso[0]                 = 3.5;
    loosePhIso.phIso[0]                 = 1.3;

    loosePhIso.chIso[1]                 = 2.3;
    loosePhIso.nhIso[1]                 = 2.9;
    loosePhIso.phIso[1]                 = 99999;

    mediumPhID.cutName                    = "mediumPhID";
    mediumPhID.PassedEleSafeVeto[0]       = 1;
    mediumPhID.HadOverEm[0]               = 0.05;
    mediumPhID.sigmaIetaIeta[0]           = 0.011;

    mediumPhID.PassedEleSafeVeto[1]       = 1;
    mediumPhID.HadOverEm[1]               = 0.05;
    mediumPhID.sigmaIetaIeta[1]           = 0.033;

    mediumPhIso.cutName                   = "mediumPhIso";
    mediumPhIso.chIso[0]                = 1.5;
    mediumPhIso.nhIso[0]                = 1.0;
    mediumPhIso.phIso[0]                = 0.7;

    mediumPhIso.chIso[1]                = 1.2;
    mediumPhIso.nhIso[1]                = 1.5;
    mediumPhIso.phIso[1]                = 1.0;

    preSelPhID.cutName                    = "preSelPhID";
    preSelPhID.PassedEleSafeVeto[0]       = 1;
    preSelPhID.HadOverEm[0]               = 0.082;
    preSelPhID.sigmaIetaIeta[0]           = 0.014;
    preSelPhID.HcalIso[0]                 = 50;
    preSelPhID.TrkIso[0]                  = 50;
    preSelPhID.ChPfIso[0]                 = 4;

    preSelPhID.PassedEleSafeVeto[1]       = 1;
    preSelPhID.HadOverEm[1]               = 0.075;
    preSelPhID.sigmaIetaIeta[1]           = 0.034;
    preSelPhID.HcalIso[1]                 = 4;
    preSelPhID.TrkIso[1]                  = 4;
    preSelPhID.ChPfIso[1]                 = 4;

    catPhMVAID.cutName                    = "catPhMVAID";
    catPhMVAID.mvaValCat1                 = 0.126;
    catPhMVAID.mvaValCat2                 = 0.107;
    catPhMVAID.mvaValCat3                 = 0.126;
    catPhMVAID.mvaValCat4                 = 0.135;

    vbfJetID.cutName                        = "vbfJetID";
    vbfJetID.betaStarC[0]                   = 0.2;
    vbfJetID.dR2Mean[0]                     = 0.06;
    vbfJetID.betaStarC[1]                   = 0.3;
    vbfJetID.dR2Mean[1]                     = 0.05;
    vbfJetID.dR2Mean[2]                     = 0.05;
    vbfJetID.dR2Mean[3]                     = 0.055;

    looseJetID.cutName                      = "looseJetID";
    looseJetID.NHF                          = 0.99;
    looseJetID.NEMF                         = 0.99;
    looseJetID.NumConst                     = 1;
    looseJetID.MUF                          = 0.8;
    looseJetID.CHF                          = 0;
    looseJetID.CHM                          = 0;
    looseJetID.CEMF                         = 0.99;
    looseJetID.CSV                          = -99;

    bJetID.cutName                          = "bJetID";
    bJetID.NHF                              = 0.99;
    bJetID.NEMF                             = 0.99;
    bJetID.NumConst                         = 1;
    bJetID.MUF                              = 0.8;
    bJetID.CHF                              = 0;
    bJetID.CHM                              = 0;
    bJetID.CEMF                             = 0.99;
    bJetID.CSV                              = 0.898;  // medium WP
}
