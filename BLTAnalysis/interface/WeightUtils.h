/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

// c++ libraries
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        float   GetPUWeight(float);
        pair<float, float>   GetTriggerEffWeight(string, TLorentzVector&) const;
        float   GetMuonIDEff(TLorentzVector&) const; 
        float   GetMuonISOEff(TLorentzVector&) const; 
        float   GetElectronRecoIdEff(TLorentzVector&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        TGraph  *_puReweight;
        TGraphAsymmErrors *_eff_IsoMu24_DATA[4]; 
        TGraphAsymmErrors *_muSF_ID_DATA[4], *_muSF_ID_MC[4]; 
        TGraphAsymmErrors *_muSF_ISO_DATA[4], *_muSF_ISO_MC[4]; 

        TGraphErrors *_eleSF_RECO, *_eleSF_ID[5];

        //TGraph *_muSF_ID_err[4], *_muSF_ISO_err[4];

        //TH2D    *h2_MuTriggerSFs[2]; // Good for Mu17_Mu8 or Mu17_TkMu8
        //TH2D    *h2_EleMVASF;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
