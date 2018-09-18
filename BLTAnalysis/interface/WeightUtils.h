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

// BaconAna class definitions
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"

using namespace std;


class EfficiencyContainer: public TObject{
    public:
        EfficiencyContainer();
        EfficiencyContainer(float, float, float, float);
        virtual ~EfficiencyContainer() {};
        void SetData(float, float, float, float);
        pair<double, double> GetEff() {return make_pair(_dataEff, _mcEff);};
        pair<double, double> GetErr() {return make_pair(_dataErr, _mcErr);};
        float GetSF() {return (_dataEff/_mcEff);};
        float GetVar() {return pow(_dataEff/_mcEff, 2)*(pow(_dataErr/_dataEff, 2) + pow(_mcErr/_mcEff, 2));};

    private:
        float _dataEff, _mcEff;
        float _dataErr, _mcErr;
};


class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        float               GetPUWeight(float);
        EfficiencyContainer GetSingleMuonTriggerEff(TLorentzVector&) const;
        EfficiencyContainer GetSingleElectronTriggerEff(TLorentzVector&) const;
        EfficiencyContainer GetDoubleMuonTriggerEff(TLorentzVector&, int) const;
        EfficiencyContainer GetDoubleElectronTriggerEff(const baconhep::TElectron*, int) const;
        EfficiencyContainer GetTriggerEff(string, TLorentzVector&) const;
        EfficiencyContainer GetHZZMuonIDEff(TLorentzVector&) const;
        EfficiencyContainer GetHZZElectronIDRecoEff(const baconhep::TElectron*) const;

        ClassDef(WeightUtils, 0);

    private:
        // Input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // Pileup
        TH1D *_puReweight;

        // Triggers
        TGraphErrors *_sf_IsoMu24_Eta2p1_data[3], *_sf_IsoMu24_Eta2p1_mc[3];


        // Lepton ID
        TH2D *_hzz_muIdSF;
        TH2F *_hzz_eleIdSF;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
