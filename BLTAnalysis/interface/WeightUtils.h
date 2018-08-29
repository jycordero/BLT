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
        EfficiencyContainer GetDoubleMuonTriggerEff(string, int, TLorentzVector&) const;
        EfficiencyContainer GetDoubleElectronTriggerEff(string, int, const baconhep::TElectron*) const;
        EfficiencyContainer GetTriggerEff(string, TLorentzVector&) const;
        EfficiencyContainer GetMuonIDEff(TLorentzVector&) const; 
        EfficiencyContainer GetLooseMuonIDEff(TLorentzVector&) const;
        EfficiencyContainer GetMuonTightISOEff(TLorentzVector&) const; 
        EfficiencyContainer GetMuonLooseISOEff(TLorentzVector&) const; 
        EfficiencyContainer GetHZZMuonIDEff(TLorentzVector&) const;
        EfficiencyContainer GetElectronRecoEff(TLorentzVector&) const;
        EfficiencyContainer GetHZZElectronRecoEff(const baconhep::TElectron*) const;

        ClassDef(WeightUtils, 0);

    private:
        // Input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // Pileup
        TGraph  *_puReweight;

        // Muon triggers, ID, iso
        TGraphAsymmErrors *_eff_doubleMu_leg1_DATA[4], *_eff_doubleMu_leg1_MC[4];
        TGraphAsymmErrors *_eff_doubleMu_leg2_DATA[4], *_eff_doubleMu_leg2_MC[4];
        TH2F *_hzz_muIdSF, *_hzz_muIdErr;

        // Electron RECO/ID scale factors (what about ISO?)
        TH2F *_eleSF_RECO, *_eleSF_RECO_lowEt;
        TH2F *_hzz_eleIdSF, *_hzz_eleIdSF;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
