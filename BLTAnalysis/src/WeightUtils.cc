#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;
    std::string fileName;

    const std::string cmssw_base = getenv("CMSSW_BASE");




    //////////////////
    //    PILEUP    //
    //////////////////

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/PUWeights_2012.root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TH1D*)puFile->Get("pileup");




    //////////////////
    //   TRIGGERS   //
    //////////////////

    //--- SINGLE LEPTON ---//

    // Muon
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root";
    TFile* triggerFile = new TFile(fileName.c_str(), "OPEN");
    _sf_IsoMu24_Eta2p1_data[0] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_data[1] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_data[2] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");

    _sf_IsoMu24_Eta2p1_mc[0] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_mc[1] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_mc[2] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");


    //--- DOUBLE LEPTON ---//

    // Muon
//  fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root";



    //////////////////
    //      ID      //
    //////////////////

    //--- HZZ MUON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/MuonScaleFactors_2011_2012.root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2D*) f_hzz_muIdSF->Get("TH2D_ALL_2012");



    //--- HZZ ELECTRON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/CombinedMethod_ScaleFactors_RecoIdIsoSip.root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleIdSF = (TH2F*) f_hzz_eleIdSF->Get("h_electronScaleFactor_RecoIdIsoSip");
}





//--- FUNCTIONS ---//



//////////////////
//  UTILITIES   //
//////////////////


void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}



//////////////////
//    PILEUP    //
//////////////////

float WeightUtils::GetPUWeight(float nPU)
{
    return _puReweight->GetBinContent(_puReweight->FindBin(nPU)); 
}


//////////////////
//   TRIGGERS   //
//////////////////


//--- DOUBLE LEPTON ---//


// Muons
// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs#2012_data
// FIXME apparently these histograms are not suitable for 4l final state
// https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2012#Reference_Efficiencies_for_2012
EfficiencyContainer WeightUtils::GetDoubleMuonTriggerEff(TLorentzVector &muon, int leg) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}


// Electrons
// FIXME
EfficiencyContainer WeightUtils::GetDoubleElectronTriggerEff(const baconhep::TElectron* electron, int leg) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



//--- SINGLE LEPTON ---//

// Muon
EfficiencyContainer WeightUtils::GetSingleMuonTriggerEff(TLorentzVector &muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    float binningEta[] = {0., 0.9, 1.2, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 3; i++)
    {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1])
        {
            etaBin = i;
            break;
        }
    }

    if (muon.Pt() < 300.)
    {
        effMC   = _sf_IsoMu24_Eta2p1_data[etaBin]->Eval(muon.Pt());
        effData = _sf_IsoMu24_Eta2p1_mc[etaBin]->Eval(muon.Pt());
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}


// Electron
// FIXME
EfficiencyContainer WeightUtils::GetSingleElectronTriggerEff(TLorentzVector &muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}




//////////////////
//      ID      //
//////////////////


//--- HZZ MUON ---//

// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Lepton_efficiency_scale_factors
EfficiencyContainer WeightUtils::GetHZZMuonIDEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float minPt = 5, maxPt = 100, maxEta = 2.4;
    if (fabs(muon.Eta() < maxEta))
    {
        int bin;

        if (muon.Pt() > maxPt)
            bin = _hzz_muIdSF->FindBin(0.99 * maxPt, muon.Eta());
        else if (muon.Pt() < minPt)
            bin = _hzz_muIdSF->FindBin(1.01 * minPt, muon.Eta());
        else
            bin = _hzz_muIdSF->FindBin(muon.Pt(), muon.Eta());

        effData = _hzz_muIdSF->GetBinContent(bin);
        errData = _hzz_muIdSF->GetBinError(bin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



//--- HZZ ELECTRON ---//

// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Lepton_efficiency_scale_factors
EfficiencyContainer WeightUtils::GetHZZElectronIDRecoEff(const baconhep::TElectron* electron) const 
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float minPt = 7, maxPt = 200, maxEta = 2.5;
    int bin;
    if (fabs(electron->scEta) < maxEta)     // I am actually not sure if this is really supposed to be scEta :')
    {   
        if (electron->ptHZZ4l > maxPt)
            bin = _hzz_eleIdSF->FindBin(0.99 * maxPt, electron->scEta);
        else if (electron->ptHZZ4l < minPt)
            bin = _hzz_eleIdSF->FindBin(1.01 * minPt, electron->scEta);
        else
            bin = _hzz_eleIdSF->FindBin(electron->ptHZZ4l, electron->scEta);

        effData = _hzz_eleIdSF->GetBinContent(bin);
        errData = _hzz_eleIdSF->GetBinError(bin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}




//////////////////
//   EFF CONT   //
//////////////////

EfficiencyContainer::EfficiencyContainer()
{
    _dataEff = 0.;
    _mcEff   = 0.;
    _dataErr = 0.;
    _mcErr   = 0.;
};

EfficiencyContainer::EfficiencyContainer(float dataEff, float mcEff, float dataErr, float mcErr)
{
    _dataEff = dataEff;
    _mcEff   = mcEff;
    _dataErr = dataErr;
    _mcErr   = mcErr;
};

void EfficiencyContainer::SetData(float dataEff, float mcEff, float dataErr, float mcErr)
{
    _dataEff = dataEff;
    _mcEff   = mcEff;
    _dataErr = dataErr;
    _mcErr   = mcErr;
};
