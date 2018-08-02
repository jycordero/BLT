#include "BLT/BLTAnalysis/interface/WeightUtils.h"


// helper functions

template <class Graph>
int GetBinNumber(Graph graph, float x0)
{
    Double_t x,y;
    float xmin = 1e9;
    for (int i = 0; i < graph->GetN(); i++) {
        graph->GetPoint(i, x, y);
        float diff = fabs(x - x0);
        if (diff > xmin) {
            return i - 1;
        } else {
            xmin = diff;
        }
    }
    return graph->GetN() - 1; 
}



// definitions for WeightUtils


WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;
    std::string fileName;

    const std::string cmssw_base = getenv("CMSSW_BASE");

    // PU weights
    std::string puFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup_sf_2016_full.root";
    TFile* puFile = new TFile(puFileName.c_str(), "OPEN");
    _puReweight = (TGraph*)puFile->Get("pileup_sf");

    // muon trigger sf (BCDEF)
//  std::string triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_BCDEF.root";
    std::string triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_RunBtoF.root";
    TFile* triggerFile_BCDEF = new TFile(triggerFileName.c_str(), "OPEN");

    std::string filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _muSF_IsoMu24_DATA_BCDEF[0] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[1] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[2] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[3] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/";
    _muSF_IsoMu24_MC_BCDEF[0] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[1] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[2] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[3] = (TGraphAsymmErrors*)triggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());

    // muon trigger sf (GH)
//  triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_GH.root";
    triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_Period4.root";
    TFile* triggerFile_GH = new TFile(triggerFileName.c_str(), "OPEN");

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _muSF_IsoMu24_DATA_GH[0] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[1] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[2] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[3] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/";
    _muSF_IsoMu24_MC_GH[0] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[1] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[2] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[3] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _muSF_IsoMu24_DATA_GH[0] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[1] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[2] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[3] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/";
    _muSF_IsoMu24_MC_GH[0] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[1] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[2] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[3] = (TGraphAsymmErrors*)triggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());


    // muon tight ID sf (BCDEF)
    std::string idFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF2012_ID = new TFile(idFileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // muon tight ID sf (GH)
    idFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_ID_GH = new TFile(idFileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // muon loose ID sf (BCDEF)
    idFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF_Loose_ID_BCDEF = new TFile(idFileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // muon loose ID sf (GH)
    idFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_Loose_ID_GH = new TFile(idFileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // loose muon ISO sf (BCDEF)
    std::string isoFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF2012_ISO = new TFile(isoFileName.c_str(), "OPEN"); 

    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesDATA/";
    _muSF2012_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_DATA").c_str());

    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesMC/";
    _muSF2012_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_MC").c_str());
  
    // hzz muon id efficiencies
    // (from topic_wbranch)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/hzz_muon_id_sf.root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2F*)f_hzz_muIdSF->Get("FINAL");
    _hzz_muIdErr = (TH2F*)f_hzz_muIdSF->Get("ERROR");
  
    // tight muon ISO sf (BCDEF)
    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ISO->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

    // tight muon ISO sf (GH)
    isoFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_ISO_GH = new TFile(isoFileName.c_str(), "OPEN"); 

    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

    // loose muon ISO sf (BCDEF)
    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesDATA/";
    _muSF2012_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_DATA").c_str());
    _muSF2012_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_DATA").c_str());

    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesMC/";
    _muSF2012_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_MC").c_str());
    _muSF2012_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_MC").c_str());

    // electron reco efficiencies
    // (from topic_wbranch)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_reco_2016.root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");

    // electron id efficiencies
    // (from topic_wbranch)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");

    // hzz electron id efficiencies
    // (from topic_jbueghly)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egamma_eff_hzz_ID_2016.root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleSF_ID[0] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_0");
    _hzz_eleSF_ID[1] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_1");
    _hzz_eleSF_ID[2] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_2");
    _hzz_eleSF_ID[3] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_3");
    _hzz_eleSF_ID[4] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_4");
    _hzz_eleSF_ID[5] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_5");
    _hzz_eleSF_ID[6] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_6");
    _hzz_eleSF_ID[7] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_7");
    _hzz_eleSF_ID[8] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_8");
    _hzz_eleSF_ID[9] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_9");
    _hzz_eleSF_ID[10] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_10");
    _hzz_eleSF_ID[11] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_11");
    _hzz_eleSF_ID[12] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_12");
}



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

EfficiencyContainer WeightUtils::GetPUEff(float nPU)
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    effData = _puReweight->Eval(nPU);
    errData = sqrt(0.01 * abs(effData));

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}




EfficiencyContainer WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    if (triggerName == "HLT_IsoMu24_v*")
    {
        float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
        int etaBin = 0;
        for (int i = 0; i < 4; i++)
        {
            if (fabs(lepton.Eta()) > binningEta[i] && fabs(lepton.Eta()) <= binningEta[i+1])
            {
                etaBin = i;     break;
            }
        }
        if (_dataPeriod == "2016BtoF")
        {
            effData = _muSF_IsoMu24_DATA_BCDEF[etaBin]->Eval(lepton.Pt());
            effMC   = _muSF_IsoMu24_MC_BCDEF[etaBin]->Eval(lepton.Pt());

            int ptBin = GetBinNumber(_muSF_IsoMu24_DATA_BCDEF[etaBin], lepton.Pt()); 
            errData = _muSF_IsoMu24_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_IsoMu24_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        }
        else if (_dataPeriod == "2016GH")
        {
            effData = _muSF_IsoMu24_DATA_GH[etaBin]->Eval(lepton.Pt());
            effMC   = _muSF_IsoMu24_MC_GH[etaBin]->Eval(lepton.Pt());

            int ptBin = GetBinNumber(_muSF_IsoMu24_DATA_GH[etaBin], lepton.Pt()); 
            errData = _muSF_IsoMu24_DATA_GH[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_IsoMu24_MC_GH[etaBin]->GetErrorY(ptBin);
        }
    }

    else if (triggerName == "HLT_Ele27_WPTight_Gsf_v*")
    {
        int etaBin = 0;
        int ptBin = 0;
        for (int i = 0; i < 13; i++)
        {
            if (lepton.Eta() > _eleEtaBins[i] && lepton.Eta() <= _eleEtaBins[i+1])
            {
                etaBin = i;     break;
            }
        }
        for (int i = 0; i < 8; i++)
        {
            if (lepton.Pt() > _elePtBins[i] && lepton.Pt() <= _elePtBins[i+1])
            {
                ptBin = i;      break;
            }
        }
        effData = _ele_trigEff_data[etaBin][ptBin];
        effMC   = _ele_trigEff_mc[etaBin][ptBin];
        errData = 0.005*_ele_trigEff_data[etaBin][ptBin];
        errMC   = 0.005*_ele_trigEff_mc[etaBin][ptBin];
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}


EfficiencyContainer WeightUtils::GetMuonIDEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; i++)
    {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1])
        {
            etaBin = i;     break;
        }
    }

    if (_dataPeriod == "2016BtoF")
    {
        effData = _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_BCDEF[etaBin], muon.Pt());
        errData = _muSF_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }

    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_GH[etaBin], muon.Pt());
        errData = _muSF_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_ID_MC_GH[etaBin]->GetErrorY(ptBin);
    }
    
    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



// (from topic_jbueghly)
EfficiencyContainer WeightUtils::GetLooseMuonIDEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; i++)
    {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1])
        {
            etaBin = i;     break;
        }
    }

    if (muon.Pt() < 200.)
    {
        if (_dataPeriod == "2016BtoF")
        {
            effData = _muSF_Loose_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
            effMC   = _muSF_Loose_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());

            int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Loose_ID_DATA_BCDEF[etaBin], muon.Pt());
            errData = _muSF_Loose_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_Loose_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        }
        else if (_dataPeriod == "2016GH")
        {
            effData = _muSF_Loose_ID_DATA_GH[etaBin]->Eval(muon.Pt());
            effMC   = _muSF_Loose_ID_MC_GH[etaBin]->Eval(muon.Pt());

            int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Loose_ID_DATA_GH[etaBin], muon.Pt());
            errData = _muSF_Loose_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_Loose_ID_MC_GH[etaBin]->GetErrorY(ptBin);
        }
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}


// "Corrected"
EfficiencyContainer WeightUtils::GetHZZMuonIDEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningEta[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};
    float binningPt[] = {5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60, 80, 9999};
    unsigned etaBin = 0, ptBin = 0;
    for (unsigned i = 0; i < 15; i++)
    {
        if (muon.Eta() >= binningEta[i] && muon.Eta() < binningEta[i+1])
        {
            etaBin = i + 1;
            break;
        }
    }
    for (unsigned i = 0; i < 15; i++)
    {
        if (muon.Pt() >= binningPt[i] && muon.Pt() < binningPt[i+1])
        {
            ptBin = i + 1;
            break;
        }
    }
    effData = _hzz_muIdSF->GetBinContent(etaBin, ptBin);
    errData = _hzz_muIdErr->GetBinContent(etaBin, ptBin);

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



EfficiencyContainer WeightUtils::GetMuonLooseISOEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; i++)
    {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1])
        {
            etaBin = i;     break;
        }
    }

    if (_dataPeriod == "2016BtoF")
    {
        effData = _muSF2012_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF2012_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF2012_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        errData = _muSF2012_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF2012_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }
    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF2012_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF2012_ISO_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF2012_ISO_DATA_GH[etaBin], muon.Pt());
        errData = _muSF2012_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF2012_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetMuonTightISOEff(TLorentzVector& muon) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; i++)
    {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1])
        {
            etaBin = i;     break;
        }
    }

    if (_dataPeriod == "2016BtoF")
    {
        effData = _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        errData = _muSF_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }
    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_GH[etaBin], muon.Pt());
        errData = _muSF_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronRecoEff(TLorentzVector& electron) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningPt[] = {10., 20., 35., 50., 90., 500.};
    int ptBin = 0;
    for (int i = 0; i < 5; i++)
    {
        if (fabs(electron.Pt()) > binningPt[i] && fabs(electron.Pt()) <= binningPt[i+1])
        {
            ptBin = i;      break;
        }
    }

    float sfReco = _eleSF_RECO->Eval(electron.Eta());
    float sfId   = _eleSF_ID[ptBin]->Eval(electron.Eta());
    effData = sfReco * sfId;

    int etaBin;
    etaBin = GetBinNumber<TGraphErrors*>(_eleSF_RECO, electron.Eta()); 
    float errReco = _eleSF_RECO->GetErrorY(etaBin);

    etaBin = GetBinNumber<TGraphErrors*>(_eleSF_ID[ptBin], electron.Eta()); 
    float errId  = _eleSF_ID[ptBin]->GetErrorY(etaBin);
    errData = sfReco * sfId * (pow(errReco/sfReco, 2) + pow(errId/sfId, 2));

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetHZZElectronRecoEff(TElectron& electron) const 
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    float binningPt[] = {7., 15., 20., 30., 40., 50., 60., 70., 80., 100., 120., 140., 160., 9999}; 
    int ptBin = 0;
    for (int i = 0; i < 13; i++)
    {
        if (fabs(electron.calibPt) > binningPt[i] && fabs(electron.calibPt) <= binningPt[i+1])
        {
            ptBin = i;      break;
        }
    }

    float sfReco = _eleSF_RECO->Eval(electron.scEta);
    float sfId   = _hzz_eleSF_ID[ptBin]->Eval(electron.scEta);
    effData = sfReco * sfId;

    int etaBin;
    etaBin = GetBinNumber<TGraphErrors*>(_eleSF_RECO, electron.scEta); 
    float errReco = _eleSF_RECO->GetErrorY(etaBin);

    etaBin = GetBinNumber<TGraphErrors*>(_hzz_eleSF_ID[ptBin], electron.scEta); 
    float errId  = _hzz_eleSF_ID[ptBin]->GetErrorY(etaBin);
    errData = sfReco*sfId*(pow(errReco/sfReco, 2) + pow(errId/sfId, 2));

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



//
// Definitions for EfficiencyContainer
//

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
