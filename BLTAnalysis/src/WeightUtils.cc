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




    //////////////////
    //    PILEUP    //
    //////////////////


    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup_sf_2016_full.root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TGraph*)puFile->Get("pileup_sf");




    //////////////////
    //   TRIGGERS   //
    //////////////////


    //--- SINGLE MUON ---//

    // BCDEF
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_RunBtoF.root";
    TFile* triggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");

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


    // GH
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_Period4.root";
    TFile* triggerFile_GH = new TFile(fileName.c_str(), "OPEN");

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


    
    //--- DOUBLE MUON ---//

    // Leg 1
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu17Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg1_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_data");
    _eff_doubleMu_leg1_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu17Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg1_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_data");
    _eff_doubleMu_leg1_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu17Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg1_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_data");
    _eff_doubleMu_leg1_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu17Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg1_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_data");
    _eff_doubleMu_leg1_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_mc"); 


    // Leg 2
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu8Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg2_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_data");
    _eff_doubleMu_leg2_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu8Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg2_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_data");
    _eff_doubleMu_leg2_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu8Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg2_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_data");
    _eff_doubleMu_leg2_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/sf_Mu8Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg2_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_data");
    _eff_doubleMu_leg2_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_mc"); 



    //--- DOUBLE ELECTRON ---//

    // Leg 1
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/SFs_Leg1_Ele23_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleEle_leg1_DATA = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffData2D");
    _eff_doubleEle_leg1_MC   = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffMC2D"); 

    // Leg 2
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/SFs_Leg2_Ele12_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleEle_leg2_DATA = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffData2D");
    _eff_doubleEle_leg2_MC   = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffMC2D");




    //////////////////
    //      ID      //
    //////////////////

/*
    //--- TIGHT MUON ---//

    // BCDEF
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_BCDEF.root";
    TFile* f_muSF_ID_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Tight_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Tight_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Tight_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Tight_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Tight_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Tight_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Tight_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Tight_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());


    // GH
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_GH.root";
    TFile* f_muSF_ID_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Tight_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Tight_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Tight_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Tight_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Tight_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Tight_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Tight_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Tight_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());



    //--- LOOSE MUON ---//

    // BCDEF
    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

  
    // GH
    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());
*/


    //--- HZZ MUON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/hzz_muon_id_sf.root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2F*)f_hzz_muIdSF->Get("FINAL");
    _hzz_muIdErr = (TH2F*)f_hzz_muIdSF->Get("ERROR");


/*
    //--- (TIGHT?) ELECTRON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");
*/


    //--- HZZ ELECTRON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/hzz_electron_id_sf.root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleIdSF = (TH2F*)f_hzz_eleIdSF->Get("EGamma_SF2D");




    //////////////////
    //     RECO     //
    //////////////////


    //--- ELECTRON ---//

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egammaEffi.txt_EGM2D.root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN");
    _eleSF_RECO = (TH2F*) f_eleRecoSF->Get("EGamma_SF2D");




    //////////////////
    //     ISO      //
    //////////////////


    //--- TIGHT MUON ---//

    // BCDEF
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_BCDEF.root";
    TFile* f_muSF_ISO_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_Tight_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_Tight_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());


    // GH
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_GH.root";
    TFile* f_muSF_ISO_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_Tight_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_Tight_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_Tight_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_Tight_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());



    //--- LOOSE MUON ---//

    // BCDEF
    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesDATA/";
    _muSF_Loose_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_DATA").c_str());

    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesMC/";
    _muSF_Loose_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muSF_ISO_BCDEF->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_MC").c_str());
 

    // GH 
    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesDATA/";
    _muSF_Loose_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_DATA").c_str());

    filePath = "tkLooseISO_highptID_newpt_eta/efficienciesMC/";
    _muSF_Loose_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin0_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin1_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin2_&_HighPt_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muSF_ISO_GH->Get((filePath + "pair_newTuneP_probe_pt_PLOT_abseta_bin3_&_HighPt_pass_MC").c_str());
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
    float weight = 1;

    if (_isRealData)
        return weight;

    weight = _puReweight->Eval(nPU);

    return weight;
}




//////////////////
//   TRIGGERS   //
//////////////////


//--- SINGLE LEPTON ---//

EfficiencyContainer WeightUtils::GetTriggerEff(string triggerName, TLorentzVector &lepton) const
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



//--- DOUBLE LEPTON ---//

// Muons
// FIXME this is fake
EfficiencyContainer WeightUtils::GetDoubleMuonTriggerEff(TLorentzVector &muon, int leg) const
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
            etaBin = i;
            break;
        }
    }

    if (leg == 1)
    {
        effData = _eff_doubleMu_leg1_DATA[etaBin]->Eval(muon.Pt());
        effMC   = _eff_doubleMu_leg1_MC[etaBin]->Eval(muon.Pt());
    }
    else if (leg == 2)
    {
        effData = _eff_doubleMu_leg2_DATA[etaBin]->Eval(muon.Pt());
        effMC   = _eff_doubleMu_leg2_MC[etaBin]->Eval(muon.Pt());
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}


// Electrons
// FIXME this is fake
EfficiencyContainer WeightUtils::GetDoubleElectronTriggerEff(string triggerName, int leg, const baconhep::TElectron* electron) const
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }

    if (leg == 1)
    {
        effData = _eff_doubleEle_leg1_DATA->GetBinContent(_eff_doubleEle_leg1_DATA->FindBin(electron->scEta, electron->pt));
        effMC   = _eff_doubleEle_leg1_MC->GetBinContent(_eff_doubleEle_leg1_MC->FindBin(electron->scEta, electron->pt));
    }
    else if (leg == 2)
    {
        effData = _eff_doubleEle_leg2_DATA->GetBinContent(_eff_doubleEle_leg2_DATA->FindBin(electron->scEta, electron->pt));
        effMC   = _eff_doubleEle_leg2_MC->GetBinContent(_eff_doubleEle_leg2_MC->FindBin(electron->scEta, electron->pt));
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}




//////////////////
//      ID      //
//////////////////

/*
//--- TIGHT MUON ---//

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
        effData = _muSF_Tight_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Tight_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Tight_ID_DATA_BCDEF[etaBin], muon.Pt());
        errData = _muSF_Tight_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Tight_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }

    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF_Tight_ID_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Tight_ID_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Tight_ID_DATA_GH[etaBin], muon.Pt());
        errData = _muSF_Tight_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Tight_ID_MC_GH[etaBin]->GetErrorY(ptBin);
    }
    
    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



//--- LOOSE MUON ---//

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
*/


//--- HZZ MUON ---//

// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2017#Muon_scale_factors
// Using Moriond 17
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



//--- HZZ ELECTRON ---//

// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2017#Electron_scale_factors
EfficiencyContainer WeightUtils::GetHZZElectronRecoEff(const baconhep::TElectron* electron) const 
{
    float effData = 1, errData = 0, effMC = 1, errMC = 0;

    if (_isRealData)
    {
        EfficiencyContainer effCont(effData, effMC, errData, errMC);
        return effCont;
    }


    // ID
    float maxPt = 500, maxEta = 2.5;
    float sfID = 1, errID = 0;
    int bin;
    if (fabs(electron->scEta) < maxEta)
    {   
        if (kFALSE) // (electron->isGap)
        {
            if (electron->calibPt > maxPt)
                bin = _hzz_eleIdSF_gap->FindBin(electron->scEta, 0.99 * maxPt);
            else
                bin = _hzz_eleIdSF_gap->FindBin(electron->scEta, electron->calibPt);

            sfID = _hzz_eleIdSF_gap->GetBinContent(bin);
            errID = _hzz_eleIdSF_gap->GetBinError(bin);
        }
        else
        {
            if (electron->calibPt > maxPt)
                bin = _hzz_eleIdSF->FindBin(electron->scEta, 0.99 * maxPt);
            else
                bin = _hzz_eleIdSF->FindBin(electron->scEta, electron->calibPt);

            sfID = _hzz_eleIdSF->GetBinContent(bin);
            errID = _hzz_eleIdSF->GetBinError(bin);
        }
    }

    // Reco
    float minPt = 25, lowPt = 20, highPt = 75;
    float sfReco = 1, errReco = 0;
    if (fabs(electron->scEta) < maxEta)
    {
        if (electron->calibPt > maxPt)
            bin = _eleSF_RECO->FindBin(electron->scEta, 0.99 * maxPt);
        else if (electron->calibPt < minPt)
            bin = _eleSF_RECO->FindBin(electron->scEta, 1.01 * minPt);
        else
            bin = _eleSF_RECO->FindBin(electron->scEta, electron->calibPt);

        sfReco = _eleSF_RECO->GetBinContent(bin);
        errReco = _eleSF_RECO->GetBinError(bin);

        // Additional 1% error
        if (electron->calibPt < lowPt || electron->calibPt > highPt)
            errReco += sfReco * 0.01;
    }

    effData = sfID * sfReco;
    errData = sqrt(errID * errID + errReco * errReco);

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



/*
//////////////////
//     ISO      //
//////////////////


//--- TIGHT MUON ---//

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
        effData = _muSF_Tight_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Tight_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Tight_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        errData = _muSF_Tight_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Tight_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }
    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF_Tight_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Tight_ISO_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Tight_ISO_DATA_GH[etaBin], muon.Pt());
        errData = _muSF_Tight_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Tight_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}



//--- LOOSE MUON ---//

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
        effData = _muSF_Loose_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Loose_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Loose_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        errData = _muSF_Loose_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Loose_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);
    }
    else if (_dataPeriod == "2016GH")
    {
        effData = _muSF_Loose_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC   = _muSF_Loose_ISO_MC_GH[etaBin]->Eval(muon.Pt());

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_Loose_ISO_DATA_GH[etaBin], muon.Pt());
        errData = _muSF_Loose_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC   = _muSF_Loose_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}
*/



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
