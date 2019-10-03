#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
	_dataPeriod = dataPeriod;
	_selection  = selection;
	_isRealData = isRealData;

	rng = new TRandom3();


	//const std::string _cmssw_base = getenv("CMSSW_BASE");
	_cmssw_base = getenv("CMSSW_BASE");
    
	
	if(_dataPeriod == "2016ReReco")
		PeriodFolder = "ReReco2016";
	else if (dataPeriod == "2016Legacy")
		PeriodFolder = "Legacy2016";
	else if (dataPeriod == "2017Rereco")
		PeriodFolder = "ReReco2017";
	


	//SetTriggerWeights(PeriodFolder);
	SetMuonIDWeights     (PeriodFolder);
	SetMuonISOWeights    (PeriodFolder);



	PeriodFolder = "ReReco2016";
	SetMuonTriggerWeights(PeriodFolder);

	SetElectronTriggerWeights(PeriodFolder);	
	SetElectronIDWeights(PeriodFolder);	
	SetElectronISOWeights(PeriodFolder);	

	SetPhotonIDWeights(PeriodFolder);

	// PU weights
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/pileup_sf_2016_full.root";
	TFile* puFile = new TFile(_fileName.c_str(), "OPEN");
	_puReweight = (TGraph*)puFile->Get("pileup_sf");

	// photon r9 reweighting
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/photon_corrections/photon_r9_reweighting_2016.root";
	TFile* f_photon_r9 = new TFile(_fileName.c_str(), "OPEN"); 
	_photon_r9_barrel = (TGraph *)f_photon_r9->Get("transffull5x5R9EB");
	_photon_r9_endcap = (TGraph *)f_photon_r9->Get("transffull5x5R9EE");

	// photon r9 reweighting
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/photon_corrections/transformation_pho_presel_BDTUpto6000.root";
	TFile* f_photon_shower = new TFile(_fileName.c_str(), "OPEN"); 
	_photon_etawidth_barrel = (TGraph *)f_photon_shower->Get("transfEtaWidthEB");
	_photon_etawidth_endcap = (TGraph *)f_photon_shower->Get("transfEtaWidthEE");
	_photon_phiwidth_barrel = (TGraph *)f_photon_shower->Get("transfPhiWidthEB");
	_photon_phiwidth_endcap = (TGraph *)f_photon_shower->Get("transfPhiWidthEE");

	_photon_sieie_barrel    = (TGraph *)f_photon_shower->Get("transfSieieEB");
	_photon_sieie_endcap    = (TGraph *)f_photon_shower->Get("transfSieieEE");

	_photon_sieip_barrel    = (TGraph *)f_photon_shower->Get("transfSieipEB");
	_photon_sieip_endcap    = (TGraph *)f_photon_shower->Get("transfSieipEE");

	_photon_s4_barrel       = (TGraph *)f_photon_shower->Get("transfS4EB");
	_photon_s4_endcap       = (TGraph *)f_photon_shower->Get("transfS4EE");
}

void WeightUtils::SetMuonTriggerWeights(std::string PeriodFolder ){
	TFile* triggerFile;
	// muon trigger efficiencies 
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_trigger/EfficienciesAndSF_BCDEF.root";
	triggerFile = new TFile(_fileName.c_str(), "OPEN");

	_filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
	_eff_IsoMu24_DATA[0] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[1] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[2] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[3] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());


	triggerFile = new TFile(_fileName.c_str(), "OPEN");

	_filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
	_eff_IsoMu24_DATA[0] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[1] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[2] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
	_eff_IsoMu24_DATA[3] = (TGraphAsymmErrors*)triggerFile->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());


	// double muon trigger efficiencies

	// leg 1
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu17Leg_Eta0to09.root";
	TFile* f_DoubleMuTrigSF_leg1_0 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg1_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_data");
	_eff_doubleMu_leg1_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_mc"); 
	_sf_doubleMu_leg1[0] = (TH1F *)f_DoubleMuTrigSF_leg1_0->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu17Leg_Eta09to12.root";
	TFile* f_DoubleMuTrigSF_leg1_1 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg1_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_data");
	_eff_doubleMu_leg1_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_mc"); 
	_sf_doubleMu_leg1[1] = (TH1F *)f_DoubleMuTrigSF_leg1_1->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu17Leg_Eta12to21.root";
	TFile* f_DoubleMuTrigSF_leg1_2 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg1_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_data");
	_eff_doubleMu_leg1_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_mc"); 
	_sf_doubleMu_leg1[2] = (TH1F *)f_DoubleMuTrigSF_leg1_2->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu17Leg_Eta21to24.root";
	TFile* f_DoubleMuTrigSF_leg1_3 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg1_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_data");
	_eff_doubleMu_leg1_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_mc"); 
	_sf_doubleMu_leg1[3] = (TH1F *)f_DoubleMuTrigSF_leg1_3->Get("scale_factor");

	// leg 2 
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu8Leg_Eta0to09.root";
	TFile* f_DoubleMuTrigSF_leg2_0 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg2_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_data");
	_eff_doubleMu_leg2_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_mc"); 
	_sf_doubleMu_leg2[0] = (TH1F *)f_DoubleMuTrigSF_leg2_0->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu8Leg_Eta09to12.root";
	TFile* f_DoubleMuTrigSF_leg2_1 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg2_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_data");
	_eff_doubleMu_leg2_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_mc"); 
	_sf_doubleMu_leg2[1] = (TH1F *)f_DoubleMuTrigSF_leg2_1->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu8Leg_Eta12to21.root";
	TFile* f_DoubleMuTrigSF_leg2_2 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg2_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_data");
	_eff_doubleMu_leg2_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_mc"); 
	_sf_doubleMu_leg2[2] = (TH1F *)f_DoubleMuTrigSF_leg2_2->Get("scale_factor");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/doublemuon_trigger/sf_Mu8Leg_Eta21to24.root";
	TFile* f_DoubleMuTrigSF_leg2_3 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleMu_leg2_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_data");
	_eff_doubleMu_leg2_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_mc"); 
	_sf_doubleMu_leg2[3] = (TH1F *)f_DoubleMuTrigSF_leg2_3->Get("scale_factor");

}


void WeightUtils::SetMuonIDWeights(std::string PeriodFolder){

	if(_dataPeriod== "2016ReReco"){
		// muon tight ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
		_muSF_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
		_muSF_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
		_muSF_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
		_muSF_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

		_filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
		_muSF_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_MC").c_str());
		_muSF_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_MC").c_str());
		_muSF_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_MC").c_str());
		_muSF_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_MC").c_str());

		// muon tight ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
		_muSF_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
		_muSF_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
		_muSF_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
		_muSF_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

		_filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
		_muSF_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin0_MC").c_str());
		_muSF_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin1_MC").c_str());
		_muSF_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin2_MC").c_str());
		_muSF_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin3_MC").c_str());

		// muon loose ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_Loose_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
		_muSF_Loose_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
		_muSF_Loose_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
		_muSF_Loose_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
		_muSF_Loose_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

		_filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
		_muSF_Loose_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_MC").c_str());
		_muSF_Loose_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_MC").c_str());
		_muSF_Loose_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_MC").c_str());
		_muSF_Loose_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_MC").c_str());

		// muon loose ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_Loose_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
		_muSF_Loose_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
		_muSF_Loose_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
		_muSF_Loose_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
		_muSF_Loose_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

		_filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
		_muSF_Loose_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin0_MC").c_str());
		_muSF_Loose_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin1_MC").c_str());
		_muSF_Loose_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin2_MC").c_str());
		_muSF_Loose_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((_filePath + "pt_PLOT_abseta_bin3_MC").c_str());

		// hzz muon id efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/hzz_muon_id_sf.root";
		TFile* f_hzz_muIdSF = new TFile(_fileName.c_str(), "OPEN");
		_hzz_muIdSF = (TH2F*)f_hzz_muIdSF->Get("FINAL");

	}
	else if(_dataPeriod == "2016Legacy"){
		// muon tight ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ID_BCDEF = (TH2F*)f_muRecoSF_ID_BCDEF->Get("NUM_TightID_DEN_genTracks_eta_pt");

		// muon tight ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ID_GH = (TH2F*)f_muRecoSF_ID_GH->Get("NUM_TightID_DEN_genTracks_eta_pt");


		
		// muon loose ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_Loose_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ID_BCDEF = (TH2F*)f_muRecoSF_Loose_ID_BCDEF->Get("NUM_LooseID_DEN_genTracks_eta_pt");

		// muon loose ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_Loose_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ID_GH = (TH2F*)f_muRecoSF_Loose_ID_GH->Get("NUM_LooseID_DEN_genTracks_eta_pt");
	}
	else if(_dataPeriod == "2017ReReco"){
		// muon tight ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ID_BCDEF = (TH2F*)f_muRecoSF_ID_BCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta");

		// muon tight ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ID_GH = (TH2F*)f_muRecoSF_ID_GH->Get("NUM_TightID_DEN_genTracks_pt_abseta");


		
		// muon loose ID sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_BCDEF.root";
		TFile* f_muRecoSF_Loose_ID_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ID_BCDEF = (TH2F*)f_muRecoSF_Loose_ID_BCDEF->Get("NUM_LooseID_DEN_genTracks_pt_abseta");

		// muon loose ID sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_id/EfficienciesAndSF_GH.root";
		TFile* f_muRecoSF_Loose_ID_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ID_GH = (TH2F*)f_muRecoSF_Loose_ID_GH->Get("NUM_LooseID_DEN_genTracks_pt_abseta");
	}
}
void WeightUtils::SetMuonISOWeights(std::string PeriodFolder){

	if(_dataPeriod == "2016ReReco"){
		// tight muon ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
		_muSF_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

		_filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
		_muSF_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

		// tight muon ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
		_muSF_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
		_muSF_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

		_filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
		_muSF_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
		_muSF_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

		// loose muon ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_Loose_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "LooseISO_LooseID_pt_eta/efficienciesDATA/";
		_muSF_Loose_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_&_PF_pass_DATA").c_str());

		_filePath = "LooseISO_LooseID_pt_eta/efficienciesMC/";
		_muSF_Loose_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin1_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin2_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((_filePath + "pt_PLOT_abseta_bin3_&_PF_pass_MC").c_str());

		// loose muon ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_Loose_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_filePath = "LooseISO_LooseID_pt_eta/efficienciesDATA/";
		_muSF_Loose_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
		_muSF_Loose_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());

		_filePath = "LooseISO_LooseID_pt_eta/efficienciesMC/";
		_muSF_Loose_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
		_muSF_Loose_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((_filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
	}
	else if(_dataPeriod == "2016Legacy"){
		// muon tight ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ISO_BCDEF = (TH2F*)f_muRecoSF_ISO_BCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
		
		// muon tight ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ISO_GH = (TH2F*)f_muRecoSF_ISO_GH->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");



		// muon loose ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_Loose_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ISO_BCDEF = (TH2F*)f_muRecoSF_Loose_ISO_BCDEF->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt");
		
		// muon loose ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_Loose_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ISO_GH = (TH2F*)f_muRecoSF_Loose_ISO_GH->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt");
	}
	else if(_dataPeriod == "2017ReReco"){
		// muon tight ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ISO_BCDEF = (TH2F*)f_muRecoSF_ISO_BCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
		
		// muon tight ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_tight_ISO_GH = (TH2F*)f_muRecoSF_ISO_GH->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");



		// muon loose ISO sf (BCDEF)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_BCDEF.root";
		 TFile* f_muRecoSF_Loose_ISO_BCDEF = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ISO_BCDEF = (TH2F*)f_muRecoSF_Loose_ISO_BCDEF->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt");
		
		// muon loose ISO sf (GH)
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/muon_iso/EfficienciesAndSF_GH.root";
		 TFile* f_muRecoSF_Loose_ISO_GH = new TFile(_fileName.c_str(), "OPEN"); 

		_muSF_loose_ISO_GH = (TH2F*)f_muRecoSF_Loose_ISO_GH->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt");
		
	}
}


void WeightUtils::SetElectronTriggerWeights(std::string PeriodFolder ){

	// double electron trigger efficiencies
	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/doubleg_trigger/SFs_Leg1_Ele23_HZZSelection_Tag35.root";
	TFile* f_elTrigSF_leg1 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleg_leg1_DATA = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffData2D");
	_eff_doubleg_leg1_MC   = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffMC2D"); 
	_sf_doubleg_leg1 = (TH2F *)f_elTrigSF_leg1->Get("EGamma_SF2D");

	_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/doubleg_trigger/SFs_Leg2_Ele12_HZZSelection_Tag35.root";
	TFile* f_elTrigSF_leg2 = new TFile(_fileName.c_str(), "OPEN");
	_eff_doubleg_leg2_DATA = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffData2D");
	_eff_doubleg_leg2_MC   = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffMC2D");
	_sf_doubleg_leg2 = (TH2F *)f_elTrigSF_leg2->Get("EGamma_SF2D");
}
void WeightUtils::SetElectronRECOWeights(std::string PeriodFolder){
	if(PeriodFolder == "ReReco2016"){
		// electron reco efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/electron_reco/egamma_eff_reco_2016.root";
		TFile* f_eleRecoSF = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");
		_eleSF_RECO_2D = (TH2F *)f_eleRecoSF->Get("EGamma_SF2D");
	} else if(PeriodFolder == "Legacy2016"){
		// electron reco efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/electron_reco/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root";
		TFile* f_eleRecoSF = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");
		_eleSF_RECO_2D = (TH2F *)f_eleRecoSF->Get("EGamma_SF2D");
	} else if(PeriodFolder == "ReReco2017"){
		// electron reco efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/electron_reco/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root";
		TFile* f_eleRecoSF = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");
		_eleSF_RECO_2D = (TH2F *)f_eleRecoSF->Get("EGamma_SF2D");
	}

}
void WeightUtils::SetElectronIDWeights(std::string PeriodFolder ){
	if(PeriodFolder == "ReReco2016"){
		// electron id efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/egamma_eff_ID_2016.root";
		TFile* f_eleIdSF = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
		_eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
		_eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
		_eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
		_eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");

		// hzz electron id efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/electron_id/egamma_eff_hzz_ID_2016.root";
		TFile* f_hzz_eleIdSF = new TFile(_fileName.c_str(), "OPEN"); 
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

		_hzz_eleSF_ID_2D = (TH2F *)f_hzz_eleIdSF->Get("EGamma_SF2D");
	} else if(PeriodFolder == "ReReco2017"){
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root";
		TFile* f_eleIdSF_tight = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_RECO_2D = (TH2F *)f_eleIdSF_tight->Get("EGamma_SF2D");
	} else if(PeriodFolder == "Legacy2016"){
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/2016LegacyReReco_ElectronTight.root";
		TFile* f_eleIdSF_tight = new TFile(_fileName.c_str(), "OPEN"); 
		_eleSF_RECO_2D = (TH2F *)f_eleIdSF_tight->Get("EGamma_SF2D");
	}

}
void WeightUtils::SetElectronISOWeights(std::string PeriodFolder ){}

void WeightUtils::SetPhotonIDWeights(std::string PeriodFolder ){
	if(PeriodFolder == "ReReco2016"){
		// photon mva id (90%) efficiencies
		_fileName = _cmssw_base + "/src/BLT/BLTAnalysis/data/"+ PeriodFolder+"/photon_id/photon_mva_id_2016.root";
		TFile* f_mva_gammaIdSF = new TFile(_fileName.c_str(), "OPEN");
		_mva_gammaSF_ID[0] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_0");
		_mva_gammaSF_ID[1] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_1");
		_mva_gammaSF_ID[2] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_2");
		_mva_gammaSF_ID[3] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_3");
		_mva_gammaSF = (TH2F *)f_mva_gammaIdSF->Get("EGamma_SF2D");
	}
	else if(PeriodFolder == "Legacy2016"){
	}
	else if(PeriodFolder == "ReReco2017"){
	}
}
void WeightUtils::SetPhotonISOWeights(std::string PeriodFolder ){}

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

float WeightUtils::GetPUWeight(float nPU)
{
    return _puReweight->Eval(nPU); 
}

std::pair<float,float> WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(lepton.Eta()) > binningEta[i] && fabs(lepton.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }
    
    float effMC   = 1;
    float effData = 1;
    if (triggerName == "HLT_IsoMu24_v*") {
        if (lepton.Pt() < 500.) {
            effData = _eff_IsoMu24_DATA[etaBin]->Eval(lepton.Pt());
        }
    }
    
    return std::make_pair(effData, effMC);
}

/*std::pair<float,float> WeightUtils::GetDoubleEGTriggerEffWeight(string triggerName, TElectron &electron) const
{
    float effData = 1;
    float effMC = 1;

    if (electron.calibPt < 200.) {
        if (triggerName == "HLT_DoubleEG_leg1") {
            //effData = _eff_doubleg_leg1_DATA->Interpolate(lepton.Eta(), lepton.Pt());
            //effMC   = _eff_doubleg_leg1_MC->Interpolate(lepton.Eta(), lepton.Pt());
            effData = _eff_doubleg_leg1_DATA->GetBinContent(_eff_doubleg_leg1_DATA->FindBin(electron.scEta, electron.calibPt));
            effMC = _eff_doubleg_leg1_MC->GetBinContent(_eff_doubleg_leg1_MC->FindBin(electron.scEta, electron.calibPt));
        }
        else if (triggerName == "HLT_DoubleEG_leg2") {
            //effData = _eff_doubleg_leg2_DATA->Interpolate(lepton.Eta(), lepton.Pt());
            //effMC   = _eff_doubleg_leg2_MC->Interpolate(lepton.Eta(), lepton.Pt());
            effData = _eff_doubleg_leg2_DATA->GetBinContent(_eff_doubleg_leg2_DATA->FindBin(electron.scEta, electron.calibPt));
            effMC = _eff_doubleg_leg2_MC->GetBinContent(_eff_doubleg_leg2_MC->FindBin(electron.scEta, electron.calibPt));
        }
    }

    if (effMC == 0) {
        cout << "zero value for effMC" << endl;
        effMC = 1;
    }

    return std::make_pair(effData, effMC);
}*/

float WeightUtils::GetDoubleEGTriggerEffWeight(string triggerName, TElectron &electron) const
{
    float weight = 1.;
    float tmpElePt = (electron.pt < 200.) ? electron.pt : 199.;

    if (triggerName == "HLT_DoubleEG_leg1") 
        weight *= _sf_doubleg_leg1->GetBinContent(_sf_doubleg_leg1->FindBin(electron.scEta, tmpElePt));
    else if (triggerName == "HLT_DoubleEG_leg2") 
        weight *= _sf_doubleg_leg2->GetBinContent(_sf_doubleg_leg2->FindBin(electron.scEta, tmpElePt));

    return weight;
}

/*std::pair<float,float> WeightUtils::GetDoubleMuonTriggerEffWeight(string triggerName, TMuon &muon) const
{
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;
    float effData = 1;
    float effMC = 1;

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.eta) > binningEta[i] && fabs(muon.eta) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (triggerName == "HLT_DoubleMuon_leg1") {
        effData = _eff_doubleMu_leg1_DATA[etaBin]->Eval(muon.pt);
        effMC = _eff_doubleMu_leg1_MC[etaBin]->Eval(muon.pt);
        //effData = _eff_doubleMu_leg1_DATA[etaBin]->FindBin(tmpMuPt);
        //effMC = _eff_doubleMu_leg1_MC[etaBin]->FindBin(tmpMuPt);
    }
    else if (triggerName == "HLT_DoubleMuon_leg2") {
        effData = _eff_doubleMu_leg2_DATA[etaBin]->Eval(muon.pt);
        effMC = _eff_doubleMu_leg2_MC[etaBin]->Eval(muon.pt);
        //effData = _eff_doubleMu_leg2_DATA[etaBin]->FindBin(tmpMuPt);
        //effMC = _eff_doubleMu_leg2_MC[etaBin]->FindBin(tmpMuPt);
    }

    if (effMC == 0) {
        cout << "zero value for effMC" << endl;
        effMC = 1;
    }

    return std::make_pair(effData, effMC);
}*/

float WeightUtils::GetDoubleMuonTriggerEffWeight(string triggerName, TMuon &muon) const
{
    float weight = 1.;
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.eta) > binningEta[i] && fabs(muon.eta) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (triggerName == "HLT_DoubleMuon_leg1") 
        weight *= _sf_doubleMu_leg1[etaBin]->GetBinContent(_sf_doubleMu_leg1[etaBin]->FindBin(tmpMuPt));
    else if (triggerName == "HLT_DoubleMuon_leg2") 
        weight *= _sf_doubleMu_leg2[etaBin]->GetBinContent(_sf_doubleMu_leg2[etaBin]->FindBin(tmpMuPt));

    return weight;
}

float WeightUtils::GetMuonIDEff(TLorentzVector& muon) const
{
	float weight = 1;
	if(_dataPeriod == "2016ReReco"){
	    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
	    int etaBin = 0;
	    for (int i = 0; i < 4; ++i) {
		if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
		    etaBin = i;
		    break;
		}
	    }
	    float random = rng->Rndm();
	    if (muon.Pt() < 200.) {
		if (random > 0.468) {
		    weight   *= _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
		} else {
		    weight   *= _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());
		}
	    }
	}
	else if(_dataPeriod == "2016Legacy"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_tight_ID_BCDEF->GetBinContent(_muSF_tight_ID_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_tight_ID_GH->GetBinContent(_muSF_tight_ID_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }
	}
	else if(_dataPeriod == "2017ReReco"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_tight_ID_BCDEF->GetBinContent(_muSF_tight_ID_BCDEF->FindBin(muon.Pt(),fabs(muon.Eta()) ));
		} else {
		    weight   *= _muSF_tight_ID_GH->GetBinContent(_muSF_tight_ID_GH->FindBin(muon.Pt(),fabs(muon.Eta())));
		}
	    }
	}

	return weight;
}

float WeightUtils::GetLooseMuonIDEff(TLorentzVector& muon) const
{

	float weight = 1;
	if(_dataPeriod == "2016ReReco"){
		float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
		int etaBin = 0;
		for (int i = 0; i < 4; ++i) {
			if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
			    etaBin = i;
			    break;
			}
		}

		float random = rng->Rndm();
		if (muon.Pt() < 200.) {
			if (random > 0.468) {
			    weight   *= _muSF_Loose_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_Loose_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
			} else {
			    weight   *= _muSF_Loose_ID_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_Loose_ID_MC_GH[etaBin]->Eval(muon.Pt());
			}
		}
	}
	else if(_dataPeriod == "2016Legacy"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_loose_ID_BCDEF->GetBinContent(_muSF_loose_ID_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_loose_ID_GH->GetBinContent(_muSF_loose_ID_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }
	}
	else if(_dataPeriod == "2017ReReco"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_loose_ID_BCDEF->GetBinContent(_muSF_loose_ID_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_loose_ID_GH->GetBinContent(_muSF_loose_ID_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }
	}

	return weight;
}

float WeightUtils::GetMuonISOEff(TLorentzVector& muon) const
{

	float weight = 1;
	if(_dataPeriod == "2016ReReco"){
	    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
	    int etaBin = 0;
	    for (int i = 0; i < 4; ++i) {
		if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
		    etaBin = i;
		    break;
		}
	    }

	    float random = rng->Rndm();
	    if (muon.Pt() < 200.) {
		if (random > 0.468) {
		    weight   *= _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
		} else {
		    weight   *= _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());
		}
	    }
	}
	else if(_dataPeriod == "2016Legacy"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_tight_ISO_BCDEF->GetBinContent(_muSF_tight_ISO_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_tight_ISO_GH->GetBinContent(_muSF_tight_ISO_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }

	}
	else if(_dataPeriod == "2017ReReco"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_tight_ISO_BCDEF->GetBinContent(_muSF_tight_ISO_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_tight_ISO_GH->GetBinContent(_muSF_tight_ISO_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }

	}
	return weight;
}

float WeightUtils::GetLooseMuonISOEff(TLorentzVector& muon) const
{

	float weight = 1;
	if(_dataPeriod == "2016ReReco"){
		float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
		int etaBin = 0;
		for (int i = 0; i < 4; ++i) {
			if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
			    etaBin = i;
			    break;
			}
		}

		float random = rng->Rndm();
		if (muon.Pt() < 200.) {
			if (random > 0.468) {
			    weight   *= _muSF_Loose_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_Loose_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
			} else {
			    weight   *= _muSF_Loose_ISO_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_Loose_ISO_MC_GH[etaBin]->Eval(muon.Pt());
			}
		}
	}
	else if(_dataPeriod == "2016Legacy"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_loose_ISO_BCDEF->GetBinContent(_muSF_loose_ISO_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_loose_ISO_GH->GetBinContent(_muSF_loose_ISO_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }

	}
	else if(_dataPeriod == "2017ReReco"){
	    float random = rng->Rndm();
	    if (muon.Pt() < 120.) {
		if (random > 0.468) {
		    weight   *= _muSF_loose_ISO_BCDEF->GetBinContent(_muSF_loose_ISO_BCDEF->FindBin(fabs(muon.Eta()), muon.Pt()));
		} else {
		    weight   *= _muSF_loose_ISO_GH->GetBinContent(_muSF_loose_ISO_GH->FindBin(fabs(muon.Eta()), muon.Pt()));
		}
	    }

	}

	return weight;
}

float WeightUtils::GetHZZMuonIDEff(TMuon& muon) const
{
    float weight = 1;
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;
    
    //weight *= _hzz_muIdSF->Interpolate(muon.eta, muon.pt);
    weight *= _hzz_muIdSF->GetBinContent(_hzz_muIdSF->FindBin(muon.eta, tmpMuPt));
    
    return weight;
}



float WeightUtils::GetElectronRecoIdEff(TLorentzVector& electron) const
{
    float binningPt[] = {10., 20., 35., 50., 90., 500.};
    int ptBin = 0;
    for (int i = 0; i < 5; ++i) {
        if (fabs(electron.Pt()) > binningPt[i] && fabs(electron.Pt()) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }

    float weight = 1;
    if (electron.Pt() < 500.) {
        weight *= _eleSF_RECO->Eval(electron.Eta());
        weight *= _eleSF_ID[ptBin]->Eval(electron.Eta());
    }
    
    return weight;
}

float WeightUtils::GetHZZElectronRecoIdEff(TElectron& electron) const 
{
    //float tmpElePt = (electron.calibPt < 200.) ? electron.calibPt : 199.;
    float tmpElePt = (electron.pt < 200.) ? electron.pt : 199.;
    /*float binningPt[] = {7., 15., 20., 30., 40., 50., 60., 70., 80., 100., 120., 140., 160., 200.}; 
    int ptBin = 0;
    for (int i = 0; i < 13; ++i) {
        if (tmpElePt > binningPt[i] && tmpElePt <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }*/

    float weight = 1;
    weight *= _eleSF_RECO_2D->GetBinContent(_eleSF_RECO_2D->FindBin(electron.scEta, 50.));
    weight *= _hzz_eleSF_ID_2D->GetBinContent(_hzz_eleSF_ID_2D->FindBin(electron.scEta, tmpElePt));
    //weight *= _eleSF_RECO->Eval(electron.scEta);
    //weight *= _hzz_eleSF_ID[ptBin]->Eval(electron.scEta);
    
    return weight;
}

float WeightUtils::GetPhotonMVAIdEff(TPhoton& photon) const
{
    /*float binningPt[] = {20., 35., 50., 90., 150.};
    int ptBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(photon.calibPt) > binningPt[i] && fabs(photon.calibPt) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }*/
    float tmpPhotonPt = (photon.pt < 150.) ? photon.pt : 149.;
    float weight = 1;
    //if (photon.calibPt < 150.) {
    //    weight *= _mva_gammaSF_ID[ptBin]->Eval(photon.scEta);
    //}
    
    //weight *= _mva_gammaSF->GetBinContent(_mva_gammaSF->FindBin(photon.scEta, photon.pt));
    weight *= _mva_gammaSF->GetBinContent(_mva_gammaSF->FindBin(photon.scEta, tmpPhotonPt));

    // electron veto scale factor
    if (fabs(photon.scEta) <= 1.49)
        weight *= 0.9938;
    else if (fabs(photon.scEta) > 1.49)
        weight *= 0.9875;

    if (weight == 0)
	weight = 1.; 
    return weight;
}
/////////////////////////////////////////////////////
//
//     Correction Functions 
float WeightUtils::GetCorrectedPhotonR9(TPhoton& photon) const 
{
    float r9 = photon.r9;
    if (fabs(photon.scEta) < 1.444)
        r9 = _photon_r9_barrel->Eval(photon.r9);
    else if (fabs(photon.scEta) > 1.566)
        r9 = _photon_r9_endcap->Eval(photon.r9);
    else
        std::cout << "bad value of photon scEta: returning original r9" << std::endl;
    
    return r9;
}
float WeightUtils::GetCorrectedPhotonEtaWidth(TPhoton& photon) const 
{
    float etawidth = photon.scEtaWidth;
    if (fabs(photon.scEta) < 1.444)
        etawidth = _photon_etawidth_barrel->Eval(photon.scEtaWidth);
    else if (fabs(photon.scEta) > 1.566)
        etawidth = _photon_etawidth_endcap->Eval(photon.scEtaWidth);
    else
        std::cout << "bad value of photon scEta: returning original etawidth" << std::endl;
    
    return etawidth;
}
float WeightUtils::GetCorrectedPhotonPhiWidth(TPhoton& photon) const 
{
    float phiwidth = photon.scPhiWidth;
    if (fabs(photon.scEta) < 1.444)
        phiwidth = _photon_phiwidth_barrel->Eval(photon.scPhiWidth);
    else if (fabs(photon.scEta) > 1.566)
        phiwidth = _photon_phiwidth_endcap->Eval(photon.scPhiWidth);
    else
        std::cout << "bad value of photon scEta: returning original phiwidth" << std::endl;
    
    return phiwidth;
}
float WeightUtils::GetCorrectedPhotonSieie(TPhoton& photon) const 
{
    float sieie = photon.sieie;
    if (fabs(photon.scEta) < 1.444)
        sieie = _photon_sieie_barrel->Eval(photon.sieie);
    else if (fabs(photon.scEta) > 1.566)
        sieie = _photon_sieie_endcap->Eval(photon.sieie);
    else
        std::cout << "bad value of photon scEta: returning original sieie" << std::endl;
    
    return sieie;
}
float WeightUtils::GetCorrectedPhotonSieip(TPhoton& photon) const 
{
    float sieip = photon.sieip;
    if (fabs(photon.scEta) < 1.444)
        sieip = _photon_sieip_barrel->Eval(photon.sieip);
    else if (fabs(photon.scEta) > 1.566)
        sieip = _photon_sieip_endcap->Eval(photon.sieip);
    else
        std::cout << "bad value of photon scEta: returning original sieip" << std::endl;
    
    return sieip;
}
float WeightUtils::GetCorrectedPhotonS4(TPhoton& photon) const 
{
    float s4 = photon.e2x2/photon.e5x5;
    if (fabs(photon.scEta) < 1.444)
        s4 = _photon_s4_barrel->Eval(s4);
    else if (fabs(photon.scEta) > 1.566)
        s4 = _photon_s4_endcap->Eval(s4);
    else
        std::cout << "bad value of photon scEta: returning original s4" << std::endl;
    
    return s4;
}
float WeightUtils::GetCorrectedPhotonRho(TPhoton& photon, TEventInfo& Info) const 
{
    float rho = Info.rhoJet;
    if (fabs(photon.scEta) < 1.444)
        rho = _photon_rho_barrel->Eval(Info.rhoJet);
    else if (fabs(photon.scEta) > 1.566)
        rho = _photon_rho_endcap->Eval(Info.rhoJet);
    else
        std::cout << "bad value of photon scEta: returning original rho" << std::endl;
    
    return rho;
}
