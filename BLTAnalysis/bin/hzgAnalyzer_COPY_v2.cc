#include "hzgAnalyzer.hh"
#include <map>
#include <fstream>
#include <math.h>

#include <TSystem.h>
#include <TF2.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}

hzgAnalyzer::hzgAnalyzer() : BLTSelector()
{

}

hzgAnalyzer::~hzgAnalyzer()
{

}

void hzgAnalyzer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
            std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));
    nTrig = 0;
    noPass = 0;
    if (params->selection == "mugmjetjet"){
	// SINGLE ELECTRON TRIGGERS
	//triggerNames.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*");

	/// SINGLE MUON TRIGGERS
	//triggerNames.push_back("HLT_IsoTkMu24_v*");
	triggerNames.push_back("HLT_IsoMu24_v*");
	
	/// DOUBLE MUON TRIGGERS
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
 	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
	
	/// MUON + GAMMA TRIGGER
	//triggerNames.push_back("HLT_Mu17_Photon22_CaloIdL_L1ISO_v*");

	/// JET TRIGGERS
	//triggerNames.push_back("HLT_BTagMu_DiJet20_Mu5_v*");
	//triggerNames.push_back("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v*");
    }
    if (params->selection == "mugm"){
	/// SINGLE MUON TRIGGERS
	//triggerNames.push_back("HLT_IsoTkMu24_v*");
	//triggerNames.push_back("HLT_IsoMu24_v*");
	
	/// DOUBLE MUON TRIGGERS
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
 	//triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
	
	/// MUON + GAMMA TRIGGER
	//triggerNames.push_back("HLT_Mu17_Photon22_CaloIdL_L1ISO_v*");

	/// JET TRIGGERS
	//triggerNames.push_back("HLT_BTagMu_DiJet20_Mu5_v*");
	//triggerNames.push_back("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v*");
    }
    if (params->selection == "mumu" || params->selection == "mumug") {
	//triggerNames.push_back("HLT_Mu17_Photon22_CaloIdL_L1ISO_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    }
    else if (params->selection == "ee" || params->selection == "elelg") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
    }
    else if (params->selection == "tautaug") { // select one muon plus one hadronic tau (for now)
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    }

    // Set the cuts
    cuts.reset(new Cuts()); 
    particleSelector.reset(new ParticleSelector(*params, *cuts));
    // Pre/Selection
    if(triggerNames[0] == "HLT_IsoMu24_v*"){
    	PreSel_pt_mu        = 23;
    	Sel_pt_mu1          = 25;
    	PreSel_pt_ph        = 0; 
    	Sel_pt_ph           = 20;//15;
    }
    else {
    	PreSel_pt_mu 	= 16;
	Sel_pt_mu1      = 18;
    	PreSel_pt_ph    = 21;
    	Sel_pt_ph       = 23;
    }
    Sel_pt_mu2          = 18;
    PreSel_eta_mu       = 2.4;


    PreSel_pt_el	= 22;
    PreSel_eta_el	= 2.4;
    Sel_pt_el1          = 24;

    PreSel_eta_ph       = 2.4;

    PreSel_pt_jet  	= 30;
    PreSel_eta_jet 	= 4.7;

    leptonTagCharge = 1;
    nJetPass_Sel  = 0;


    nPre_Mu = 0;
    nPre_El = 0;
    nPre_Jet = 0;
    nPre_Ph = 0;

    nSel_Mu = 0;
    nSel_El = 0;
    nSel_Jet = 0;
    nSel_Ph = 0;

    // Gen level counts 
    nGen = 0;
    nB =0;
    nPh = 0;
    nLep = 0; 

    nEl_elel =0; nEl_mumu =0; nEl_tautau =0; nEl_hh =0;
    nMu_elel =0; nMu_mumu =0; nMu_tautau =0; nMu_hh =0;
    nTau_elel =0; nTau_mumu =0; nTau_tautau =0; nTau_hh =0;
    nHad_elel =0; nHad_mumu =0; nHad_tautau =0; nHad_hh =0;

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");
    rng = new TRandom3();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");
    string outGenTreeName = params->get_output_treename("tree");
    outGenTreeName = "Gen" + outGenTreeName;

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");
    outGenTree = new TTree(outGenTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("xPV", &xPV);
    outTree->Branch("yPV", &yPV);
    outTree->Branch("zPV", &zPV);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("metNC", &metNC);
    outTree->Branch("metPhiNC", &metPhiNC);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);
    outTree->Branch("flagGen", &flagGen);

    // weights
    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("elIDWeightOne", &elIDWeightOne);
    outTree->Branch("elIDWeightTwo", &elIDWeightTwo);
    outTree->Branch("elTrigWeightOne", &elTrigWeightOne);
    outTree->Branch("elTrigWeightTwo", &elTrigWeightTwo);
    outTree->Branch("muonIDWeightOne", &muonIDWeightOne);
    outTree->Branch("muonIDWeightTwo", &muonIDWeightTwo);
    outTree->Branch("muonTrigWeightOne", &muonTrigWeightOne);
    outTree->Branch("muonTrigWeightTwo", &muonTrigWeightTwo);
    outTree->Branch("photonIDWeight", &photonIDWeight);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonTwoP4", &leptonTwoP4);

    outTree->Branch("leptonOnePt", &leptonOnePt);
    outTree->Branch("leptonOneEta", &leptonOneEta);
    outTree->Branch("leptonOnePhi", &leptonOnePhi);
    outTree->Branch("leptonOnePtKin", &leptonOnePtKin);
    //outTree->Branch("leptonOneP4KinFit", &leptonOneP4KinFit);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);

    outTree->Branch("leptonTwoPt", &leptonTwoPt);
    outTree->Branch("leptonTwoEta", &leptonTwoEta);
    outTree->Branch("leptonTwoPhi", &leptonTwoPhi);
    outTree->Branch("leptonTwoPtKin", &leptonTwoPtKin);
    //outTree->Branch("leptonTwoP4KinFit", &leptonTwoP4KinFit);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
 
    outTree->Branch("tauDecayMode", &tauDecayMode);
    outTree->Branch("tauMVA", &tauMVA);

    outTree->Branch("isLeptonTag", &isLeptonTag);
    outTree->Branch("isDijetTag", &isDijetTag);
    outTree->Branch("isTightDijetTag", &isTightDijetTag);

    // photons
    outTree->Branch("photonOneP4", &photonOneP4);
    outTree->Branch("photonOnePt", &photonOnePt);
    outTree->Branch("photonOneEta", &photonOneEta);
    outTree->Branch("photonOnePhi", &photonOnePhi);
    outTree->Branch("photonOneR9", &photonOneR9);
    outTree->Branch("photonOneMVA", &photonOneMVA);
    outTree->Branch("photonOneERes", &photonOneERes);
    outTree->Branch("passElectronVeto", &passElectronVeto);

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetTwoP4", &jetTwoP4);

    outTree->Branch("jetOnePt", &jetOnePt);
    outTree->Branch("jetOneEta", &jetOneEta);
    outTree->Branch("jetOnePhi", &jetOnePhi);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoPt", &jetTwoPt);
    outTree->Branch("jetTwoEta", &jetTwoEta);
    outTree->Branch("jetTwoPhi", &jetTwoPhi);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    // gen level objects 	

    outGenTree->Branch("genLeptonP4", &genLeptonP4);
    outGenTree->Branch("genPhotonP4", &genPhotonP4);
    outGenTree->Branch("genJetOneP4", &genJetOneP4);
    outGenTree->Branch("genJetTwoP4", &genJetTwoP4);

    outTree->Branch("genLeptonP4", &genLeptonP4);
    outTree->Branch("genPhotonP4", &genPhotonP4);
    outTree->Branch("genJetOneP4", &genJetOneP4);
    outTree->Branch("genJetTwoP4", &genJetTwoP4);

    outTree->Branch("genLeptonOnePt", &genLeptonOnePt);
    outTree->Branch("genLeptonOneEta", &genLeptonOneEta);
    outTree->Branch("genLeptonOnePhi", &genLeptonOnePhi);
    outTree->Branch("genLeptonOneId", &genLeptonOneId);
    outTree->Branch("genLeptonTwoPt", &genLeptonTwoPt);
    outTree->Branch("genLeptonTwoEta", &genLeptonTwoEta);
    outTree->Branch("genLeptonTwoPhi", &genLeptonTwoPhi);
    outTree->Branch("genLeptonTwoId", &genLeptonTwoId);
    outTree->Branch("genPhotonPt", &genPhotonPt);
    outTree->Branch("genPhotonEta", &genPhotonEta);
    outTree->Branch("genPhotonPhi", &genPhotonPhi);
    outTree->Branch("genPhotonFHPFS", &genPhotonFHPFS);
    outTree->Branch("genPhotonIPFS", &genPhotonIPFS);
    outTree->Branch("vetoDY", &vetoDY);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nCentralJets", &nCentralJets);
    outTree->Branch("nBJets", &nBJets);
    
    // dilepton
    outTree->Branch("dileptonPt", &dileptonPt);
    outTree->Branch("dileptonEta", &dileptonEta);
    outTree->Branch("dileptonPhi", &dileptonPhi);
    outTree->Branch("dileptonM", &dileptonM);
    outTree->Branch("dileptonDEta", &dileptonDEta);
    outTree->Branch("dileptonDPhi", &dileptonDPhi);
    outTree->Branch("dileptonDR", &dileptonDR);
    outTree->Branch("dileptonMKin", &dileptonMKin);

    // dilepton vertices
    //outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    //outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    //outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    //outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);
    
    // dijet
    outTree->Branch("dijetPt", &dijetPt);
    outTree->Branch("dijetEta", &dijetEta);
    outTree->Branch("dijetPhi", &dijetPhi);
    outTree->Branch("dijetM", &dijetM);
    outTree->Branch("dijetDEta", &dijetDEta);
    outTree->Branch("dijetDPhi", &dijetDPhi);
    outTree->Branch("dijetDR", &dijetDR);

    // jet, lepton
    outTree->Branch("l1j1DEta", &l1j1DEta);
    outTree->Branch("l1j1DPhi", &l1j1DPhi);
    outTree->Branch("l1j1DR", &l1j1DR);
    outTree->Branch("l1j2DEta", &l1j2DEta);
    outTree->Branch("l1j2DPhi", &l1j2DPhi);
    outTree->Branch("l1j2DR", &l1j2DR);
    outTree->Branch("l2j1DEta", &l2j1DEta);
    outTree->Branch("l2j1DPhi", &l2j1DPhi);
    outTree->Branch("l2j1DR", &l2j1DR);
    outTree->Branch("l2j2DEta", &l2j2DEta);
    outTree->Branch("l2j2DPhi", &l2j2DPhi);
    outTree->Branch("l2j2DR", &l2j2DR);

    // jet, photon
    outTree->Branch("j1PhotonDEta", &j1PhotonDEta);
    outTree->Branch("j1PhotonDPhi", &j1PhotonDPhi);
    outTree->Branch("j1PhotonDR", &j1PhotonDR);
    outTree->Branch("j2PhotonDEta", &j2PhotonDEta);
    outTree->Branch("j2PhotonDPhi", &j2PhotonDPhi);
    outTree->Branch("j2PhotonDR", &j2PhotonDR);
    outTree->Branch("jPhotonDRMax", &jPhotonDRMax);
    outTree->Branch("jPhotonDRMin", &jPhotonDRMin);

    // three body
    outTree->Branch("llgPt", &llgPt);
    outTree->Branch("llgEta", &llgEta);
    outTree->Branch("llgPhi", &llgPhi);
    outTree->Branch("llgM", &llgM);
    outTree->Branch("llgPtOverM", &llgPtOverM);
    outTree->Branch("llgMKin", &llgMKin);
    outTree->Branch("l1PhotonDEta", &l1PhotonDEta);
    outTree->Branch("l1PhotonDPhi", &l1PhotonDPhi);
    outTree->Branch("l1PhotonDR", &l1PhotonDR);
    outTree->Branch("l2PhotonDEta", &l2PhotonDEta);
    outTree->Branch("l2PhotonDPhi", &l2PhotonDPhi);
    outTree->Branch("l2PhotonDR", &l2PhotonDR);
    outTree->Branch("lPhotonDRMax", &lPhotonDRMax);
    outTree->Branch("lPhotonDRMin", &lPhotonDRMin);
    outTree->Branch("dileptonPhotonDEta", &dileptonPhotonDEta);
    outTree->Branch("dileptonPhotonDPhi", &dileptonPhotonDPhi);
    outTree->Branch("dileptonPhotonDR", &dileptonPhotonDR);
    outTree->Branch("ptt", &ptt);
    outTree->Branch("zgBigTheta", &zgBigTheta);
    outTree->Branch("zgLittleTheta", &zgLittleTheta);
    outTree->Branch("zgPhi", &zgPhi);
    outTree->Branch("genBigTheta", &genBigTheta);
    outTree->Branch("genLittleTheta", &genLittleTheta);
    outTree->Branch("genPhi", &genPhi);

    // other 
    outTree->Branch("llgJJDEta", &llgJJDEta);
    outTree->Branch("llgJJDPhi", &llgJJDPhi);
    outTree->Branch("llgJJDR", &llgJJDR);
    outTree->Branch("zepp", &zepp);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    ReportPostBegin();
}

Bool_t hzgAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);
    
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    
    genWeight = 1;
    if (!isData) {
        if (fGenEvtInfo->weight < 0) {
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
            hTotalEvents->Fill(maxBin);
        }
    }

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
   
    //bool sync_print = false;
    const bool noTrig = false;//true;
    bool sync_print_precut = false;

    if (sync_print_precut) {          
        //ULong64_t myEventNumber = 4368674045;
        if (!(
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 405) || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 419) || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 413) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 443) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 447) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 409) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 421) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 430) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 435) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 3 && fInfo->evtNum == 438) 
                )
            ) return kTRUE;

        //cout << "run, lumi, event" << endl;
        //cout << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << endl;
    }
          
    ///////////////////////
    // Generator objects // 
    ///////////////////////
    vector<TGenParticle*> genLeptons;
    vector<TGenParticle*> genPhotons;
    int ZtoQuarkID = 5;
	//cout << "/////////////////////////////////////////////////\n";
    fB = false;
    fPh = false;
    fLep = false;
    if (!isData) {
        unsigned count = 0;             
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {       
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
		if(particle->parent >0){	
               		TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
			if( (particle->pdgId ==  1 
			||  particle->pdgId ==  2
			||  particle->pdgId ==  5 
			||  particle->pdgId ==  3 
			||  particle->pdgId ==  4 
			||  particle->pdgId ==  5) 
			&&  mother->pdgId == 23) nB +=1;
			if(particle->pdgId == -leptonTagCharge*13 && fabs(mother->pdgId) == 24) nLep +=1;
			if(particle->pdgId == 22 && mother->pdgId == 25) nPh +=1;	

			if( (	particle->pdgId == 1	 
				||  particle->pdgId ==  2
				||  particle->pdgId ==  5 
				||  particle->pdgId ==  3 
				||  particle->pdgId ==  4 
				||  particle->pdgId ==  5) 
				&& mother->parent > 0){
				if (fabs(mother->pdgId) == 23){
					TGenParticle* momZ = (TGenParticle*) fGenParticleArr->At(mother->parent);
					for(int j =0; j < fGenParticleArr->GetEntries(); j++){
						TGenParticle* particle2 = (TGenParticle*) fGenParticleArr->At(j);
						if(particle2->parent > 0 && momZ->parent > 0){
							TGenParticle* mom2 = (TGenParticle*) fGenParticleArr->At(particle2->parent);
							if( mom2->parent > 0){

    								//nEl_elel =0; nEl_mumu =0; nEl_tautau =0; nEl_hh =0;
								//nMu_elel =0; nMu_mumu =0; nMu_tautau =0; nMu_hh =0;
    								//nTau_elel =0; nTau_mumu =0; nTau_tautau =0; nTau_hh =0;
			    					//nHad_elel =0; nHad_mumu =0; nHad_tautau =0; nHad_hh =0;
								if(
									particle2->parent == particle->parent
									//&& particle2->pdgId == -ZtoQuarkID
									&&  (particle->pdgId ==  1
						                        ||  particle->pdgId ==  2
						                        ||  particle->pdgId ==  5
						                        ||  particle->pdgId ==  3
						                        ||  particle->pdgId ==  4
						                        ||  particle->pdgId ==  5)
									){
										fB = true;
										//nB += 1;
										if(particle->pt <  particle2->pt){
											genJetOneP4.SetPtEtaPhiM(particle2->pt, particle2->eta, particle2->phi, particle2->mass);	
											genJetTwoP4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);	
										}
										else{	
											genJetOneP4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);	
											genJetTwoP4.SetPtEtaPhiM(particle2->pt, particle2->eta, particle2->phi, particle2->mass);
										}	
								}
								if( 
									momZ->pdgId == 25 
									&& particle2->parent == mother->parent
									){
										if(particle2->pdgId == 22){
											fPh = true;
											//nPh += 1;
											genPhotonP4.SetPtEtaPhiM(particle2->pt, particle2->eta, particle2->phi, particle2->mass);	
										}	
									}
								if(
									(particle2->pdgId == -leptonTagCharge*13)// || particle2->pdgId == 11)
									&& fabs(mom2->pdgId) == 24
									) {		
										fLep = true;
										//nLep += 1;
										genLeptonP4.SetPtEtaPhiM(particle2->pt, particle2->eta, particle2->phi, particle2->mass);	
								}
							}
						}
					}
				}
			} 
		}
	}
	flagGen = ( fB && fPh && fLep );
	if(flagGen){
		//if(genPhotonP4.Pt() > Sel_pt_ph && genLeptonP4.Pt() > Sel_pt_mu1){
    			outGenTree->Fill();
			nGen++;
		//}
	}

	nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

    }else {
       	nPartons = 0;
    }

    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    vector<string> passTriggerNames;
    if(!noTrig){
	    for (unsigned i = 0; i < triggerNames.size(); ++i) {
       		bool triggered = false;
	        triggered = trigger->pass(triggerNames[i], fInfo->triggerBits);
	        passTrigger |= triggered;

	        if (triggered) 
	        	passTriggerNames.push_back(triggerNames[i]);
	    }
    }
    else{
	passTrigger = true;
    }

    if (!passTrigger)// && isData)
        return kTRUE;
    else 
	nTrig += 1;
    hTotalEvents->Fill(3);

    /////////////////////
    /////////////////////
    // Fill event info //
    /////////////////////
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        nPU = fInfo->nPUmean;
        puWeight = weights->GetPUWeight(nPU); // pileup reweighting
        eventWeight *= puWeight;
    } else {
        nPU = 0;
    }    
    ///////////////////
    ///////////////////
    // Select objects//
    ///////////////////
    ///////////////////



    
    //////////////
    // Vertices //
    //////////////
    //TVertex* thePV;
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        //thePV = (TVertex *)fPVArr->At(0);
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        xPV = pv.X();
        yPV = pv.Y();
        zPV = pv.Z();
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /*if (sync_print_precut) {
        cout << "pvx, pvy, pvz, ndof" << endl;
        cout << thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << thePV->ndof << endl;
    }*/





    ///////////
    // MUONS //
    ///////////
    vector<TMuon*> muons;
    vector<TLorentzVector> veto_muons;

    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                    0, 0);
        }
        muon->pt = muonSF*muon->pt; 
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

	//cout << "mu\n";
	// Selection
	if( 
		//particleSelector->PassMuonID(muon, cuts->looseMuID)
		particleSelector->PassMuonID(muon, cuts->tightMuID)
		&& muon->pt 	   > PreSel_pt_mu
		&& fabs(muon->eta) < PreSel_eta_mu
                && GetMuonIsolation(muon)/muon->pt < 0.15
		){ 
		muons.push_back(muon);
		nPre_Mu++;
	}  
        // muons for jet veto
        if ( 
		particleSelector->PassMuonID(muon, cuts->vetoMuID) 
		&& muon->pt        > PreSel_pt_mu
		&& fabs(muon->eta) < PreSel_eta_mu
              	//&& GetMuonIsolation(muon)/muonP4.Pt() < 0.15
           ) {
            veto_muons.push_back(muonP4);
        }
	//cout << "mu passed\n";
    }
    //if(muons.size()==0 )
    //	noPass++; 

    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    ////////////////////////////////////////////////////////////////////////////



    ///////////////
    // ELECTRONS //
    ///////////////
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->calibPt, electron->eta, electron->phi, ELE_MASS);
        
        if ( 
		particleSelector->PassElectronID( electron, cuts->tightElID)
		&& electron->pt > PreSel_pt_el
		&& electron->eta < PreSel_eta_el
	) {
		electrons.push_back(electron);
		nPre_El++;
        }

	if (
		particleSelector->PassElectronID( electron, cuts->vetoElID)
	){
		veto_electrons.push_back(electronP4);
	}
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    ////////////////////////////////////////////////////////////////////////////






    //////////
    // TAUS //
    //////////
    vector<TTau*> taus;
    vector<TLorentzVector> veto_taus;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
        TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons and jets
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (tauP4.DeltaR(mu) < 0.3) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (tauP4.DeltaR(el) < 0.3) {
                elOverlap = true;
                break;
            }
        }

        // apply tau energy scale correction (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_energy_scale)
        if (!isData) {
            if (tau->decaymode == 0) {
                tau->pt *= 0.995;
            } else if (tau->decaymode == 1) {
                tau->pt *= 1.01;
            } else if (tau->decaymode == 10) {
                tau->pt *= 1.006;
            }
        }

        if ( 
                tau->pt > 18
                && abs(tau->eta) < 2.3
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);
    ////////////////////////////////////////////////////////////////////////////





    /////////////
    // PHOTONS //
    /////////////
    vector <TPhoton*> photons;
    vector<TLorentzVector> veto_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        if ( 
		//particleSelector->PassPhotonID(photon, cuts->preSelPhID)
		particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
		&& photon->pt        > PreSel_pt_ph
		&& fabs(photon->eta) < PreSel_eta_ph
	){
            photons.push_back(photon);
            veto_photons.push_back(photonP4);
	    nPre_Ph++;
        }
    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);  
    ////////////////////////////////////////////////////////////////////////////




    //////////
    // JETS //
    //////////
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        if (
		particleSelector->PassJetID(jet, cuts->looseJetID)
		//particleSelector->PassJetID(jet, cuts->tightJetID)
                && jet->pt        > PreSel_pt_jet
                && fabs(jet->eta) < PreSel_eta_jet
		//&& particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))
                //&& !muOverlap 
                //&& !elOverlap
                //&& !phoOverlap
           ) {
            
            jets.push_back(jet);
            nPre_Jet++;

            if (fabs(jet->eta) <= 2.4) { 
                hadronicP4 += jetP4;
                sumJetPt += jetP4.Pt();

                if (isData) {
                    if (jet->bmva > 0.9432) { 
                        ++nBJets;
                    } else {
                        ++nCentralJets;
                    }
                } else {
                    if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))) { 
                        ++nBJets;
			jet->bTagged = 1;
                    } else {
                        ++nCentralJets;
                    }
                }
            } else {
                if (fabs(jet->eta) > 2.5) {
                    hadronicP4 += jetP4;
                    sumJetPt += jetP4.Pt();
                    ++nFwdJets;
                }
            }
        }
    }

    //std::sort(jets.begin(), jets.end(), sort_by_btag);
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;
    metNC  = fInfo->pfMET;
    metPhiNC = fInfo->pfMETphi;

    /*if (sync_print_precut) {
        std::cout << "met, metPhi, metNC, metPhiNC" << std::endl;
        std::cout << met << ", " << metPhi << ", " << metNC << ", " << metPhiNC << std::endl;
    }*/

    /* HT */
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();
    /////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////
    //  Overlap Removal //
    //////////////////////
    if(muons.size() > 0 && photons.size() > 0 && jets.size() > 0){  	
    	//cout << muons.size() << " " << photons.size() << " " << jets.size() << endl;    	
	TLorentzVector tempJet, tempMuon,tempPhoton;

	for (unsigned int i = 0; i < muons.size(); ++i) {
		for (unsigned int j = 0; j < jets.size(); ++j) {
			tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
			tempJet.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
			
			if(tempJet.DeltaR(tempMuon) < 0.4)
				muons.erase(muons.begin()+i);
			
		}
	}

    		
	    for (unsigned int i = 0; i < photons.size(); ++i) {
		for (unsigned int j = 0; j < jets.size(); ++j) {
			tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0);
			tempJet.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
			if(tempJet.DeltaR(tempPhoton) < 0.4)
				photons.erase(photons.begin()+i);
    		}
    	}
		
    }

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nTaus      = taus.size();
    nPhotons   = photons.size();

    if(params->selection == "mugmjetjet"){
	if (jets.size() < 2 )
		return kTRUE;
	hTotalEvents->Fill(6);
	if (muons.size() < 1 && electrons.size() < 1)
            	return kTRUE;
        hTotalEvents->Fill(7);
        if (photons.size() < 1)
            	return kTRUE;
        hTotalEvents->Fill(8);


	////////////
	// LEPTON //
	////////////
	unsigned int electronOneIndex = 0;
	bool hasValidElectron = false;

        unsigned int muonOneIndex = 0;
	bool hasValidMuon = false;

	bool ElFlag = false;
	bool MuFlag = false;
	/*
	if(muons.size() == 0 && electrons.size() == 0 ) return kTRUE;
	if(muons.size() == 0 && electrons.size() > 1 ) ElFlag = 1;
	if(muons.size() > 1 && electrons.size() == 0 ) MuFlag = 1;
	if(muons.size() > 1 && electrons.size() > 1 ){
		if(muons[0]->q == leptonTagCharge)
			MuFlag = true;
		else{
			if( electrons[0]->q == leptonTagCharge)
				ElFlag = true;
			else return kTRUE;
		}
	}
	*/
	MuFlag = true;
	if(MuFlag){
		//////////
		// MUON //
		//////////
		for(unsigned int i = 0; i < muons.size(); i++){
			if(
				muons[i]->q == leptonTagCharge
				&& muons[i]->pt > Sel_pt_mu1
			){
			muonOneIndex = i;
			hasValidMuon = true;	
			nSel_Mu++;
			break;
			}
		}
		if(!hasValidMuon)
			return kTRUE;
	        hTotalEvents->Fill(9);
		leptonOneP4.SetPtEtaPhiM(muons[muonOneIndex]->pt, muons[muonOneIndex]->eta, muons[muonOneIndex]->phi, MUON_MASS);
	}else if(ElFlag){
		//////////////
		// ELECTRON //
		//////////////
		for(unsigned int i =0; i < electrons.size(); i++){
			if(
				electrons[i]->q == leptonTagCharge
				&& electrons[i]->pt > Sel_pt_el1
			){
				electronOneIndex = i;
				hasValidElectron = true;	
				nSel_El++;
				break;
			}
		}
	       	 if (!hasValidElectron)
       		    return kTRUE;
	        hTotalEvents->Fill(9);
		leptonOneP4.SetPtEtaPhiM(electrons[electronOneIndex]->pt, electrons[electronOneIndex]->eta, electrons[electronOneIndex]->phi, ELE_MASS);
	}
	////////////////////////////////////////////////////////////////////////	


	//////////
	// JETS //
	//////////
	
        // checking for dijet tag
        isDijetTag = false;
        isTightDijetTag = false;
        //unsigned int jetOneIndex = 0;
        //unsigned int jetTwoIndex = 0;
	TLorentzVector storedJetOne, storedJetTwo;
	if(jets.size() > 1){
                for (unsigned int i = 0; i < jets.size(); ++i) {
			for (unsigned int j = i+1; j < jets.size(); ++j) {
                        	TLorentzVector tempJetOne, tempJetTwo;
                        
	                        tempJetOne.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
	                        tempJetTwo.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);

	                        TLorentzVector tempDijet = tempJetOne + tempJetTwo;
	                        //float zeppen = llgP4.Eta() - (tempJetOne.Eta() + tempJetTwo.Eta())/2.;
	                        if ( 
					tempDijet.M() >= 50  && tempDijet.M() <= 120
	                                //&& tempJetOne.DeltaR(leptonOneP4) >= 0.4
					//&& tempJetTwo.DeltaR(leptonOneP4) >= 0.4
	                                //&& tempJetOne.DeltaR(photonOneP4) >= 0.4 
					//&& tempJetTwo.DeltaR(photonOneP4) >= 0.4
					//&&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
					//&&  fabs(zeppen) <= 2.5 
                            		//&&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
	                           ) {
					//if(jets[i]->bTagged || jets[j]->bTagged){
					//	jets[i]->bTagged = false;
					//	jets[j]->bTagged = false;

	                        		//storedDijet = tempJetOne + tempJetTwo;
			                        isDijetTag = true;
	        	                    	isTightDijetTag = true;
	                        	    	//jetOneIndex = i;
		                            	//jetTwoIndex = j;
						nJetPass_Sel +=1;
    						jetOneP4.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
    						jetTwoP4.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
 
						break;
						//nSel_Jet++;
					//}		
				}
		     	}	
		if(isDijetTag) break;
		}
	}
	//cout << "JET " << nJetPass_Sel << endl;
	//nJetPass_Sel = 0;
	if(!isDijetTag)
		return kTRUE;

    	//jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
    	//jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
        
	hTotalEvents->Fill(11);
	nSel_Jet++;


	////////////
	// Photon //
	////////////
        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            	TLorentzVector tempPhoton;
		TLorentzVector tempDiJet;
		TLorentzVector tempJJG;
            	//TLorentzVector tempDilepton;
	        //TLorentzVector tempLLG;
            	tempPhoton.SetPtEtaPhiM(photons[i]->calibPt, photons[i]->eta, photons[i]->phi, 0.);
            	tempDiJet = jetOneP4 + jetTwoP4;
		tempJJG = tempDiJet + tempPhoton;

            	//tempDilepton = leptonOneP4 + leptonTwoP4;
            	//tempLLG = leptonOneP4 + leptonTwoP4 + tempPhoton;
            	//tempLLG = dileptonP4 + tempPhoton;
            	//float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            	//float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);

            	float this_dr1 = jetOneP4.DeltaR(tempPhoton);
            	float this_dr2 = jetTwoP4.DeltaR(tempPhoton);
            	if (
                	tempPhoton.Pt() > Sel_pt_ph 
			//////// CHANGE UPCOMING LINES
                	//&& tempPhoton.Et()/tempJJG.M() > (15.0/110.0) 
                	// && tempDilepton.M() + tempLLG.M() > 185.0 
                	// && dileptonP4.M() + tempLLG.M() > 185.0 
                	&& tempJJG.M() > 100. && tempJJG.M() < 140. 
                	//&& this_dr1 >= 0.4 && this_dr2 >= 0.4
                	) {
                	hasValidPhoton = true;
                	photonIndex = i;
			nSel_Ph++;
			//cout << "Photon ind::"<< photonIndex << endl;
	    		//cout << "Pt:: " << photons[photonIndex]->calibPt << endl;
                	break;
            		}
        }

        if (!hasValidPhoton)
        	return kTRUE;
        hTotalEvents->Fill(12);

        /*if (sync_print_precut)
            cout << "passed the photon selection" << endl;*/

        //TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->calibPt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
       	if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(13);

        //TLorentzVector llgP4 = dileptonP4 + photonOneP4;

        // DY photon overlap removal
        vetoDY = false;
        for (unsigned int i = 0; i < genPhotons.size(); ++i) {
            TGenParticle *pho = genPhotons.at(i); 
            if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
                TLorentzVector thisGenPhotonP4;
                thisGenPhotonP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
                if (thisGenPhotonP4.DeltaR(photonOneP4) < 0.1) {
                    vetoDY = true;
                    break;
                }
            }
        }
	////////////////////////////////////////////////////////////////////////////////
	
        if (isDijetTag) {
            //jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
            //jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
            jetOnePt = jetOneP4.Pt();
            jetOneEta = jetOneP4.Eta();
            jetOnePhi = jetOneP4.Phi();
            jetOneM   = jetOneP4.M();
            
            jetTwoPt = jetTwoP4.Pt();
            jetTwoEta = jetTwoP4.Eta();
            jetTwoPhi = jetTwoP4.Phi();
            jetTwoM   = jetTwoP4.M();

            TLorentzVector dijet = jetOneP4 + jetTwoP4;
            dijetPt = dijet.Pt();
            dijetEta = dijet.Eta();
            dijetPhi = dijet.Phi();
            dijetM = dijet.M();
            dijetDEta = fabs(jetOneP4.Eta() - jetTwoP4.Eta());
            dijetDPhi = fabs(jetOneP4.DeltaPhi(jetTwoP4));
            dijetDR = jetOneP4.DeltaR(jetTwoP4);

            l1j1DEta = fabs(leptonOneP4.Eta() - jetOneP4.Eta());
            l1j1DPhi = fabs(leptonOneP4.DeltaPhi(jetOneP4));
            l1j1DR = leptonOneP4.DeltaR(jetOneP4);
            l1j2DEta = fabs(leptonOneP4.Eta() - jetTwoP4.Eta());
            l1j2DPhi = fabs(leptonOneP4.DeltaPhi(jetTwoP4));
            l1j2DR = leptonOneP4.DeltaR(jetTwoP4);
            l2j1DEta = fabs(leptonTwoP4.Eta() - jetOneP4.Eta());
            l2j1DPhi = fabs(leptonTwoP4.DeltaPhi(jetOneP4));
            l2j1DR = leptonTwoP4.DeltaR(jetOneP4);
            l2j2DEta = fabs(leptonTwoP4.Eta() - jetTwoP4.Eta());
            l2j2DPhi = fabs(leptonTwoP4.DeltaPhi(jetTwoP4));
            l2j2DR = leptonTwoP4.DeltaR(jetTwoP4);

            j1PhotonDEta = fabs(jetOneP4.Eta() - photonOneP4.Eta());
            j1PhotonDPhi = fabs(jetOneP4.DeltaPhi(photonOneP4));
            j1PhotonDR = jetOneP4.DeltaR(photonOneP4);
            j2PhotonDEta = fabs(jetTwoP4.Eta() - photonOneP4.Eta());
            j2PhotonDPhi = fabs(jetTwoP4.DeltaPhi(photonOneP4));
            j2PhotonDR = jetTwoP4.DeltaR(photonOneP4);
        
            if (j1PhotonDR > j2PhotonDR) {
                jPhotonDRMax = j1PhotonDR;
                jPhotonDRMin = j2PhotonDR; 
            }
            else {
                jPhotonDRMax = j2PhotonDR;
                jPhotonDRMin = j1PhotonDR;
            }
            
            //zepp = llgP4.Eta() - (jetOneP4.Eta() + jetTwoP4.Eta())/2.;
            //llgJJDEta = fabs(llgP4.Eta() - dijet.Eta());
            //llgJJDPhi = fabs(llgP4.DeltaPhi(dijet));
            //llgJJDR = llgP4.DeltaR(dijet);

        }
 



    //}
   } // end mugm selection

    ///////////////////
    // Fill jet info //
    ///////////////////
    

    if (!isDijetTag) {
        if (jets.size() > 0) {
            jetOnePt = jets[0]->pt;
            jetOneEta = jets[0]->eta;
            jetOnePhi = jets[0]->phi;
            jetOneM = jets[0]->mass;
            jetOneTag    = jets[0]->csv;
        } else {
            jetOnePt = 0.;
            jetOneEta = 0.;
            jetOnePhi = 0.;
            jetOneM = 0.;
            jetOneTag    = 0.;
        }

        if (jets.size() > 1) {
            jetTwoPt = jets[1]->pt;
            jetTwoEta = jets[1]->eta;
            jetTwoPhi = jets[1]->phi;
            jetTwoM = jets[1]->mass;
            jetTwoTag    = jets[1]->csv;
        } else {
            jetTwoPt = 0.;
            jetTwoEta = 0.;
            jetTwoPhi = 0.;
            jetTwoM = 0.;
            jetTwoTag    = 0.;
        }
    }


    if (!isData && genLeptons.size() == 2) {
        genLeptonOneId = genLeptons[0]->pdgId;
        genLeptonOnePt = genLeptons[0]->pt;
        genLeptonOneEta = genLeptons[0]->eta;
        genLeptonOnePhi = genLeptons[0]->phi;
        genLeptonTwoId = genLeptons[1]->pdgId;
        genLeptonTwoPt = genLeptons[1]->pt;
        genLeptonTwoEta = genLeptons[1]->eta;
        genLeptonTwoPhi = genLeptons[1]->phi;  
    }
    else {
        genLeptonOnePt = 0.;
        genLeptonOneEta = 0.;
        genLeptonOnePhi = 0.;
        genLeptonTwoPt = 0.;
        genLeptonTwoEta = 0.;
        genLeptonTwoPhi = 0.;
    }

    TLorentzVector genPhotonP4;
    if (!isData && genPhotons.size() > 0) {
        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photonOnePt, photonOneEta, photonOnePhi, 0.);
        float min_phot_dr = 1000.;
        for (unsigned int i = 0; i < genPhotons.size(); i++) {
            TLorentzVector tmpGenPhot;
            tmpGenPhot.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
            float this_dr = tmpGenPhot.DeltaR(photonOneP4);
            if (this_dr < min_phot_dr) {
                genPhotonP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                genPhotonFHPFS = genPhotons[i]->fromHardProcessFinalState;
                genPhotonIPFS = genPhotons[i]->isPromptFinalState;
                min_phot_dr = this_dr;
            }
        }
        genPhotonPt = genPhotonP4.Pt();
        genPhotonEta = genPhotonP4.Eta();
        genPhotonPhi = genPhotonP4.Phi();
    }
    else {
        genPhotonPt = 0.;
        genPhotonEta = 0.;
        genPhotonPhi = 0.;
    }
        
    // gen angles
    if (!isData && genLeptons.size() == 2 && genPhotons.size() > 0) {        
        TLorentzVector l_minus, l_plus; 
        if (genLeptonOneId > 0) {
            l_minus.SetPtEtaPhiM(genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi, genLeptons[0]->mass);
            l_plus.SetPtEtaPhiM(genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi, genLeptons[1]->mass);
        }
        else {
            l_minus.SetPtEtaPhiM(genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi, genLeptons[1]->mass);
            l_plus.SetPtEtaPhiM(genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi, genLeptons[0]->mass);
        }
     
        TLorentzVector dileptonP4Gen, llgP4Gen;
        dileptonP4Gen = l_minus + l_plus;
        llgP4Gen = dileptonP4Gen + genPhotonP4;

        TVector3 llgFrame = -1*llgP4Gen.BoostVector();
        dileptonP4Gen.Boost(llgFrame);
        l_minus.Boost(llgFrame);
        l_minus.Boost(dileptonP4Gen.BoostVector());
        genLittleTheta = cos(dileptonP4Gen.Angle(l_minus.Vect()));
        genBigTheta = cos(dileptonP4Gen.Angle(llgP4Gen.Vect()));
        
        TVector3 ppAxis(0, 0, 1);
        TVector3 zAxis = dileptonP4Gen.Vect().Unit();
        TVector3 yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rotation;
        rotation = rotation.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4Gen.Transform(rotation);
        l_minus.Transform(rotation);
        genPhi = l_minus.Phi();
    }
    else {
        genLittleTheta = -9.;
        genBigTheta = -9.;
        genPhi = -9.;
    }

    outTree->Fill();
    this->passedEvents++;

    /*if (sync_print_precut) {
        cout << "event should have been filled" << endl;
        }*/
	/*
    if(nGen == 0){
	genJetOneP4.SetPtEtaPhiM(0, 0, 0, 0);	
	genJetTwoP4.SetPtEtaPhiM(0, 0, 0, 0);
	genPhotonP4.SetPtEtaPhiM(0, 0, 0, 0);
	genLeptonP4.SetPtEtaPhiM(0, 0, 0, 0);	
    	outGenTree->Fill();

	jetOneP4.SetPtEtaPhiM(0, 0, 0, 0);	
	jetTwoP4.SetPtEtaPhiM(0, 0, 0, 0);
	photonP4.SetPtEtaPhiM(0, 0, 0, 0);
	leptonOneP4.SetPtEtaPhiM(0, 0, 0, 0);	
    	outTree->Fill();
    }
	*/
    return kTRUE;
}

void hzgAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void hzgAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void hzgAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
    std::cout << "NTRIG :: " << nTrig << " NO PASS " << noPass << endl << endl;
    std::cout << " =============================================================" << std::endl;
    std::cout << " W\\Z     ||  elel   ||  mumu   ||   tt   ||    h h   ||\n";
    std::cout << " ---------------------------------------------------------\n"; 
    std::cout << "   e     ||  " << nEl_elel << "   ||   " << nEl_mumu << "   ||   " << nEl_tautau << "   ||   " << nEl_hh << "   ||  " << endl;
    std::cout << "  mu     ||  " << nMu_elel << "   ||   " << nMu_mumu << "   ||   " << nMu_tautau << "   ||   " << nMu_hh << "   ||  " << endl;
    std::cout << " tau     ||  " << nTau_elel << "   ||   " << nTau_mumu << "   ||   " << nTau_tautau << "   ||   " << nTau_hh << "   ||  " << endl;
    std::cout << "   h     ||  " << nHad_elel << "   ||   " << nHad_mumu << "   ||   " << nHad_tautau << "   ||   " << nHad_hh << "   ||  " << endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "nGen " << nGen << " |||| nB " << nB << " || nPh " << nPh << " || nLep " << nLep << endl;

    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  Pre |||| Mu:: " << nPre_Mu << " || El:: " << nPre_El << " || Jet:: " << nPre_Jet << " || Ph:: " << nPre_Ph << endl;
    std::cout << "  Sel |||| Mu:: " << nSel_Mu << " || El:: " << nSel_El << " || Jet:: " << nSel_Jet << " || Ph:: " << nSel_Ph<< endl;
    std::cout << "  ===========================================================" << std::endl;
	cout <<  " JETS IN LOOP " << nJetPass_Sel << endl;
    std::cout << "  ===========================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<hzgAnalyzer> selector(new hzgAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float hzgAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float hzgAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = el->chHadIso + std::max(0., (double)el->neuHadIso + el->gammaIso - rho*effArea[iEta]);

    return combIso;
}

float hzgAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(pho->scEta) > etaBins[i] && fabs(pho->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = pho->chHadIso + std::max(0., (double)pho->neuHadIso + pho->gammaIso - rho*effArea[iEta]);

    return combIso;
}

void hzgAnalyzer::EvalMuonEnergyResolution(std::map<string, float> mva_input_floats, std::map<string, int> mva_input_ints, float &mean, float &sigma, float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;
     // one-time MVA initialization: get trainings
    if (!mvaMean) {
        // current working directory
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");
         // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;
         RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_muons"));
        if (!ws) FATAL("TFile::Get() failed");
         invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");
         mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");
     } // one-time initialization
     // load necessary tree branches
    float rho                      = mva_input_floats["rho"];
    float muEnergy                 = mva_input_floats["muEnergy"];
    float muEta                    = mva_input_floats["muEta"];
    float muChi2NDF                = mva_input_floats["muTkChi2"];
    Int_t muNumberOfValidTrkLayers = mva_input_ints["muNumberOfValidTrkLayers"];
    Int_t muNumberOfValidPixelHits = mva_input_ints["muNumberOfValidPixelHits"];
    Int_t muNumberOfValidMuonHits  = mva_input_ints["muNumberOfValidMuonHits"];
    Int_t muStations               = mva_input_ints["muStations"];
    float muPFIsoR04_CH            = mva_input_floats["muPFIsoR04_CH"];
    float muPFIsoR04_NH            = mva_input_floats["muPFIsoR04_NH"];
    float muPFIsoR04_Pho           = mva_input_floats["muPFIsoR04_Pho"];
    float muPFIsoR04_PU            = mva_input_floats["muPFIsoR04_PU"];
     // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = muEnergy;
    *invar[2] = muEta;
    *invar[3] = muChi2NDF;
    *invar[4] = muNumberOfValidTrkLayers;
    *invar[5] = muNumberOfValidPixelHits;
    *invar[6] = muNumberOfValidMuonHits;
    *invar[7] = muStations;
    *invar[8] = muPFIsoR04_CH;
    *invar[9] = muPFIsoR04_NH;
    *invar[10] = muPFIsoR04_Pho;
    *invar[11] = muPFIsoR04_PU;
     mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();
     // NOTE: negative = infinite; powers were fixed at the training level
    powerL = -1;
    powerR = -1;
 }
 void hzgAnalyzer::EvalElectronEnergyResolution(std::map<string, float> mva_inputs, float &mean, float &sigma, 
                                                    float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;
    static RooAbsReal* mvaPowerL;
    static RooAbsReal* mvaPowerR;
     // one-time MVA initialization: get trainings
    if (!mvaMean) {
        // current working directory
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");
         // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;
         RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_electrons"));
        if (!ws) FATAL("TFile::Get() failed");
         invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");
        invar[12] = ws->var("var13");
        invar[13] = ws->var("var14");
        invar[14] = ws->var("var15");
         mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");
        mvaPowerL = ws->function("limPowerL");
        mvaPowerR = ws->function("limPowerR");
     } // one-time initialization
     // load necessary tree branches
    float rho = mva_inputs["rho"];
    float elEnergy = mva_inputs["elEnergy"];
    float elScEta = mva_inputs["elScEta"];
    float elScPhi = mva_inputs["elScPhi"];
    float elR9 = mva_inputs["elR9"];
    float elE1x5OverE = mva_inputs["elE1x5OverE"];
    float elE2x5OverE = mva_inputs["elE2x5OverE"];
    float elE5x5OverE = mva_inputs["elE5x5OverE"];
    float elFBrem = mva_inputs["elFBrem"];
    float elHOverE = mva_inputs["elHOverE"];
    float elSigmaIEtaIEta = mva_inputs["elSigmaIEtaIEta"];
    float elPFIsoR04_CH = mva_inputs["elPFIsoR04_CH"];
    float elPFIsoR04_NH = mva_inputs["elPFIsoR04_NH"];
    float elPFIsoR04_Pho = mva_inputs["elPFIsoR04_Pho"];
    float elPFIsoR04_PU = mva_inputs["elPFIsoR04_PU"];
     // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = elEnergy;
    *invar[2] = elScEta;
    *invar[3] = elScPhi;
    *invar[4] = elR9;
    *invar[5] = elE1x5OverE;
    *invar[6] = elE2x5OverE;
    *invar[7] = elE5x5OverE;
    *invar[8] = elFBrem;
    *invar[9] = elHOverE;
    *invar[10] = elSigmaIEtaIEta;
    *invar[11] = elPFIsoR04_CH;
    *invar[12] = elPFIsoR04_NH;
    *invar[13] = elPFIsoR04_Pho;
    *invar[14] = elPFIsoR04_PU;
     mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();
    powerL = mvaPowerL->getVal();
    powerR = mvaPowerR->getVal();
}

 void hzgAnalyzer::find_optimized(double* p, double &e1, double& e2)
{
   // Returns best-fitted (most probable) energies.
    static TF2* fun = NULL;
    // one-time initialization
   if (!fun) {
      fun = new TF2("fun", NegativeProbability, 1, 2000, 1, 2000, 16);
      fun->SetNpx(50);
      fun->SetNpy(50);
   }
    for (int i = 0; i < 16; i++)
      fun->SetParameter(i, p[i]);
    // limit the search range by +-3sigma regions
   fun->SetRange(p[2] * (1 - 3 * p[5]), p[3] * (1 - 3 * p[11]),
                 p[2] * (1 + 3 * p[5]), p[3] * (1 + 3 * p[11]));
    fun->GetMinimumXY(e1, e2);
}
