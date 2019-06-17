#include "zgAnalyzer_MuEff.hh"
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

zgAnalyzer_MuEff::zgAnalyzer_MuEff() : BLTSelector()
{

}

zgAnalyzer_MuEff::~zgAnalyzer_MuEff()
{

}

void zgAnalyzer_MuEff::Begin(TTree *tree)
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

    if (params->selection == "elgmjetjet"){
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    }
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
    }
    if (params->selection == "mumugm"){
	/// SINGLE MUON TRIGGERS
	triggerNames.push_back("HLT_IsoTkMu24_v*");
	triggerNames.push_back("HLT_IsoMu24_v*");
	
	/// DOUBLE MUON TRIGGERS
	triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
	triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
	triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
 	triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    }
    if (params->selection == "elelgm"){
	/// SINGLE ELECTRON TRIGGERS
	//triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
	
	/// DOUBLE ELECTRON TRIGGERS
	triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }

    ///////////////////////////////////
    ////  Some Selection variables  ///
    ///////////////////////////////////
    
    // Debugging
    debugFlag = false;
    //debugFlag = true;

    // Set the cuts
    cuts.reset(new Cuts()); 
    particleSelector.reset(new ParticleSelector(*params, *cuts));
    
    if(params->selection=="mugmjetjet"){
	// Pre/Selection
	if(triggerNames[0] == "HLT_IsoMu24_v*"){
		PreSel_pt_mu        = 23;
		Sel_pt_mu1          = 25;
		PreSel_pt_ph        = 13; 
		//Sel_pt_ph           = 20;//15;
		Sel_pt_ph = 15;
	}
	else {
		PreSel_pt_mu 	= 16;
		Sel_pt_mu1      = 18;
		PreSel_pt_ph    = 21;
		Sel_pt_ph       = 23;
	}
    }
    if(params->selection=="mumugm"){
	// Pre/Selection
	PreSel_pt_mu        = 10;
	Sel_pt_mu1          = 18;
    	Sel_pt_mu2          = 9;

	PreSel_pt_ph        = 10; 
	Sel_pt_ph 	    = 10;
    }
    if(params->selection=="elelgm"){
	// Pre/Selection
	PreSel_pt_el        = 15;
	Sel_pt_el1          = 24;
	Sel_pt_el2          = 13;

	PreSel_pt_ph        = 10; 
	Sel_pt_ph = 10;
    }
    
    nMuProbe = 0;
    nMuProbePass = 0;

    PreSel_eta_mu       = 2.4;
    PreSel_eta_el       = 2.4;
    PreSel_eta_ph       = 2.4;

    PreSel_pt_jet  	= 30;
    PreSel_eta_jet 	= 4.7;

    leptonTagCharge = -1;

    nPre_Mu = 0;
    nPre_El = 0;
    nPre_Jet = 0;
    nPre_Ph = 0;
    
    nPartCount = 0;

    nSel_Mu = 0;
    nSel_El = 0;
    nSel_Jet = 0;
    nSel_Ph = 0;

    // Gen level counts 
    nGen = 0;
    nB =0;
    nPh = 0;
    nLep = 0;
    nLeptonMom = 0;
    LeptonMom = 0; 

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
    string outBugTreeName = params->get_output_treename("tree");
    outGenTreeName = "Gen" + outGenTreeName;
    outBugTreeName = "Bug" + outBugTreeName;

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");
    outGenTree = new TTree(outGenTreeName.c_str(), "bltTree");
    outBugTree = new TTree(outBugTreeName.c_str(), "bltTree");

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

    // Efficientcy
    outTree->Branch("nMuProbe"    , &nMuProbe);
    outTree->Branch("nMuProbePass", &nMuProbePass);

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
    outTree->Branch("jetThreeP4", &jetThreeP4);

    outTree->Branch("jetOnePt", &jetOnePt);
    outTree->Branch("jetOneEta", &jetOneEta);
    outTree->Branch("jetOnePhi", &jetOnePhi);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoPt", &jetTwoPt);
    outTree->Branch("jetTwoEta", &jetTwoEta);
    outTree->Branch("jetTwoPhi", &jetTwoPhi);
    outTree->Branch("jetTwoTag", &jetTwoTag);


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
    

    // Debug
    outBugTree->Branch("MuonMult",&nLepton);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    ReportPostBegin();
}

Bool_t zgAnalyzer_MuEff::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    eventStep = 1;
    this->totalEvents++;
    hTotalEvents->Fill(eventStep);
    
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    
    genWeight = 1;
    if (!isData) {
        if (fGenEvtInfo->weight < 0) {
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
	    //cout << "MaxBin::  " << maxBin << endl;
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
    
    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    eventStep++;
    hTotalEvents->Fill(eventStep);

    /* Trigger selection */
    bool passTrigger = false;
    vector<string> passTriggerNames;
    if(!noTrig){
	    for (unsigned i = 0; i < triggerNames.size(); ++i) {
       		bool triggered = false;
	        triggered = trigger->pass(triggerNames[i], fInfo->triggerBits);
	        passTrigger |= triggered;

	        if (triggered) {
			//cout << "Triggered\n";
	        	passTriggerNames.push_back(triggerNames[i]);
		}
	    }
    }
    else{
	passTrigger = true;
    }

    if (!passTrigger)// && isData)
        return kTRUE;
    else 
	nTrig += 1;
    eventStep++;
    hTotalEvents->Fill(eventStep);

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
    eventStep++;
    hTotalEvents->Fill(eventStep);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /*if (sync_print_precut) {
        cout << "pvx, pvy, pvz, ndof" << endl;
        cout << thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << thePV->ndof << endl;
    }*/

    if(debugFlag)
	cout << "--- Vertices" << endl;   



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
	        && trigger->passObj("HLT_IsoMu24_v*",1,muon->hltMatchBits)
                //&& GetMuonIsolation(muon)/muon->pt < 0.10
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
   
    
    vector<TMuon*> muonsProbe;

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

	//cout << "mu\n";
	// Selection
	//if( 
	//	//particleSelector->PassMuonID(muon, cuts->looseMuID)
	//	particleSelector->PassMuonID(muon, cuts->tightMuID)
	//	&& muon->pt 	   > PreSel_pt_mu
	//	&& fabs(muon->eta) < PreSel_eta_mu
        //        //&& GetMuonIsolation(muon)/muon->pt < 0.15
        //        && GetMuonIsolation(muon)/muon->pt < 0.10
	//	){ 
	//	muons.push_back(muon);
	//	nPre_Mu++;
	//}
	if(
		muon->pt 	   > PreSel_pt_mu
		&& fabs(muon->eta) < PreSel_eta_mu
	){
		muonsProbe.push_back(muon); 
	}
    }
    sort(muonsProbe.begin(), muonsProbe.end(), sort_by_higher_pt<TMuon>);
    ////////////////////////////////////////////////////////////////////////////
    if(debugFlag)
	cout << "--- Muons" << endl;   
    


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
		&& particleSelector->PassElectronIso( electron, cuts->tightElIso) 
		//cuts->tightElIso is not really used all the values are hardcoded in the function itself
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

    if(debugFlag)
    	cout << "-------- Electrons\n";

    /////////////
    // PHOTONS //
    /////////////
    float EAPho[7][3] = {{0.0360,0.0597,0.1210},
			 {0.0377,0.0807,0.1107},
			 {0.0306,0.0629,0.0699},
			 {0.0283,0.0197,0.1056},
			 {0.0254,0.0184,0.1457},
			 {0.0217,0.0284,0.1719},
			 {0.0167,0.0591,0.1998}
			};
    vector <TPhoton*> photons;
    vector<TLorentzVector> veto_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        if ( 
		//particleSelector->PassPhotonID(photon, cuts->preSelPhID)
		//particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
		particleSelector->PassPhotonID(photon, cuts->mediumPhID)
		//&& particleSelector->PassPhotonIso(photon, cuts->mediumPhIso,EAPho)
		//&& GetPhotonIsolation(photon,cuts->mediumPhIso,fInfo->rhoJet)
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

    if(debugFlag)
    	cout << "-------- Photons\n";



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
			//jet->bTagged = 1;
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

    if(debugFlag)
    	cout << "-------- Jets\n";
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
			
			//if(tempJet.DeltaR(tempMuon) < 0.4 && muons.size() != 0)
			//	muons.erase(muons.begin()+i);
			//cout << "Muon " << j << endl;
			if(tempJet.DeltaR(tempMuon) < 0.4 && jets.size() != 0)
				jets.erase(jets.begin()+j);

			
		}
	}

	if(jets.size() != 0){
    		for (unsigned int i = 0; i < photons.size(); ++i) {
			for (unsigned int j = 0; j < jets.size(); ++j) {
				tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0);
				tempJet.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
				//if(tempJet.DeltaR(tempPhoton) < 0.4 && photons.size() != 0)
					//photons.erase(photons.begin()+i);

				//cout << "Photon " << i << " - " << j << endl;
				if(tempJet.DeltaR(tempPhoton) < 0.4 && jets.size() != 0)
					jets.erase(jets.begin()+j);
    			}
    		}
	}
		
    }
    
    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nPhotons   = photons.size();
    if(debugFlag)
	cout << "--------start\n";
    if(params->selection == "mumugm"){
	if(muons.size() == 0)
		return kTRUE;
	////////////
	// LEPTON //
	////////////
	unsigned int electronOneIndex = 0;
	unsigned int electronTwoIndex = 0;
	bool hasValidElectron = false;

        unsigned int muonOneIndex = 0;
        unsigned int muonTwoIndex = 0;
	bool hasValidMuon = false;
	//bool hasValidPair = false;

	bool ElFlag = false;
	bool MuFlag = false;
	
	nMuProbe = 0;
	nMuProbePass = 0;	

	// Efficientcy 	//
	//int isGlobalOne  = (muons[0]->typeBits & baconhep::kGlobal) > 0;
	//int isGlobalTwo  = (muons[1]->typeBits & baconhep::kGlobal) > 0;
	//int isTrackerOne = (muons[0]->typeBits & baconhep::kTraker) > 0;
	//int isTrackerTwo = (muons[1]->typeBits & baconhep::kTraker) > 0;
	
	if(debugFlag)
		cout << "-------- Trigger Matching and count\n";
	for(unsigned int i = 0;i < muonsProbe.size(); i++){
		TLorentzVector tempLepOneP4,tempLepTwoP4, tempDILEP;
		tempLepOneP4.SetPtEtaPhiM(muons[muonOneIndex]->pt,muons[muonOneIndex]->eta, muons[muonOneIndex]->phi,MUON_MASS);	
		tempLepTwoP4.SetPtEtaPhiM(muonsProbe[i]->pt,muonsProbe[i]->eta, muonsProbe[i]->phi,MUON_MASS);
		tempDILEP = tempLepOneP4 + tempLepTwoP4;
		
		if(
			muonsProbe[i]->q == -muons[0]->q
               		&& GetMuonIsolation(muonsProbe[i])/muonsProbe[i]->pt < 0.10
			&& tempDILEP.M() < 105 && tempDILEP.M() > 80
		){	
			nMuProbe = 1;
			hasValidMuon =true;		
			leptonOneP4.SetPtEtaPhiM(muonsProbe[i]->pt, muonsProbe[i]->eta, muonsProbe[i]->phi, MUON_MASS);
	        	if( trigger->passObj("HLT_IsoMu24_v*",1,muonsProbe[i]->hltMatchBits)){
				nMuProbePass = 1;
				leptonTwoP4.SetPtEtaPhiM(muonsProbe[i]->pt, muonsProbe[i]->eta, muonsProbe[i]->phi, MUON_MASS);
				//hasValidPair = true;
				break;
			}
			else 
				leptonTwoP4.SetPtEtaPhiM(0,0,0, MUON_MASS);
		}
			
	}	
	if(!hasValidMuon)
		return kTRUE;
	//if(!hasValidPair)
	//	return kTRUE;
	

	
	/*	
	////////////
	// Photon //
	////////////
        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            	TLorentzVector tempPhoton;
		//TLorentzVector tempDiJet;
		//TLorentzVector tempJJG;
            	//TLorentzVector tempDilepton;
	        //TLorentzVector tempLLG;
            	
		//tempPhoton.SetPtEtaPhiM(photons[i]->calibPt, photons[i]->eta, photons[i]->phi, 0.);
		tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0.);
            	//cout << "Photon PT ::" << photons[i]->pt << endl;
            	//cout << "Photon calibPT ::" << photons[i]->calibPt << endl;
		//tempDiJet = jetOneP4 + jetTwoP4;
		//tempJJG = tempDiJet + tempPhoton;

            	//tempDilepton = leptonOneP4 + leptonTwoP4;
            	//tempLLG = leptonOneP4 + leptonTwoP4 + tempPhoton;
            	//tempLLG = dileptonP4 + tempPhoton;
            	//float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            	//float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);

            	//float this_dr1 = jetOneP4.DeltaR(tempPhoton);
            	//float this_dr2 = jetTwoP4.DeltaR(tempPhoton);
            	if (
                	tempPhoton.Pt() > Sel_pt_ph 
			//////// CHANGE UPCOMING LINES
                	//&& tempPhoton.Et()/tempJJG.M() > (15.0/110.0) 
                	// && tempDilepton.M() + tempLLG.M() > 185.0 
                	// && dileptonP4.M() + tempLLG.M() > 185.0 
                	//&& tempJJG.M() > 100. && tempJJG.M() < 140. 
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

	if(debugFlag)
		cout << "-------- Photons selection\n";
        if (!hasValidPhoton)
		cout << "";
        //	return kTRUE;
        
	if(debugFlag)
		cout << "-------- Photons Before Fill\n";
        	
	eventStep++;
        hTotalEvents->Fill(eventStep);

	if(debugFlag)
		cout << "-------- Photons Save\n";
	cout    << photonIndex << endl;
	cout	<< photons[photonIndex]->pt  << endl;
	cout	<< photons[photonIndex]->eta << endl;
	cout	<< photons[photonIndex]->phi << endl;
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->pt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);


	if(debugFlag)
		cout << "-------- Photons End Selection\n";
	*/
    //}
   } // end mujetjetgm selection


   if(debugFlag)
   	cout << "-------- Fill\n";

    outTree->Fill();
    this->passedEvents++;

    return kTRUE;
}

void zgAnalyzer_MuEff::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zgAnalyzer_MuEff::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zgAnalyzer_MuEff::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
    std::cout << "  W\\Z     ||  elel   ||  mumu   ||   tt   ||    h h   ||\n";
    std::cout << "  ---------------------------------------------------------\n"; 
    std::cout << "    e    ||  " << nEl_elel << "   ||   " << nEl_mumu << "   ||   " << nEl_tautau << "   ||   " << nEl_hh << "   ||  " << endl;
    std::cout << "   mu    ||  " << nMu_elel << "   ||   " << nMu_mumu << "   ||   " << nMu_tautau << "   ||   " << nMu_hh << "   ||  " << endl;
    std::cout << "  tau    ||  " << nTau_elel << "   ||   " << nTau_mumu << "   ||   " << nTau_tautau << "   ||   " << nTau_hh << "   ||  " << endl;
    std::cout << "    h    ||  " << nHad_elel << "   ||   " << nHad_mumu << "   ||   " << nHad_tautau << "   ||   " << nHad_hh << "   ||  " << endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  nGen " << nGen << " |||| nB " << nB << " || nPh " << nPh << " || nLep " << nLep << endl;
    std::cout << "  ============================================================" << std::endl;
    std::cout << "  ============================================================" << std::endl;
    std::cout << "  NTRIG :: " << nTrig << " NO PASS :: " << noPass << " Lepton :: " << nLeptonMom << endl;
    std::cout << "  nPartCount :: " << nPartCount << endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  Pre |||| Mu:: " << nPre_Mu << " || El:: " << nPre_El << " || Ph:: " << nPre_Ph << " || Jet:: " << nPre_Jet << endl;
    std::cout << "  Sel |||| Mu:: " << nSel_Mu << " || El:: " << nSel_El << " || Ph:: " << nSel_Ph << " || Jet:: " << nSel_Jet <<  endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  Probe     :: " << nMuProbe     << std::endl;
    std::cout << "  ProbePass :: " << nMuProbePass << std::endl;
    std::cout << "  ===========================================================" << std::endl;
    std::cout << "  ===========================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<zgAnalyzer_MuEff> selector(new zgAnalyzer_MuEff());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zgAnalyzer_MuEff::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float zgAnalyzer_MuEff::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float zgAnalyzer_MuEff::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    ///float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    float effArea[8] = {0.1210,0.1107,0.0699,0.1056,0.1457,0.1719,0.1998}; //Spring2016
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(pho->scEta) > etaBins[i] && fabs(pho->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }
    float combIso = pho->chHadIso + std::max(0., (double)pho->neuHadIso + pho->gammaIso - rho*effArea[iEta]);

    return combIso;
}

void zgAnalyzer_MuEff::EvalMuonEnergyResolution(std::map<string, float> mva_input_floats, std::map<string, int> mva_input_ints, float &mean, float &sigma, float &alphaL, float &powerL, float &alphaR, float &powerR) 
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
 void zgAnalyzer_MuEff::EvalElectronEnergyResolution(std::map<string, float> mva_inputs, float &mean, float &sigma, 
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

 void zgAnalyzer_MuEff::find_optimized(double* p, double &e1, double& e2)
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
