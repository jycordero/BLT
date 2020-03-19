#include "zgAnalyzer.h"
#include <map>
#include <fstream>
#include <math.h>

#include <TSystem.h>
#include <TF2.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
//#include "KinZfitter/KinZfitter/interface/KinZfitter.h"

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

zgAnalyzer::zgAnalyzer() : BLTSelector()
{

}

zgAnalyzer::~zgAnalyzer()
{

}

void zgAnalyzer::Begin(TTree *tree)
{
    bool confDebug = true;
    //bool confDebug = false;
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
            std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    string PeriodFolder; 

    if(confDebug){
	    cout << options.at(0) << endl
		 << options.at(1) << endl
		 << options.at(2) << endl
		 << options.at(3) << endl
    		 << options.at(4) << endl;
    }


 // Param Formating
    if(params->period == "2016")
	params->period = "2016Legacy";
    else if (params->period == "2017")
	params->period = "2017ReReco";

    if(params->period == "2016Legacy")
	PeriodFolder= "Legacy2016";
    else if (params->period == "2017ReReco")
	PeriodFolder = "ReReco2017";
    else if (params->period == "2016ReReco")
	PeriodFolder = "ReReco2016";

    
    if(confDebug)
	    cout << "--- Cuts and Particle Selector \n";

    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    if(confDebug)
    	cout << "--- Trigger \n";

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));
    if(params->period == "2016" || params->period ==  "2016ReReco" || params->period == "2016Legacy"){
	    if (params->selection == "mumu" || params->selection == "mumug") {
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
	    }
	    else if (params->selection == "ee" || params->selection == "elelg") {
		triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
		triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
	    }
	    else if (params->selection == "tautaug") { // select one muon plus one hadronic tau (for now)
		triggerNames.push_back("HLT_IsoMu24_v*");
		triggerNames.push_back("HLT_IsoTkMu24_v*");
	    }
    }
    else if(params->period == "2017" || params->period ==  "2017ReReco" || params->period == "2017Legacy"){
	if (params->selection == "mumu" || params->selection == "mumug") {
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
		cout << "\n\n\nTHE RIGHT TRIGGERS\n\n\n";
	}
	else if (params->selection == "ee" || params->selection == "elelg") {
		triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
		triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
	}
    }
 
    if(confDebug)
	cout << "--- Weights \n";
    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask


    if(confDebug)
    	cout << "--- LumiMask \n";
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if( params->period == "2016Legacy"){
        string lumiFile = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/lumiMask/" + lumiFile;

        muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/ReReco2016/roccor.Run2.v3/RoccoR2016.txt");
    }
    else if(params->period == "2016ReReco"){
        string lumiFile = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/lumiMask/" + lumiFile;

        muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/roccor.Run2.v3/RoccoR2017.txt");
    }
    else if(params->period == "2017ReReco"){
        string lumiFile = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/lumiMask/" + lumiFile;

        muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/roccor.2017.v0/RoccoR2017v0.txt");
    }
    else if(params->period == "20178ReReco"){
        string lumiFile = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/lumiMask/" + lumiFile;

        muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/roccor.2018.v3/RoccoR2018.txt");
    }


    lumiMask.AddJSONFile(jsonFileName);

    if(confDebug)
    	cout << "--- Declarations \n";
    // muon momentum corrections
    //muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/ReReco2016/rcdata.2016.v3");
    rng = new TRandom3();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    Sgen = 0;
    SgenAccep = 0;

    // event data
    outTree->Branch("runNumber"    , &runNumber);
    outTree->Branch("evtNumber"    , &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection"  , &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("nPV"          , &nPV);
    outTree->Branch("nPU"          , &nPU);
    outTree->Branch("nPartons"     , &nPartons);
    outTree->Branch("xPV"          , &xPV);
    outTree->Branch("yPV"          , &yPV);
    outTree->Branch("zPV"          , &zPV);
    outTree->Branch("Rho"          , &Rho);

    outTree->Branch("met"     , &met);
    outTree->Branch("metPhi"  , &metPhi);
    outTree->Branch("metNC"   , &metNC);
    outTree->Branch("metPhiNC", &metPhiNC);
    outTree->Branch("ht"      , &ht);
    outTree->Branch("htPhi"   , &htPhi);
    outTree->Branch("htSum"   , &htSum);

    // weights
    outTree->Branch("genWeight"          , &genWeight);
    outTree->Branch("eventWeight"        , &eventWeight);
    outTree->Branch("puWeight"           , &puWeight);
    outTree->Branch("triggerWeight"      , &triggerWeight);
    outTree->Branch("elIDWeightOne"      , &elIDWeightOne);
    outTree->Branch("elIDWeightTwo"      , &elIDWeightTwo);
    outTree->Branch("elTrigWeightOne"    , &elTrigWeightOne);
    outTree->Branch("elTrigWeightTwo"    , &elTrigWeightTwo);
    outTree->Branch("muonIDWeightOne"    , &muonIDWeightOne);
    outTree->Branch("muonIDWeightTwo"    , &muonIDWeightTwo);
    outTree->Branch("muonISOWeightOne"   , &muonISOWeightOne);
    outTree->Branch("muonISOWeightTwo"   , &muonISOWeightTwo);
    outTree->Branch("muonTrigWeightOne"  , &muonTrigWeightOne);
    outTree->Branch("muonTrigWeightTwo"  , &muonTrigWeightTwo);
    outTree->Branch("photonIDWeight"     , &photonIDWeight);
    outTree->Branch("photonIsConvWeight" , &photonIsConvWeight);
    
    outTree->Branch("Sgen", &Sgen);
    outTree->Branch("SgenAccep", &SgenAccep);

    // leptons
    outTree->Branch("leptonOnePt"        , &leptonOnePt);
    outTree->Branch("leptonOneEta"       , &leptonOneEta);
    outTree->Branch("leptonOnePhi"       , &leptonOnePhi);
    outTree->Branch("leptonOnePtKin"     , &leptonOnePtKin);
    outTree->Branch("leptonOnePtKinJames", &leptonOnePtKinJames);
    outTree->Branch("leptonOneIso"       , &leptonOneIso);
    outTree->Branch("leptonOneFlavor"    , &leptonOneFlavor);
    outTree->Branch("leptonOneMother"    , &leptonOneMother);
    outTree->Branch("leptonOneD0"        , &leptonOneD0);
    outTree->Branch("leptonOneDZ"        , &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
    outTree->Branch("leptonOneECALDriven", &leptonOneECALDriven);
    outTree->Branch("leptonOneCharge"    , &leptonOneCharge);
    outTree->Branch("leptonOneTag"       , &leptonOneTag);

    outTree->Branch("leptonTwoPt"        , &leptonTwoPt);
    outTree->Branch("leptonTwoEta"       , &leptonTwoEta);
    outTree->Branch("leptonTwoPhi"       , &leptonTwoPhi);
    outTree->Branch("leptonTwoPtKin"     , &leptonTwoPtKin);
    outTree->Branch("leptonTwoPtKinJames", &leptonTwoPtKinJames);
    outTree->Branch("leptonTwoIso"       , &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor"    , &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother"    , &leptonTwoMother);
    outTree->Branch("leptonTwoD0"        , &leptonTwoD0);
    outTree->Branch("leptonTwoDZ"        , &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
    outTree->Branch("leptonTwoECALDriven", &leptonTwoECALDriven);
    outTree->Branch("leptonTwoCharge"    , &leptonTwoCharge);
    outTree->Branch("leptonTwoTag"       , &leptonTwoTag);
 
    outTree->Branch("tauDecayMode", &tauDecayMode);
    outTree->Branch("tauMVA"      , &tauMVA);

    outTree->Branch("isLeptonTag"    , &isLeptonTag);
    outTree->Branch("isDijetTag"     , &isDijetTag);
    outTree->Branch("isTightDijetTag", &isTightDijetTag);

    // photons
    outTree->Branch("photonOnePt"     , &photonOnePt);
    outTree->Branch("photonOneEta"    , &photonOneEta);
    outTree->Branch("photonOnePhi"    , &photonOnePhi);
    outTree->Branch("photonOneR9"     , &photonOneR9);
    outTree->Branch("photonOneMVA"    , &photonOneMVA);
    outTree->Branch("photonOneERes"   , &photonOneERes);
    outTree->Branch("passElectronVeto", &passElectronVeto);

    outTree->Branch("photonOneSieie" , &photonOneSieie); 
    outTree->Branch("photonOneHoverE", &photonOneHoverE);
    outTree->Branch("photonOneIneu"  , &photonOneIneu);
    outTree->Branch("photonOneIph"   , &photonOneIph);
    outTree->Branch("photonOneIch"   , &photonOneIch);

    outTree->Branch("photonOneSieip"     , &photonOneSieip);      
    outTree->Branch("photonOneSipip"     , &photonOneSipip);     
    outTree->Branch("photonOneSrr"       , &photonOneSrr);       
    outTree->Branch("photonOneE2x2"      , &photonOneE2x2);      
    outTree->Branch("photonOneE5x5"      , &photonOneE5x5);      
    outTree->Branch("photonOneScEtaWidth", &photonOneScEtaWidth); 
    outTree->Branch("photonOneScPhiWidth", &photonOneScPhiWidth); 
    outTree->Branch("photonOneScRawE"    , &photonOneScRawE); 
    outTree->Branch("photonOnePreShowerE", &photonOnePreShowerE); 
    outTree->Branch("photonOneScBrem"    , &photonOneScBrem); 

    // jets
    outTree->Branch("jetOnePt" , &jetOnePt);
    outTree->Branch("jetOneEta", &jetOneEta);
    outTree->Branch("jetOnePhi", &jetOnePhi);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoPt" , &jetTwoPt);
    outTree->Branch("jetTwoEta", &jetTwoEta);
    outTree->Branch("jetTwoPhi", &jetTwoPhi);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    // gen level objects 
    outTree->Branch("genLeptonOnePt" , &genLeptonOnePt);
    outTree->Branch("genLeptonOneEta", &genLeptonOneEta);
    outTree->Branch("genLeptonOnePhi", &genLeptonOnePhi);
    outTree->Branch("genLeptonOneId" , &genLeptonOneId);
    outTree->Branch("genLeptonTwoPt" , &genLeptonTwoPt);
    outTree->Branch("genLeptonTwoEta", &genLeptonTwoEta);
    outTree->Branch("genLeptonTwoPhi", &genLeptonTwoPhi);
    outTree->Branch("genLeptonTwoId" , &genLeptonTwoId);
    outTree->Branch("genPhotonPt"    , &genPhotonPt);
    outTree->Branch("genPhotonEta"   , &genPhotonEta);
    outTree->Branch("genPhotonPhi"   , &genPhotonPhi);
    outTree->Branch("genPhotonFHPFS" , &genPhotonFHPFS);
    outTree->Branch("genPhotonIPFS"  , &genPhotonIPFS);
    outTree->Branch("vetoDY"         , &vetoDY);
    outTree->Branch("genIsoPass"     , &genIsoPass);

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
    outTree->Branch("dileptonPt"       , &dileptonPt);
    outTree->Branch("dileptonEta"      , &dileptonEta);
    outTree->Branch("dileptonPhi"      , &dileptonPhi);
    outTree->Branch("dileptonM"        , &dileptonM);
    outTree->Branch("dileptonDEta"     , &dileptonDEta);
    outTree->Branch("dileptonDPhi"     , &dileptonDPhi);
    outTree->Branch("dileptonDR"       , &dileptonDR);
    outTree->Branch("dileptonMKin"     , &dileptonMKin);
    outTree->Branch("dileptonMKinJames", &dileptonMKinJames);

    // dilepton vertices
    //outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    //outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    //outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    //outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);
    
    // dijet
    outTree->Branch("dijetPt"  , &dijetPt);
    outTree->Branch("dijetEta" , &dijetEta);
    outTree->Branch("dijetPhi" , &dijetPhi);
    outTree->Branch("dijetM"   , &dijetM);
    outTree->Branch("dijetDEta", &dijetDEta);
    outTree->Branch("dijetDPhi", &dijetDPhi);
    outTree->Branch("dijetDR"  , &dijetDR);

    // jet, lepton
    outTree->Branch("l1j1DEta", &l1j1DEta);
    outTree->Branch("l1j1DPhi", &l1j1DPhi);
    outTree->Branch("l1j1DR"  , &l1j1DR);
    outTree->Branch("l1j2DEta", &l1j2DEta);
    outTree->Branch("l1j2DPhi", &l1j2DPhi);
    outTree->Branch("l1j2DR"  , &l1j2DR);
    outTree->Branch("l2j1DEta", &l2j1DEta);
    outTree->Branch("l2j1DPhi", &l2j1DPhi);
    outTree->Branch("l2j1DR"  , &l2j1DR);
    outTree->Branch("l2j2DEta", &l2j2DEta);
    outTree->Branch("l2j2DPhi", &l2j2DPhi);
    outTree->Branch("l2j2DR"  , &l2j2DR);

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
    outTree->Branch("llgPt"       , &llgPt);
    outTree->Branch("llgEta"      , &llgEta);
    outTree->Branch("llgPhi"      , &llgPhi);
    outTree->Branch("llgM"        , &llgM);
    outTree->Branch("llgPtOverM"  , &llgPtOverM);
    outTree->Branch("llgMKin"     , &llgMKin);
    outTree->Branch("llgMKinJames", &llgMKinJames);
    outTree->Branch("l1PhotonDEta", &l1PhotonDEta);
    outTree->Branch("l1PhotonDPhi", &l1PhotonDPhi);
    outTree->Branch("l1PhotonDR"  , &l1PhotonDR);
    outTree->Branch("l1PhotonPt"  , &l1PhotonPt);
    outTree->Branch("l1PhotonM"   , &l1PhotonM);
    outTree->Branch("l2PhotonDEta", &l2PhotonDEta);
    outTree->Branch("l2PhotonDPhi", &l2PhotonDPhi);
    outTree->Branch("l2PhotonDR"  , &l2PhotonDR);
    outTree->Branch("l2PhotonPt"  , &l2PhotonPt);
    outTree->Branch("l2PhotonM"   , &l2PhotonM);

    outTree->Branch("lPhotonDRMax"      , &lPhotonDRMax);
    outTree->Branch("lPhotonDRMin"      , &lPhotonDRMin);
    outTree->Branch("dileptonPhotonDEta", &dileptonPhotonDEta);
    outTree->Branch("dileptonPhotonDPhi", &dileptonPhotonDPhi);
    outTree->Branch("dileptonPhotonDR"  , &dileptonPhotonDR);
    outTree->Branch("ptt", &ptt);
    outTree->Branch("zgBigTheta", &zgBigTheta);
    outTree->Branch("zgLittleTheta", &zgLittleTheta);
    outTree->Branch("zgLittleThetaMY", &zgLittleThetaMY);
    outTree->Branch("zgPhi", &zgPhi);
    outTree->Branch("zgBigThetaJames", &zgBigThetaJames);
    outTree->Branch("zgLittleThetaJames", &zgLittleThetaJames);
    outTree->Branch("zgPhiJames"    , &zgPhiJames);
    outTree->Branch("genBigTheta"   , &genBigTheta);
    outTree->Branch("genLittleTheta", &genLittleTheta);
    outTree->Branch("genPhi"        , &genPhi);

    // other 
    outTree->Branch("llgJJDEta", &llgJJDEta);
    outTree->Branch("llgJJDPhi", &llgJJDPhi);
    outTree->Branch("llgJJDR"  , &llgJJDR);
    outTree->Branch("zepp"     , &zepp);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    string outHistNameGen = params->get_output_treename("TotalEventsGen");
    hTotalEventsGen = new TH1D(outHistNameGen.c_str(),"TotalEventsGen",5,0.5,5.5);

    ReportPostBegin();
}

Bool_t zgAnalyzer::Process(Long64_t entry)
{
    //bool debug = true;
    bool debug = false;
    //bool debugSelect = true;
    bool debugSelect = false;
    
    //bool debugSF = true;
    if(debug){
	cout << "=========================";
	cout << "Start ..........\n";
    }

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);
	
    if(debug)   
	    cout << "GetEntry.....\n";
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData); 
    if (debug)
	cout << "SetRealData.....\n";



    genWeight = 1;
    if (!isData) {
        
	if(debug)
		cout << "... If Gen...\n";
	if (fGenEvtInfo->weight < 0) {
	    if(debug)
		cout << "... .... If negWeight...\n";
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
            hTotalEvents->Fill(maxBin);

	    if(debug)
		cout << "... .... Fill hist...\n";
        }
    }
    
    if(debug)
	cout << "Weights....\n";
    //cout << "---------RunNum:: " << fInfo->runNum << endl; 
    if (entry%10000==0)  
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
   
    //bool sync_print = false;
    bool sync_print_precut = false;
    //bool sync_print_precut = true;

    if (sync_print_precut) {          
       
        // for checking non-overlapping events
        if (!(
                // muon channel
                (fInfo->runNum == 1 && fInfo->lumiSec == 40     && fInfo->evtNum == 7817)   || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 65     && fInfo->evtNum == 12829)  || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 75     && fInfo->evtNum == 14933)  ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 57     && fInfo->evtNum == 11340)  ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1550   && fInfo->evtNum == 309928) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1569   && fInfo->evtNum == 313636) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1569   && fInfo->evtNum == 313680) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1570   && fInfo->evtNum == 313826) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 428    && fInfo->evtNum == 85455)  ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 421    && fInfo->evtNum == 84164)  ||
                //electron channel
                (fInfo->runNum == 1 && fInfo->lumiSec == 63    && fInfo->evtNum == 12497) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 79    && fInfo->evtNum == 15658) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 775   && fInfo->evtNum == 154801) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1098  && fInfo->evtNum == 219407) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1549  && fInfo->evtNum == 309788) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1543  && fInfo->evtNum == 308557) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 1555  && fInfo->evtNum == 310804) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 379   && fInfo->evtNum == 75615) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 380   && fInfo->evtNum == 75912) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 408   && fInfo->evtNum == 81429) 
                )
            ) 
		return kTRUE;
	   

        cout << "run, lumi, event" << endl;
        cout << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << endl;
    }
          
    ///////////////////////
    // Generator objects // 
    ///////////////////////

    vector<TGenParticle*> genLeptons;
    vector<TGenParticle*> genPhotons;
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
            //cout << "\t" << particle->pt << ", " << particle->eta;
            //cout << endl;

            if (
                    particle->status == 23 
                    && (fabs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            if ((fabs(particle->pdgId) == 11 || fabs(particle->pdgId) == 13) and particle->parent > 0) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                if (fabs(mother->pdgId) == 23) 
                    genLeptons.push_back(particle);
            }
                
            // saving photons     
            if (fabs(particle->pdgId) == 22) 
                genPhotons.push_back(particle);    

        }
        ///////////////////////////
        // veto flag ar gen level 
	bool genPhotonFlag = false;
        for (unsigned int i = 0; i < genPhotons.size(); ++i) {
                TGenParticle *pho = genPhotons.at(i);
                if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
			if(pho->pdgId == 22){
				genPhotonFlag = true;
			}
                }
        }
	if(genPhotonFlag)
		hTotalEvents->Fill(29);

        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

    } else {
        nPartons = 0;
    }

    if(!isData){
	std::cout << "Enter Acceptance loop" << endl;
	if(genPhotons.size() > 0 && genLeptons.size() == 2){
		hTotalEventsGen->Fill(1);
		Sgen++;

		TLorentzVector genPhotonP4, genLeptonOneP4,genLeptonTwoP4;

		if(params->selection == "mumug"){

			genPhotonP4   .SetPtEtaPhiM(genPhotons[0]->pt, genPhotons[0]->eta, genPhotons[0]->phi, 0.);
			genLeptonOneP4.SetPtEtaPhiM(genLeptons[0]->pt, genLeptons[0]->eta, genLeptons[0]->phi, MUON_MASS);
			genLeptonTwoP4.SetPtEtaPhiM(genLeptons[1]->pt, genLeptons[1]->eta, genLeptons[1]->phi, MUON_MASS);

			if(genPhotons[0]->pt > 15 && genLeptons[0]->pt > 25 && genLeptons[1]->pt > 20){
				if((fabs(genPhotons[0]->eta) < 1.4442 || fabs(genPhotons[0]->eta) > 1.566) 
				&& fabs(genPhotons[0]->eta) < 2.4 
				&& fabs(genLeptons[0]->eta) < 2.4 
				&& fabs(genLeptons[1]->eta) < 2.4
				&& genPhotonP4.DeltaR(genLeptonOneP4) > 0.7
				&& genPhotonP4.DeltaR(genLeptonTwoP4) > 0.7
				&& GetGenIsolation(genPhotons[0]) < 5
				)

				hTotalEventsGen->Fill(4);
				hTotalEvents->Fill(15);
				SgenAccep++;

				genLeptonOnePt  = genLeptons[0]->pt;
				genLeptonOneEta = genLeptons[0]->eta;
				genLeptonOnePhi = genLeptons[0]->phi;

				genLeptonTwoPt  = genLeptons[1]->pt;
				genLeptonTwoEta = genLeptons[1]->eta;
				genLeptonTwoPhi = genLeptons[1]->phi;

				genPhotonPt  = genPhotons[0]->pt;
				genPhotonEta = genPhotons[0]->eta;
				genPhotonPhi = genPhotons[0]->phi;
			}
		}
		else if(params->selection == "elelg"){
			genPhotonP4   .SetPtEtaPhiM(genPhotons[0]->pt, genPhotons[0]->eta, genPhotons[0]->phi, 0.);
			genLeptonOneP4.SetPtEtaPhiM(genLeptons[0]->pt, genLeptons[0]->eta, genLeptons[0]->phi, ELE_MASS);
			genLeptonTwoP4.SetPtEtaPhiM(genLeptons[1]->pt, genLeptons[1]->eta, genLeptons[1]->phi, ELE_MASS);

			if(genPhotons[0]->pt > 15 && genLeptons[0]->pt > 25 && genLeptons[1]->pt > 20){
				if((fabs(genPhotons[0]->eta) < 1.4442 || fabs(genPhotons[0]->eta) > 1.566)
 				&& (fabs(genLeptons[0]->eta) < 1.4442 || fabs(genLeptons[0]->eta) > 1.566) 
 				&& (fabs(genLeptons[1]->eta) < 1.4442 || fabs(genLeptons[1]->eta) > 1.566) 
				&& fabs(genPhotons[0]->eta) < 2.5 
				&& fabs(genLeptons[0]->eta) < 2.5
				&& fabs(genLeptons[1]->eta) < 2.5
				&& genPhotonP4.DeltaR(genLeptonOneP4) > 0.7
				&& genPhotonP4.DeltaR(genLeptonTwoP4) > 0.7
				&& GetGenIsolation(genPhotons[0]) < 5
				)
				hTotalEventsGen->Fill(15);
				SgenAccep++;

				genLeptonOnePt  = genLeptons[0]->pt;
				genLeptonOneEta = genLeptons[0]->eta;
				genLeptonOnePhi = genLeptons[0]->phi;

				genLeptonTwoPt  = genLeptons[1]->pt;
				genLeptonTwoEta = genLeptons[1]->eta;
				genLeptonTwoPhi = genLeptons[1]->phi;

				genPhotonPt  = genPhotons[0]->pt;
				genPhotonEta = genPhotons[0]->eta;
				genPhotonPhi = genPhotons[0]->phi;

			}
		}
	}

	

    }

    if (debug){
	cout<< "Before Lumi mask passed\n";
	cout << "RunNum: " << fInfo->runNum << "  --  " << fInfo->lumiSec << endl;
    }
    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    if (debug)
	cout<< "Lumi mask passed\n";

    /* Trigger selection */
    bool passTrigger = false;
    vector<string> passTriggerNames;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        bool triggered = false;
        triggered = trigger->pass(triggerNames[i], fInfo->triggerBits);
        passTrigger |= triggered;
        if (triggered) { 
            passTriggerNames.push_back(triggerNames[i]);
        }
    }

    if (!passTrigger)// && isData)
        return kTRUE;
    hTotalEvents->Fill(3);

    if(debug)
	cout << "HLT Trigger passed\n" << endl;
    /////////////////////
    // Fill event info //
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
    // Select objects//
    ///////////////////

    /* Vertices */
    TVertex* thePV;
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        thePV = (TVertex *)fPVArr->At(0);
        xPV = pv.X();
        yPV = pv.Y();
        zPV = pv.Z();
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);

    if(debug)
	cout << "PV passed\n" << endl;
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);
    Rho = fInfo->rhoJet;

    if (sync_print_precut) {
        cout << "pvx, pvy, pvz, ndof" << endl;
        cout << thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << thePV->ndof << endl;
    }

    if(debug)
	cout << "Start Object Filtering\n";
    /* MUONS */
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
	if(debug)
		cout << "Muons correction\n";
        if (sync_print_precut) {
            cout << "pt_before_roccor, pt_after_roccor" << endl;
            cout << muon->pt << ", ";
        }
	
        muon->pt = muonSF*muon->pt; 
        if (sync_print_precut)
            cout << muon->pt << endl;
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

        // Remove muons with very small deltaR // THIS IS NOT YET USED
        float minDeltaR = 1e6;
        for (unsigned j=0; j < muons.size(); ++j) {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1051);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR) {
                minDeltaR = dr;
            }
        }

        if (sync_print_precut) {
           cout << "pt , eta, phi, isGlobal, isTracker, nMatchStn, nValidHits, d0, dz, sip3d, nPixHits, nTkLayers, pt_track, ptErr, iso, trkIso" << endl;
           cout << muonP4.Pt() << "," << muonP4.Eta() << "," << 
                   muonP4.Phi() << "," <<
                   (muon->typeBits & baconhep::kGlobal) << "," << 
                   (muon->typeBits & baconhep::kTracker) << "," <<  
                   muon->nMatchStn << "," << muon->nValidHits << ", " << muon->d0 << "," << 
                   muon->dz << "," << muon->sip3d << "," << 
                   muon->nPixHits << "," << muon->nTkLayers << "," << 
                   muon->pt << "," << muon->ptErr << "," << 
                   GetMuonIsolation(muon) << ", " << muon->trkIso << endl;
        }

        if (   
               // tight muon ID and ISO
               muonP4.Pt() > 10.
               && fabs(muonP4.Eta()) < 2.4
               && (muon->typeBits & baconhep::kPFMuon) 
               && (muon->typeBits & baconhep::kGlobal) 
               && muon->muNchi2    < 10.
               && muon->nMatchStn  > 1
               && muon->nPixHits   > 0
               && fabs(muon->d0)   < 0.2
               && fabs(muon->dz)   < 0.5
               && muon->nTkLayers  > 5 
               && muon->nValidHits > 0
               && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
           ){
	    if(muons.size() >= 2)
		    cout << "Muons Selected\n";
            muons.push_back(muon);
        } 
                    

        // muons for jet veto
        if (
                muonP4.Pt() > 10
                // tight muon ID and ISO
                && fabs(muonP4.Eta()) < 2.4
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
           ) {
            veto_muons.push_back(muonP4);
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    if(debug)
	cout << "----Muon collection";

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        //assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->calibPt, electron->eta, electron->phi, ELE_MASS);
        
        if(sync_print_precut || debugSelect) {
           cout << "electron pt, calibpt, eta, sc_eta, phi, d0, dz, sip3d, iso, mva, pass_mva, pass_id, pass_iso" << endl;
           //cout << electronP4.Pt() << "," << electronP4.Eta() << "," << 
           cout << electron->pt << "," << electron->calibPt << "," <<  electron->eta << ", " << 
                   electron->scEta << ", " << electron->phi << "," << electron->d0 << "," << 
                   electron->dz << "," << electron->sip3d << "," << 
                   GetElectronIsolation(electron, fInfo->rhoJet) << ", " << electron->mvaFall17V2Iso << ", " << 
                   particleSelector->PassElectronMVA(electron, cuts->hzzMVAID) << 
		   particleSelector->PassElectronID(electron, cuts->tightElID) << 
		   particleSelector->PassElectronIso(electron,cuts->tightElIso) << 
		   endl;
        } 
	//cout << "     Electron pt:: " << electron->pt << " calib:: " << electron->calibPt << endl;
	//if(particleSelector->PassElectronID(electron, cuts->tightElID))
	//	cout << "Electron ID Pass\n";
	//if(particleSelector->PassElectronIso(electron,cuts->tightElIso))
	//	cout << "Electron ISO Pass\n";
        if (
                electron->calibPt > 11
                && fabs(electron->scEta) < 2.5
                //&& particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)
                && particleSelector->PassElectronID(electron, cuts->tightElID)
		&& particleSelector->PassElectronIso(electron,cuts->tightElIso)
                //&& GetElectronIsolation(electron, fInfo->rhoJet)/electronP4.Pt() < 0.35
                //&& fabs(electron->d0) < 0.5
                //&& fabs(electron->dz) < 1.0
                //&& fabs(electron->sip3d) < 4.0 
           ) {
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    if(electrons.size() > 0)
	cout << "Electron Selected with " << electrons.size() << endl;
    if(debug)
	cout << "----Electron collection" << endl;
    /* TAUS */
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

    if(debug)
	cout << "----Tau collection" << endl;

    /* PHOTONS */
    vector <TPhoton*> photons;
    vector<TLorentzVector> veto_photons;
    float EAPho[7][3] = {
			{0.0112, 	0.0668, 	0.1113},
			{0.0108, 	0.1054, 	0.0953},
			{0.0106, 	0.0786, 	0.0619},
			{0.01002, 	0.0233, 	0.0837},
			{0.0098, 	0.0078, 	0.1070},
			{0.0089, 	0.0028, 	0.1212},
			{0.0087, 	0.0137, 	0.1466} 
			};

    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photon->calibPt, photon->eta, photon->phi, 0.);
    	//cout << " ---------PhotonP4 PT :: " << photonP4.Pt() << endl;
    	if(debug){
		cout << " ---------Photon PT   :: " << photon->pt << " cal: " << photon->calibPt << endl;
		cout << "----------Photon Eta  :: " << photon->eta << " scEta: " << photon->scEta << endl;
		//cout << "----------Photon ID   :: " << particleSelector->PassPhotonID(photon, cuts->loosePhID) << endl;
		cout << "----------Photon ID   :: " << particleSelector->PassPhotonID(photon, cuts->preSelPhID) << endl;
		cout << "----------Photon Veto :: " << photon->passElectronVeto << endl;
	}

        if (sync_print_precut) {
            cout << "photon_pt, photon_calibpt, photon_eta, photon_sc_eta, photon_phi, photon_mva, pass_electron_veto" << endl;
            cout << photon->pt << ", " << photon->calibPt << ", " << photon->eta << ", " << photon->scEta << ", " << photon->phi << ", " << photon->mvaFall17V2 
                 << ", " << photon->passElectronVeto << endl;
        }

        if (
                // ID conditions
                photon->calibPt > 10
                && fabs(photon->scEta) < 2.5 
                && (fabs(photon->scEta) <= 1.4442 || fabs(photon->scEta) >= 1.566)
                //&& particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
                
                && particleSelector->PassPhotonID (photon, cuts->preSelPhID)
                && particleSelector->PassPhotonIso(photon, cuts->preSelPhIso,EAPho)	
                //
                //&& particleSelector->PassPhotonID (photon, cuts->mediumPhID)
                //&& particleSelector->PassPhotonIso(photon, cuts->mediumPhIso,EAPho)	
		&& GetWorstChIsolation(photon) < 15
                && photon->passElectronVeto
            ) {
            photons.push_back(photon);
            veto_photons.push_back(photonP4);
        }
    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);



    if(debug)
	cout << "----Photon collection" << endl;
    /* JETS */
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

        if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
        }

        // Prevent overlap of muons and jets
        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        if (sync_print_precut) {
            double jec = 0.;
            if (isData) { // fix for broken bacon JEC
                jec = particleSelector->JetCorrector(jet, "NONE");
                //jet->pt = jet->ptRaw*jec;
            }
            std::cout << "raw jet pt, jet pt, nate_jet_pt, jet eta, jet phi" << std::endl;
            std::cout << jet->ptRaw << ", " << jetP4.Pt() << ", " << jet->ptRaw*jec << ", " << jetP4.Eta() << ", " << jetP4.Phi() << std::endl;
        } 
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (jetP4.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (jetP4.DeltaR(el) < 0.4) {
                elOverlap = true;
                break;
            }
        }
        bool phoOverlap = false;
        for (const auto& pho: veto_photons) {
            if (jetP4.DeltaR(pho) < 0.4) {
                phoOverlap = true;
                break;
            }
        }

        if (
                jet->pt > 30 
                //&& fabs(jet->eta) < 4.7
                && fabs(jet->eta) < 2.5
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !phoOverlap
           ) {
            
            jets.push_back(jet);
            ++nJets;

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

    if(debug)
	cout << "----Jet collection" << endl;
    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;
    metNC  = fInfo->pfMET;
    metPhiNC = fInfo->pfMETphi;

    if (sync_print_precut) {
        std::cout << "met, metPhi, metNC, metPhiNC" << std::endl;
        std::cout << met << ", " << metPhi << ", " << metNC << ", " << metPhiNC << std::endl;
    }

    /* HT */
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nTaus      = taus.size();
    nPhotons   = photons.size();

    if(debug)
	cout << "---- Start Selection" << endl;
    if (params->selection == "mumug") {
        if (muons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);

        if(debugSelect){
		cout << "---------------------------" << params->selection << endl;
		cout << "--------------------------- MU:: " << nMuons << endl;
		cout << "--------------------------- EL:: " << nElectrons << endl;
		cout << "--------------------------- PH:: " << nPhotons << endl;
        }
	if(debugSelect)
		cout << "------------ Multiplicity cut" << endl;
 
        TLorentzVector leptonOneP4, leptonTwoP4;
        unsigned int muonOneIndex = 0;
        unsigned int muonTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            for (unsigned int j = i+1; j < muons.size(); ++j) {
                TLorentzVector tempMuonOne, tempMuonTwo;
                tempMuonOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
                tempMuonTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
                float thisMass = (tempMuonOne + tempMuonTwo).M();
                //if (
                //        //muons[i]->q != muons[j]->q
                //        //muons[i]->q == muons[j]->q // fake selection
                //        /*&&*/ muons[i]->pt > 20.0 
                //        && muons[j]->pt > 10.0
                //        && thisMass > 50.0
                //   ) {
                if (thisMass > 50.0) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempMuonOne;
                            leptonTwoP4 = tempMuonTwo;
                            muonOneIndex = i;
                            muonTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempMuonOne;
                        leptonTwoP4 = tempMuonTwo;
                        muonOneIndex = i;
                        muonTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }


        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        //if (leptonOneP4.Pt() <= 20.0) 
        if (leptonOneP4.Pt() < 25.0) 
            return kTRUE;

        //if (leptonTwoP4.Pt() <= 10.0)
        if (leptonTwoP4.Pt() < 20.0)
            return kTRUE;

        if (sync_print_precut) {
            cout << "selected lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << endl;
        }

        TLorentzVector dileptonP4 = leptonOneP4 + leptonTwoP4;

        // trigger matching:
        bool mu1_fired_leg1 = false;
        bool mu1_fired_leg2 = false;
        bool mu2_fired_leg1 = false;
        bool mu2_fired_leg2 = false;
        for (unsigned int iT = 0; iT < triggerNames.size(); ++iT) {
            mu1_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[muonOneIndex]->hltMatchBits);
            mu1_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[muonOneIndex]->hltMatchBits);
            mu2_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[muonTwoIndex]->hltMatchBits);
            mu2_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[muonTwoIndex]->hltMatchBits);
        }
        /*if (sync_print_precut) {
            cout << "mu1_fired_leg1, mu1_fired_leg2, mu2_fired_leg1, mu2_fired_leg2" << endl;
            cout << mu1_fired_leg1 << ", " << mu1_fired_leg2 << ", " << mu2_fired_leg1 << ", " << mu2_fired_leg2 << endl;
        }*/

        // L1EMTF cut 
        if (
            fabs(leptonOneP4.DeltaPhi(leptonTwoP4)) < 70.0*(M_PI/180.0)
            && fabs(leptonOneP4.Eta()) > 1.2 
            && fabs(leptonTwoP4.Eta()) > 1.2
            && leptonOneP4.Eta()*leptonTwoP4.Eta() > 0
           )
            return kTRUE;
        hTotalEvents->Fill(8); 

        if (sync_print_precut) 
            cout << "passed the L1EMTF cut" << endl;

        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            TLorentzVector tempPhoton;
            //TLorentzVector tempDilepton;
            TLorentzVector tempLLG;
            tempPhoton.SetPtEtaPhiM(photons[i]->calibPt, photons[i]->eta, photons[i]->phi, 0.);
            //tempDilepton = leptonOneP4 + leptonTwoP4;
            //tempLLG = leptonOneP4 + leptonTwoP4 + tempPhoton;
            tempLLG = dileptonP4 + tempPhoton;
            float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);
            if (sync_print_precut) {
                cout << "photon pt, et/m, z_mass, h_mass, dr1, dr2" << endl;
                cout << tempPhoton.Pt() << ", " << tempPhoton.Et()/tempLLG.M() << ", " 
                     << dileptonP4.M() << ", " << tempLLG.M() << ", " 
                     << this_dr1 << ", " << this_dr2 << endl;
            }
	    //std::cout << " Photon Pt" << tempPhoton.Pt() << std::endl;
            if (
                tempPhoton.Pt() >= 15.0 
                //&& tempPhoton.Et()/tempLLG.M() > (15.0/110.0) 
                //&& dileptonP4.M() + tempLLG.M() > 185.0 
                //&& tempLLG.M() > 100. && tempLLG.M() < 180. 
                //&& this_dr1 > 0.4 && this_dr2 > 0.4
                ) {
                hasValidPhoton = true;
                photonIndex = i;
                break;
            }
        }
	

        if (!hasValidPhoton)
            return kTRUE;
        hTotalEvents->Fill(9);

        if (sync_print_precut)
            cout << "passed the photon selection" << endl;

        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->calibPt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(10);

        TLorentzVector llgP4 = dileptonP4 + photonOneP4;

        // DY photon overlap removal
        vetoDY = false;
        for (unsigned int i = 0; i < genPhotons.size(); ++i) {
            TGenParticle *pho = genPhotons.at(i);
            if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
		if(pho->pdgId == 22){
			TLorentzVector thisGenPhotonP4;
			thisGenPhotonP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
			if (thisGenPhotonP4.DeltaR(photonOneP4) < 0.1) {
			    vetoDY = true;
			    break;
			}
		}
            }
	    // gen Iso
	    float genIso = GetGenIsolation(pho);
	    if(genIso < 5) 
		genIsoPass = true;
        }

 
        // checking for lepton tag
        isLeptonTag = false;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            TLorentzVector tempMuon;
            tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
            if (leptonOneP4.DeltaR(tempMuon) < 0.4 ||
                leptonTwoP4.DeltaR(tempMuon) < 0.4 || 
                photonOneP4.DeltaR(tempMuon) < 0.4)
                continue;
            else
                isLeptonTag = true;
        }

        if (!isLeptonTag) {
            for (unsigned int i = 0; i < electrons.size(); ++i) {
                TLorentzVector tempElectron;
                tempElectron.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (leptonOneP4.DeltaR(tempElectron) < 0.4 ||
                    leptonTwoP4.DeltaR(tempElectron) < 0.4 || 
                    photonOneP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }

	//////////////////////////////////////////////
        for (unsigned int i = 0; i < muons.size(); ++i) {
		TLorentzVector tempMuon;
		tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
		if (leptonOneP4.DeltaR(tempMuon) < 0.1)
		    	leptonOneTag = true;
		if (leptonTwoP4.DeltaR(tempMuon) < 0.1)
		    	leptonTwoTag = true;
	}
	///////////////////////////////////
        
        // checking for dijet tag
        isDijetTag = false;
        isTightDijetTag = false;
        unsigned int jetOneIndex = 0;
        unsigned int jetTwoIndex = 0;
        unsigned int purityLevel = 0;
        TLorentzVector jetOneP4, jetTwoP4;
        //TLorentzVector llg = leptonOneP4 + leptonTwoP4 + photonOneP4;
        if (!isLeptonTag) {
            /*if (sync_print_precut) 
                cout << "event is not lepton tagged; checking dijet tag" << endl;*/
            if (jets.size() > 1)  {
                /*if (sync_print_precut)
                    cout << "event has at least two jets" << endl;*/
                for (unsigned int i = 0; i < jets.size(); ++i) {
                    for (unsigned int j = i+1; j < jets.size(); ++j) {
                        TLorentzVector tempJetOne;
                        TLorentzVector tempJetTwo;
                        tempJetOne.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
                        tempJetTwo.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
                        TLorentzVector tempDijet = tempJetOne + tempJetTwo;
                        float zeppen = llgP4.Eta() - (tempJetOne.Eta() + tempJetTwo.Eta())/2.;
                        if (    purityLevel < 5  
                            &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                            &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                            &&  fabs(zeppen) <= 2.5 && tempDijet.M() >= 500.
                            &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                           ) {
                            isDijetTag = true;
                            isTightDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 5;
                            break;
                        }
                        else if (   purityLevel < 4 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(zeppen) <= 2.5
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 4;
                        }
                        else if (   purityLevel < 3 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 3;
                        }
                        else if (   purityLevel < 2
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 2;
                        }
                        else if (   purityLevel < 1
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 1;
                        }
                                
                    }
                }
            }
        }

        if (isDijetTag) {
            jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
            jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
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
            
            zepp = llgP4.Eta() - (jetOneP4.Eta() + jetTwoP4.Eta())/2.;
            llgJJDEta = fabs(llgP4.Eta() - dijet.Eta());
            llgJJDPhi = fabs(llgP4.DeltaPhi(dijet));
            llgJJDR = llgP4.DeltaR(dijet);

        }
 
        if (sync_print_precut) {
            cout << "lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi, dilepton_mass, dr1, dr2" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << ", "
                 << (leptonOneP4 + leptonTwoP4).M() << ", " << leptonOneP4.DeltaR(photonOneP4) << ", " << leptonTwoP4.DeltaR(photonOneP4) << endl;
        } 

	if(debugSelect)
		cout << "------- Fill -------------\n";
	if(debugSelect)
		cout << "--- Lepton Fill" << endl;
        leptonOnePt     = leptonOneP4.Pt();
        leptonOneEta    = leptonOneP4.Eta();
        leptonOnePhi    = leptonOneP4.Phi();
        leptonOneIso    = GetMuonIsolation(muons[muonOneIndex]);
        leptonOneFlavor = muons[muonOneIndex]->q*13;
        leptonOneDZ     = muons[muonOneIndex]->dz;
        leptonOneD0     = muons[muonOneIndex]->d0;
        leptonOneCharge = muons[muonOneIndex]->q;
            
        leptonTwoPt     = leptonTwoP4.Pt();
        leptonTwoEta    = leptonTwoP4.Eta();
        leptonTwoPhi    = leptonTwoP4.Phi();
        leptonTwoIso    = GetMuonIsolation(muons[muonTwoIndex]);
        leptonTwoFlavor = muons[muonTwoIndex]->q*13;
        leptonTwoDZ     = muons[muonTwoIndex]->dz;
        leptonTwoD0     = muons[muonTwoIndex]->d0;
        leptonTwoCharge = muons[muonTwoIndex]->q;

	if(debugSelect)
		cout << "--- Photon Fill" << endl;
        photonOnePt     = photonOneP4.Pt();
        photonOneEta    = photonOneP4.Eta();
        photonOnePhi    = photonOneP4.Phi();
        if(params->period == "2016ReReco")
		photonOneMVA    = photons[photonIndex]->mvaSpring16;
        else if( params->period == "2016Legacy")
		photonOneMVA    = photons[photonIndex]->mvaFall17V2;
	else
		photonOneMVA    = photons[photonIndex]->mvaFall17V2;
        photonOneERes   = photons[photonIndex]->eRes;
        photonOneSieie  = photons[photonIndex]->sieie;
        photonOneHoverE = photons[photonIndex]->hovere;
        photonOneIph    = photons[photonIndex]->gammaIso;
        photonOneIneu   = photons[photonIndex]->neuHadIso;
        photonOneIch    = photons[photonIndex]->chHadIso;

	photonOneSieip       = photons[photonIndex]->sieip;
	photonOneSipip      = photons[photonIndex]->sipip;
	photonOneSrr        = photons[photonIndex]->srr;
	photonOneE2x2       = photons[photonIndex]->e2x2;
	photonOneE5x5       = photons[photonIndex]->e5x5;
	photonOneScEtaWidth = photons[photonIndex]->scEtaWidth;
	photonOneScPhiWidth = photons[photonIndex]->scPhiWidth;
	photonOneScRawE     = photons[photonIndex]->scRawE;
	photonOnePreShowerE = photons[photonIndex]->scESEn;
	photonOneScBrem     = photons[photonIndex]->scBrem;
	
	/*                    photons[photonIndex]->
	leptonOnePhotonOneDEta = fabs(leptonOneP4.Eta() - photonOneP4.Eta());
	leptonOnePhotonOneDPhi = fabs(leptonOneP4.DeltaPhi(photonOneP4));
	leptonOnePhotonOneDR   =      leptonOneP4.DeltaR(photonOneP4);

	leptonTwoPhotonOneDEta = fabs(leptonTwoP4.Eta() - photonOneP4.Eta());
	leptonTwoPhotonOneDPhi = fabs(leptonTwoP4.DeltaPhi(photonOneP4));
	leptonTwoPhotonOneDR   =      leptonTwoP4.DeltaR(photonOneP4);
	*/


	/*
	cout << "------- Fill Gen -------------\n";
	genLeptonOnePt  = genLeptons[0]->pt;
	genLeptonOneEta = genLeptons[0]->eta;
	genLeptonOnePhi = genLeptons[0]->phi;
	genLeptonTwoPt  = genLeptons[1]->pt;
	genLeptonTwoEta = genLeptons[1]->eta;
	genLeptonTwoPhi = genLeptons[1]->phi;
	genPhotonPt  = genPhotons[0]->pt;
	genPhotonEta = genPhotons[0]->eta;
	genPhotonPhi = genPhotons[0]->phi;
	*/	

        passElectronVeto = photons[photonIndex]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
        else 
            photonOneR9 = photons[photonIndex]->r9;

        if (sync_print_precut) {
            cout << "event is still alive" << endl;
        }
        
	if(debugSelect)
		cout << "--- Dilep Fill" << endl;
        dileptonPt  = dileptonP4.Pt();
        dileptonEta = dileptonP4.Eta();
        dileptonPhi = dileptonP4.Phi();
        dileptonM   = dileptonP4.M();
        //dileptonMKin = (muonOneP4KinFit + muonTwoP4KinFit).M();
        //dileptonMKinJames = (muonOneP4KinFitJames + muonTwoP4KinFitJames).M();
        dileptonDEta = fabs(leptonOneP4.Eta() - leptonTwoP4.Eta());
        dileptonDPhi = fabs(leptonOneP4.DeltaPhi(leptonTwoP4));
        dileptonDR = leptonOneP4.DeltaR(leptonTwoP4);

        llgPt  = llgP4.Pt();
        llgEta = llgP4.Eta();
        llgPhi = llgP4.Phi();
        llgM   = llgP4.M();
        llgPtOverM = llgP4.Pt()/llgP4.M();
        //llgMKin = (muonOneP4KinFit + muonTwoP4KinFit + photonOneP4).M();
        //llgMKinJames = (muonOneP4KinFitJames + muonTwoP4KinFitJames + photonOneP4).M();
        
        l1PhotonDEta = fabs(leptonOneP4.Eta() - photonOneP4.Eta());
        l1PhotonDPhi = fabs(leptonOneP4.DeltaPhi(photonOneP4));
        l1PhotonDR   = leptonOneP4.DeltaR(photonOneP4);
        l1PhotonM    = (leptonOneP4 + photonOneP4).M();
        l1PhotonPt   = (leptonOneP4 + photonOneP4).Pt();

        l2PhotonDEta = fabs(leptonTwoP4.Eta() - photonOneP4.Eta());
        l2PhotonDPhi = fabs(leptonTwoP4.DeltaPhi(photonOneP4));
        l2PhotonDR   = leptonTwoP4.DeltaR(photonOneP4);
        l2PhotonM    = (leptonTwoP4 + photonOneP4).M();
        l2PhotonPt   = (leptonTwoP4 + photonOneP4).Pt();
        
        if (l1PhotonDR > l2PhotonDR) {
            lPhotonDRMax = l1PhotonDR;
            lPhotonDRMin = l2PhotonDR; 
        }
        else {
            lPhotonDRMax = l2PhotonDR;
            lPhotonDRMin = l1PhotonDR;
        }

        dileptonPhotonDEta = fabs(dileptonP4.Eta() - photonOneP4.Eta());
        dileptonPhotonDPhi = fabs(dileptonP4.DeltaPhi(photonOneP4));
        dileptonPhotonDR = dileptonP4.DeltaR(photonOneP4);
        ptt = 2*fabs(dileptonP4.Px()*photonOneP4.Py() - photonOneP4.Px()*dileptonP4.Py())/llgP4.Pt();

        /*std::cout << "kinematics before Brian angles" << std::endl;
        leptonOneP4.Print();
        leptonTwoP4.Print();
        dileptonP4.Print();
        llgP4.Print();*/
        // calculate angles like Brian
        TVector3 Xframe = llgP4.BoostVector();
        TVector3 Z1frame = dileptonP4.BoostVector();

        // "partons"
        TLorentzVector kq, kqbar, veckq_in_Xframe, veckqbar_in_Xframe;
        kq.SetPxPyPzE(0., 0., (llgP4.E() + llgP4.Pz())/2., (llgP4.E() + llgP4.Pz())/2.);
        kqbar.SetPxPyPzE(0., 0., (llgP4.Pz() - llgP4.E())/2., (llgP4.E() - llgP4.Pz())/2.);
        veckq_in_Xframe = kq;
        veckqbar_in_Xframe = kqbar;
        veckq_in_Xframe.Boost(-1*Xframe);
        veckqbar_in_Xframe.Boost(-1*Xframe);
   
        // Z vectors
        TLorentzVector vecz_in_Xframe = dileptonP4;
        TLorentzVector vecg_in_Xframe = photonOneP4;
        TLorentzVector vecz_in_Z1frame = dileptonP4;
        vecz_in_Xframe.Boost(-1*Xframe);
        vecg_in_Xframe.Boost(-1*Xframe);
        vecz_in_Z1frame.Boost(-1*Z1frame);

        // coord system in the CM frame
        TVector3 uz_in_Xframe = vecz_in_Xframe.Vect().Unit();
        TVector3 uy_in_Xframe = (veckq_in_Xframe.Vect().Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TVector3 ux_in_Xframe = (uy_in_Xframe.Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TRotation rotation;
        rotation = rotation.RotateAxes(ux_in_Xframe, uy_in_Xframe, uz_in_Xframe).Inverse();

        // for going to the Z frames from the CM frame, boost after transform
        TLorentzVector vecz_in_Xframe_newcoords = vecz_in_Xframe;
        vecz_in_Xframe_newcoords.Transform(rotation);
        TVector3 Z1frame_from_Xframe_newcoords = vecz_in_Xframe_newcoords.BoostVector();

        // define the positive and negative leptons
        TLorentzVector l_minus_james, l_plus_james; 
        if (leptonOneFlavor > 0) {
            l_minus_james = leptonOneP4;
            l_plus_james = leptonTwoP4;
        }
        else {
            l_minus_james = leptonTwoP4;
            l_plus_james = leptonOneP4;
        }
       
        // little theta, phi in Z1 frame; first boost to CM, then redefine coords
        TLorentzVector veclm_in_Z1frame = l_minus_james;
        TLorentzVector veclp_in_Z1frame = l_plus_james;
        veclm_in_Z1frame.Boost(-1*Xframe);
        veclm_in_Z1frame.Transform(rotation);
        veclp_in_Z1frame.Boost(-1*Xframe);
        veclp_in_Z1frame.Transform(rotation);

        // then boost to Z1
        veclm_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
        veclp_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);

        // now get angles
        zgPhiJames = veclm_in_Z1frame.Phi();
        zgLittleThetaJames = veclm_in_Z1frame.CosTheta();

        if (zgPhiJames < 0) 
            zgPhiJames += 2*M_PI;

        // Big Theta in X frame
        TLorentzVector veczg_in_Xframe = llgP4;
        veczg_in_Xframe.Transform(rotation);

        TLorentzVector veczg_in_Xframe_newcoords = llgP4;
        veczg_in_Xframe_newcoords.Transform(rotation);
        zgBigThetaJames = (-1*veczg_in_Xframe_newcoords.Vect()).CosTheta();

        /////////////////////////////
        //std::cout << "kinematics after Brian angles" << std::endl;
        //leptonOneP4.Print();
        //leptonTwoP4.Print();
        //dileptonP4.Print();
        //llgP4.Print();

        // calculate angles like Ming-Yan
        TLorentzVector l_minus, l_plus; 
        if (leptonOneFlavor > 0) {
            l_minus = leptonOneP4;
            l_plus = leptonTwoP4;
        }
        else {
            l_minus = leptonTwoP4;
            l_plus = leptonOneP4;
        }
        
        TVector3 llgFrame = -1*llgP4.BoostVector();
        dileptonP4.Boost(llgFrame);
        l_minus.Boost(llgFrame);
        l_minus.Boost(dileptonP4.BoostVector());
        zgLittleTheta = cos(dileptonP4.Angle(l_minus.Vect()));
        zgBigTheta = cos(dileptonP4.Angle(llgP4.Vect()));
       
        // the way MY does it (I think wrong)
        TLorentzVector lep0 = leptonOneP4;
        lep0.Boost(llgFrame);
        lep0.Boost(dileptonP4.BoostVector());
        zgLittleThetaMY = cos(dileptonP4.Angle(lep0.Vect()));
        
        TVector3 ppAxis(0, 0, 1);
        TVector3 zAxis = dileptonP4.Vect().Unit();
        TVector3 yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rot;
        rot = rot.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4.Transform(rot);
        l_minus.Transform(rot);
        zgPhi = l_minus.Phi();
           
	if(debugSelect) 
		cout << "--- Weights Implemented" << endl;
        if (!isData) {
	    ///      ID EFFICIENTCY
	    //
            //muonIDWeightOne = weights->GetHZZMuonIDEff(*muons[muonOneIndex]); 
            //muonIDWeightTwo = weights->GetHZZMuonIDEff(*muons[muonTwoIndex]);
            if(debugSelect)
		cout << "------ Muon ID Weights" << endl;
	
            muonIDWeightOne = weights->GetMuonIDEff(leptonOneP4); 
            muonIDWeightTwo = weights->GetMuonIDEff(leptonTwoP4);
            eventWeight *= muonIDWeightOne;
            eventWeight *= muonIDWeightTwo;

	    ///      ISO EFFICIENTCY
            if(debugSelect)
		cout << "------ Muon ISO Weights" << endl;
	  
	    muonISOWeightOne = weights->GetMuonISOEff(leptonOneP4); 
	    muonISOWeightTwo = weights->GetMuonISOEff(leptonTwoP4);
            eventWeight *= muonISOWeightOne;
            eventWeight *= muonISOWeightTwo;
	    //////////////////////////// 
            float sf11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[muonOneIndex]);
            float sf12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[muonTwoIndex]);
            float sf21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[muonOneIndex]);
            float sf22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[muonTwoIndex]);
            
            if (leptonTwoPt < 20.) {
                muonTrigWeightOne = sf11;
                muonTrigWeightTwo = sf22;
            }
            else {
                float prod1 = sf11*sf22;
                float prod2 = sf21*sf12;
                if (prod1 > prod2) {
                    muonTrigWeightOne = sf11;
                    muonTrigWeightTwo = sf22;
                }
                else {
                    muonTrigWeightOne = sf21;
                    muonTrigWeightTwo = sf12;
                }
            }

            if(debugSelect)
		cout << "------ Muon Trigger Weights" << endl;
            triggerWeight = muonTrigWeightOne*muonTrigWeightTwo;
            eventWeight *= triggerWeight;



            if(debugSelect)
		cout << "------ Photon ID Weights" << endl;

	    if(params->period == "2016Legacy"){
		photonIDWeight = weights->GetPhotonMVAIdEff(*photons[photonIndex]); 
	    }
	    else if(params->period == "2017ReReco"){
		photonIDWeight = weights->GetPhotonIdEff(*photons[photonIndex]);
	    }
	    eventWeight *= photonIDWeight; 


	    if(params->period == "2016Legacy"){
		//photonIsConvWeight = weights->GetPhotonMVAIdEff(*photons[photonIndex]); 
		photonIsConvWeight = 1; 
	    }
	    else if(params->period == "2017ReReco"){
		photonIsConvWeight = weights->GetPhotonIsConvEff(*photons[photonIndex]); 
	    }
	    eventWeight *= photonIsConvWeight; 

	    cout << "Photon ID     Weight:: " << photonIDWeight     << endl
	    	 << "Photon IsConv Weight:: " << photonIsConvWeight << endl;

            //eventWeight *= photonIDWeight;
            if (sync_print_precut) {
                cout << "run,lumi,evt,puwei,totSF,trg0,trg1,id0,id1,gammaID,pt0,pt1,eta0,et1" << endl;
                cout << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << ", "
                     << puWeight << ", " << eventWeight << ", " << muonTrigWeightOne << ", "
                     << muonTrigWeightTwo << ", " << muonIDWeightOne << ", " << muonIDWeightTwo << ", "
                     << photonIDWeight << ", " << muons[muonOneIndex]->pt << ", " << muons[muonTwoIndex]->pt << ", " 
                     << muons[muonOneIndex]->eta << ", " << muons[muonTwoIndex]->eta << endl;
            }
        }

        if (sync_print_precut) {
            cout << "event still alive after event weights" << endl;
        }

	if(debugSelect) 
		cout << "--- Finish Selection" << endl;
    } // end mumug selection
    
    /*else if (params->selection == "ee") {
        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);
        TLorentzVector electronOneP4, electronTwoP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);
        if (electronOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);
        if (electronTwoP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(7);
        TLorentzVector dielectron = electronOneP4 + electronTwoP4;
        if (electrons[0]->q == electrons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        if (dielectron.M() < 80.0 || dielectron.M() > 100.0)
            return kTRUE;
        hTotalEvents->Fill(9);   
        
        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;
            
        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;
           
        if (!isData) {

            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[1]); 

            pair<float, float> eff11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[0]);
            pair<float, float> eff12 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[1]);
            pair<float, float> eff21 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[0]);
            pair<float, float> eff22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;

        }
    } */
    
    else if (params->selection == "elelg") {
	if(debugSelect){
		cout << "Start selection...\n";
		cout << "-------nEL:: "<< nElectrons << endl;
		cout << "-------nPh:: "<< nPhotons << endl;
		}
        if (electrons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
	if(debugSelect)
		cout << "-------Electron Mult Pass ---------\n";
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
	if(debugSelect)
		cout << "-------Photon Mult Pass ---------\n";
      
	if(debugSelect)
		cout << "Multiplicity selection...\n";
	//cout << "-------Multiplicity Pass ---------\n";
        TLorentzVector leptonOneP4, leptonTwoP4;
        unsigned int electronOneIndex = 0;
        unsigned int electronTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < electrons.size(); ++i) {
            for (unsigned int j = i+1; j < electrons.size(); ++j) {
                TLorentzVector tempElectronOne, tempElectronTwo;
                tempElectronOne.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                tempElectronTwo.SetPtEtaPhiM(electrons[j]->calibPt, electrons[j]->eta, electrons[j]->phi, ELE_MASS);
                float thisMass = (tempElectronOne + tempElectronTwo).M();
                if (sync_print_precut) {
                    cout << "Z pair loop for i = " << i << ", j = " << j << endl;
                    cout << "current value of zMassDiff = " << zMassDiff << endl;
                    cout << "lepton i pt, eta, phi: " << tempElectronOne.Pt() << ", " << tempElectronOne.Eta() << ", " 
                         << tempElectronOne.Phi() << endl;
                    cout << "lepton j pt, eta, phi: " << tempElectronTwo.Pt() << ", " << tempElectronTwo.Eta() << ", " 
                         << tempElectronTwo.Phi() << endl;
                    cout << "candidate Z mass = " << thisMass << endl;
                    cout << "lepton i, j q: " << electrons[i]->q << ", " << electrons[j]->q << endl;
                } 
                //if (
                //        //electrons[i]->q != electrons[j]->q
                //        /*&&*/ electrons[i]->calibPt > 25.0 
                //        //&& electrons[i]->calibPt > 30.0 
                //        && electrons[j]->calibPt > 15.0
                //        && thisMass > 50.0
                //   ) {
                if (thisMass > 50.0) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempElectronOne;
                            leptonTwoP4 = tempElectronTwo;
                            electronOneIndex = i;
                            electronTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempElectronOne;
                        leptonTwoP4 = tempElectronTwo;
                        electronOneIndex = i;
                        electronTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }

        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (leptonOneP4.Pt() < 25.0)
            return kTRUE;

        if (leptonTwoP4.Pt() < 20.0)
            return kTRUE;
        
	if(debugSelect)
		cout << "-------Lepton Pass --------\n";
        if (sync_print_precut) {
            cout << "selected lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << endl;
        }
        
        TLorentzVector dileptonP4 = leptonOneP4 + leptonTwoP4;
       
        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            TLorentzVector tempPhoton;
            TLorentzVector tempLLG;
            tempPhoton.SetPtEtaPhiM(photons[i]->calibPt, photons[i]->eta, photons[i]->phi, 0.);
            tempLLG = dileptonP4 + tempPhoton;
            float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);
            if (sync_print_precut) {
                cout << "photon pt, et/m, z_mass, h_mass, dr1, dr2" << endl;
                cout << tempPhoton.Pt() << ", " << tempPhoton.Et()/tempLLG.M() << ", " 
                     << dileptonP4.M() << ", " << tempLLG.M() << ", " 
                     << this_dr1 << ", " << this_dr2 << endl;
            }
            if (
                tempPhoton.Pt() >= 15.0 
                //&& tempPhoton.Et()/tempLLG.M() > (15.0/110.0) 
                //&& dileptonP4.M() + tempLLG.M() > 185.0 
                //&& tempLLG.M() > 100. && tempLLG.M() < 180. 
                //&& this_dr1 > 0.4 && this_dr2 > 0.4
                ) {
                hasValidPhoton = true;
                photonIndex = i;
                break;
            }
            else {
                if (sync_print_precut) {
                    cout << "photon variables in the loop" << endl;
                    cout << "pt, et/m_llg, mass_sum, llg_mass, dr1, dr2" << endl;
                    cout << tempPhoton.Pt() << ", " << tempPhoton.Et()/tempLLG.M() << ", " << 
                            dileptonP4.M() + tempLLG.M() << ", " << tempLLG.M() << ", " << 
                            this_dr1 << ", " << this_dr2 << endl;
                }
            }
        }

        if (!hasValidPhoton)
            return kTRUE; 
        hTotalEvents->Fill(8);
        if (sync_print_precut)
            cout << "valid photon" << endl;

        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->calibPt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(9);
	if(debugSelect)
		cout << "-------Photon Pass --------\n";
        
        TLorentzVector llgP4 = dileptonP4 + photonOneP4;
        
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
	    // gen Iso
	    float genIso = GetGenIsolation(pho);
	    if(genIso < 5) 
		genIsoPass = true;
        }  
        
        // checking for lepton tag
        isLeptonTag = false;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            TLorentzVector tempMuon;
            tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
            if (leptonOneP4.DeltaR(tempMuon) < 0.4 ||
                leptonTwoP4.DeltaR(tempMuon) < 0.4 || 
                photonOneP4.DeltaR(tempMuon) < 0.4)
                continue;
            else
                isLeptonTag = true;
        }

        if (!isLeptonTag) {
            for (unsigned int i = 0; i < electrons.size(); ++i) {
                TLorentzVector tempElectron;
                tempElectron.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (leptonOneP4.DeltaR(tempElectron) < 0.4 ||
                    leptonTwoP4.DeltaR(tempElectron) < 0.4 || 
                    photonOneP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }
       

	//////////////////////////////////////////////
        for (unsigned int i = 0; i < electrons.size(); ++i) {
		TLorentzVector tempElectron;
		tempElectron.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
		if (leptonOneP4.DeltaR(tempElectron) < 0.1)
		    	leptonOneTag = true;
		if (leptonTwoP4.DeltaR(tempElectron) < 0.1)
		    	leptonTwoTag = true;
	}
	///////////////////////////////////
 
        // checking for dijet tag
        isDijetTag = false;
        isTightDijetTag = false;
        unsigned int jetOneIndex = 0;
        unsigned int jetTwoIndex = 0;
        unsigned int purityLevel = 0;
        TLorentzVector jetOneP4, jetTwoP4;
        if (!isLeptonTag) {
            //if (sync_print_precut) 
             //   cout << "event is not lepton tagged; checking dijet tag" << endl;
            if (jets.size() > 1)  {
                //if (sync_print_precut)
                 //   cout << "event has at least two jets" << endl;
                for (unsigned int i = 0; i < jets.size(); ++i) {
                    for (unsigned int j = i+1; j < jets.size(); ++j) {
                        TLorentzVector tempJetOne;
                        TLorentzVector tempJetTwo;
                        tempJetOne.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
                        tempJetTwo.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
                        TLorentzVector tempDijet = tempJetOne + tempJetTwo;
                        float zeppen = llgP4.Eta() - (tempJetOne.Eta() + tempJetTwo.Eta())/2.;
                        if (    purityLevel < 5  
                            &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                            &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                            &&  fabs(zeppen) <= 2.5 && tempDijet.M() >= 500.
                            &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                           ) {
                            isDijetTag = true;
                            isTightDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 5;
                            break;
                        }
                        else if (   purityLevel < 4 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(zeppen) <= 2.5
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 4;
                        }
                        else if (   purityLevel < 3 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 3;
                        }
                        else if (   purityLevel < 2
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 2;
                        }
                        else if (   purityLevel < 1
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 1;
                        }
                                
                    }
                }
            }
        }

        
        if (isDijetTag) {
            jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
            jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
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
            
            zepp = llgP4.Eta() - (jetOneP4.Eta() + jetTwoP4.Eta())/2.;
            llgJJDEta = fabs(llgP4.Eta() - dijet.Eta());
            llgJJDPhi = fabs(llgP4.DeltaPhi(dijet));
            llgJJDR = llgP4.DeltaR(dijet);
        }
        
        if (sync_print_precut) {
            cout << "lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi, dilepton_mass, dr1, dr2" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << ", "
                 << (leptonOneP4 + leptonTwoP4).M() << ", " << leptonOneP4.DeltaR(photonOneP4) << ", " << leptonTwoP4.DeltaR(photonOneP4) << endl;
        }

	if(debugSelect)
		cout << "------- Fill -------------\n";
        leptonOnePt         = leptonOneP4.Pt();
        leptonOneEta        = leptonOneP4.Eta();
        leptonOnePhi        = leptonOneP4.Phi();
        leptonOneIso        = GetElectronIsolation(electrons[electronOneIndex], fInfo->rhoJet);
        leptonOneFlavor     = electrons[electronOneIndex]->q*11;
        leptonOneDZ         = electrons[electronOneIndex]->dz;
        leptonOneD0         = electrons[electronOneIndex]->d0;
        leptonOneECALDriven = (electrons[electronOneIndex]->typeBits & baconhep::kEcalDriven);
	leptonOneCharge     = electrons[electronOneIndex]->q;        
    
        leptonTwoPt         = leptonTwoP4.Pt();
        leptonTwoEta        = leptonTwoP4.Eta();
        leptonTwoPhi        = leptonTwoP4.Phi();
        leptonTwoIso        = GetElectronIsolation(electrons[electronTwoIndex], fInfo->rhoJet);
        leptonTwoFlavor     = electrons[electronTwoIndex]->q*11;
        leptonTwoDZ         = electrons[electronTwoIndex]->dz;
        leptonTwoD0         = electrons[electronTwoIndex]->d0;
        leptonTwoECALDriven = (electrons[electronTwoIndex]->typeBits & baconhep::kEcalDriven);
	leptonTwoCharge     = electrons[electronTwoIndex]->q;        

        photonOnePt    = photonOneP4.Pt();
        photonOneEta   = photonOneP4.Eta();
        photonOnePhi   = photonOneP4.Phi();
	if(params->period == "2016")
	        photonOneMVA   = photons[photonIndex]->mvaSpring16;
	else if(params->period == "2016Legacy")
	        photonOneMVA   = photons[photonIndex]->mvaFall17V2;
	else
	        photonOneMVA   = photons[photonIndex]->mvaFall17V2;
        photonOneERes  = photons[photonIndex]->eRes;
        photonOneSieie = photons[photonIndex]->sieie;
        photonOneIph   = photons[photonIndex]->gammaIso;
        photonOneIneu  = photons[photonIndex]->neuHadIso;
        photonOneIch   = photons[photonIndex]->chHadIso;

	photonOneSieip       = photons[photonIndex]->sieip;
	photonOneSipip      = photons[photonIndex]->sipip;
	photonOneSrr        = photons[photonIndex]->srr;
	photonOneE2x2       = photons[photonIndex]->e2x2;
	photonOneE5x5       = photons[photonIndex]->e5x5;
	photonOneScEtaWidth = photons[photonIndex]->scEtaWidth;
	photonOneScPhiWidth = photons[photonIndex]->scPhiWidth;
	photonOneScRawE     = photons[photonIndex]->scRawE;
	photonOnePreShowerE = photons[photonIndex]->scESEn;
	photonOneScBrem     = photons[photonIndex]->scBrem;


	/*
	cout << "------- Fill GEN -------------\n";
	genLeptonOnePt  = genLeptons[0]->pt;
	genLeptonOneEta = genLeptons[0]->eta;
	genLeptonOnePhi = genLeptons[0]->phi;
	genLeptonTwoPt  = genLeptons[1]->pt;
	genLeptonTwoEta = genLeptons[1]->eta;
	genLeptonTwoPhi = genLeptons[1]->phi;
	genPhotonPt  = genPhotons[0]->pt;
	genPhotonEta = genPhotons[0]->eta;
	genPhotonPhi = genPhotons[0]->phi;
	*/

        passElectronVeto = photons[photonIndex]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
        else 
            photonOneR9 = photons[photonIndex]->r9;

        //if (sync_print_precut)
         //   cout << "photon r9 = " << photonOneR9 << endl;



        dileptonPt = dileptonP4.Pt();
        dileptonEta = dileptonP4.Eta();
        dileptonPhi = dileptonP4.Phi();
        dileptonM = dileptonP4.M();
        //dileptonMKin = (electronOneP4KinFit + electronTwoP4KinFit).M();
        //dileptonMKinJames = (electronOneP4KinFitJames + electronTwoP4KinFitJames).M();
        dileptonDEta = fabs(leptonOneP4.Eta() - leptonTwoP4.Eta());
        dileptonDPhi = fabs(leptonOneP4.DeltaPhi(leptonTwoP4));
        dileptonDR = leptonOneP4.DeltaR(leptonTwoP4);

        llgPt = llgP4.Pt();
        llgEta = llgP4.Eta();
        llgPhi = llgP4.Phi();
        llgM = llgP4.M();
        llgPtOverM = llgP4.Pt()/llgP4.M();
        //llgMKin = (electronOneP4KinFit + electronTwoP4KinFit + photonOneP4).M();
        //llgMKinJames = (electronOneP4KinFitJames + electronTwoP4KinFitJames + photonOneP4).M();
        
	//cout << "------- AFTER FIT ------\n";

        l1PhotonDEta = fabs(leptonOneP4.Eta() - photonOneP4.Eta());
        l1PhotonDPhi = fabs(leptonOneP4.DeltaPhi(photonOneP4));
        l1PhotonDR   = leptonOneP4.DeltaR(photonOneP4);
        l1PhotonM    = (leptonOneP4 + photonOneP4).M();
        l1PhotonPt   = (leptonOneP4 + photonOneP4).Pt();

        l2PhotonDEta = fabs(leptonTwoP4.Eta() - photonOneP4.Eta());
        l2PhotonDPhi = fabs(leptonTwoP4.DeltaPhi(photonOneP4));
        l2PhotonDR   = leptonTwoP4.DeltaR(photonOneP4);
        l2PhotonM    = (leptonTwoP4 + photonOneP4).M();
        l2PhotonPt   = (leptonTwoP4 + photonOneP4).Pt();

        if (l1PhotonDR > l2PhotonDR) {
            lPhotonDRMax = l1PhotonDR;
            lPhotonDRMin = l2PhotonDR; 
        }
        else {
            lPhotonDRMax = l2PhotonDR;
            lPhotonDRMin = l1PhotonDR;
        }

        dileptonPhotonDEta = fabs(dileptonP4.Eta() - photonOneP4.Eta());
        dileptonPhotonDPhi = fabs(dileptonP4.DeltaPhi(photonOneP4));
        dileptonPhotonDR = dileptonP4.DeltaR(photonOneP4);
        ptt = 2*fabs(dileptonP4.Px()*photonOneP4.Py() - photonOneP4.Px()*dileptonP4.Py())/llgP4.Pt();
        
        /*std::cout << "kinematics before Brian angles" << std::endl;
        leptonOneP4.Print();
        leptonTwoP4.Print();
        dileptonP4.Print();
        llgP4.Print();*/
        // calculate angles like Brian
        TVector3 Xframe = llgP4.BoostVector();
        TVector3 Z1frame = dileptonP4.BoostVector();

        // "partons"
        TLorentzVector kq, kqbar, veckq_in_Xframe, veckqbar_in_Xframe;
        kq.SetPxPyPzE(0., 0., (llgP4.E() + llgP4.Pz())/2., (llgP4.E() + llgP4.Pz())/2.);
        kqbar.SetPxPyPzE(0., 0., (llgP4.Pz() - llgP4.E())/2., (llgP4.E() - llgP4.Pz())/2.);
        veckq_in_Xframe = kq;
        veckqbar_in_Xframe = kqbar;
        veckq_in_Xframe.Boost(-1*Xframe);
        veckqbar_in_Xframe.Boost(-1*Xframe);
   
        // Z vectors
        TLorentzVector vecz_in_Xframe = dileptonP4;
        TLorentzVector vecg_in_Xframe = photonOneP4;
        TLorentzVector vecz_in_Z1frame = dileptonP4;
        vecz_in_Xframe.Boost(-1*Xframe);
        vecg_in_Xframe.Boost(-1*Xframe);
        vecz_in_Z1frame.Boost(-1*Z1frame);

        // coord system in the CM frame
        TVector3 uz_in_Xframe = vecz_in_Xframe.Vect().Unit();
        TVector3 uy_in_Xframe = (veckq_in_Xframe.Vect().Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TVector3 ux_in_Xframe = (uy_in_Xframe.Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TRotation rotation;
        rotation = rotation.RotateAxes(ux_in_Xframe, uy_in_Xframe, uz_in_Xframe).Inverse();

        // for going to the Z frames from the CM frame, boost after transform
        TLorentzVector vecz_in_Xframe_newcoords = vecz_in_Xframe;
        vecz_in_Xframe_newcoords.Transform(rotation);
        TVector3 Z1frame_from_Xframe_newcoords = vecz_in_Xframe_newcoords.BoostVector();

        // define the positive and negative leptons
        TLorentzVector l_minus_james, l_plus_james; 
        if (leptonOneFlavor > 0) {
            l_minus_james = leptonOneP4;
            l_plus_james = leptonTwoP4;
        }
        else {
            l_minus_james = leptonTwoP4;
            l_plus_james = leptonOneP4;
        }
       
        // little theta, phi in Z1 frame; first boost to CM, then redefine coords
        TLorentzVector veclm_in_Z1frame = l_minus_james;
        TLorentzVector veclp_in_Z1frame = l_plus_james;
        veclm_in_Z1frame.Boost(-1*Xframe);
        veclm_in_Z1frame.Transform(rotation);
        veclp_in_Z1frame.Boost(-1*Xframe);
        veclp_in_Z1frame.Transform(rotation);

        // then boost to Z1
        veclm_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
        veclp_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);

        // now get angles
        zgPhiJames = veclm_in_Z1frame.Phi();
        zgLittleThetaJames = veclm_in_Z1frame.CosTheta();

        if (zgPhiJames < 0) 
            zgPhiJames += 2*M_PI;

        // Big Theta in X frame
        TLorentzVector veczg_in_Xframe = llgP4;
        veczg_in_Xframe.Transform(rotation);

        TLorentzVector veczg_in_Xframe_newcoords = llgP4;
        veczg_in_Xframe_newcoords.Transform(rotation);
        zgBigThetaJames = (-1*veczg_in_Xframe_newcoords.Vect()).CosTheta();

        /////////////////////////////
        //std::cout << "kinematics after Brian angles" << std::endl;
        //leptonOneP4.Print();
        //leptonTwoP4.Print();
        //dileptonP4.Print();
        //llgP4.Print();
    
        // calculate angles like Ming-Yan
        TLorentzVector l_minus, l_plus; 
        if (leptonOneFlavor > 0) {
            l_minus = leptonOneP4;
            l_plus = leptonTwoP4;
        }
        else {
            l_minus = leptonTwoP4;
            l_plus = leptonOneP4;
        }
        
        TVector3 llgFrame = -1*llgP4.BoostVector();
        dileptonP4.Boost(llgFrame);
        l_minus.Boost(llgFrame);
        l_minus.Boost(dileptonP4.BoostVector());
        zgLittleTheta = cos(dileptonP4.Angle(l_minus.Vect()));
        zgBigTheta = cos(dileptonP4.Angle(llgP4.Vect()));
        
        // the way MY does it (I think wrong)
        TLorentzVector lep0 = leptonOneP4;
        lep0.Boost(llgFrame);
        lep0.Boost(dileptonP4.BoostVector());
        zgLittleThetaMY = cos(dileptonP4.Angle(lep0.Vect()));
        
        TVector3 ppAxis(0, 0, 1);
        TVector3 zAxis = dileptonP4.Vect().Unit();
        TVector3 yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rot;
        rot = rot.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4.Transform(rot);
        l_minus.Transform(rot);
        zgPhi = l_minus.Phi();
            
        if (!isData) {

            elIDWeightOne = weights->GetHZZElectronRecoIdEff(*electrons[electronOneIndex]); 
            elIDWeightTwo = weights->GetHZZElectronRecoIdEff(*electrons[electronTwoIndex]); 
            eventWeight *= elIDWeightOne;
            eventWeight *= elIDWeightTwo;
            
            float sf11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[electronOneIndex]);
            float sf22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[electronTwoIndex]);

            elTrigWeightOne = sf11;
            elTrigWeightTwo = sf22;

            triggerWeight = elTrigWeightOne*elTrigWeightTwo;
            eventWeight *= triggerWeight;
           
            photonIDWeight = weights->GetPhotonMVAIdEff(*photons[photonIndex]); 
            eventWeight *= photonIDWeight;
            
            if (sync_print_precut) {
                cout << "run,lumi,evt,puwei,totSF,trg0,trg1,id0,id1,gammaID,pt0,pt1,calibpt0,calibpt1,eta0,eta1,sceta0,sceta1" << endl;
                cout << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << ", "
                     << puWeight << ", " << eventWeight << ", " << elTrigWeightOne << ", "
                     << elTrigWeightTwo << ", " << elIDWeightOne << ", " << elIDWeightTwo << ", "
                     << photonIDWeight << ", " << electrons[electronOneIndex]->pt << ", " << electrons[electronTwoIndex]->pt << ", "
                     << electrons[electronOneIndex]->calibPt << ", " << electrons[electronTwoIndex]->calibPt << ", "
                     << electrons[electronOneIndex]->eta << ", " << electrons[electronTwoIndex]->eta << ", "
                     << electrons[electronOneIndex]->scEta << ", " << electrons[electronTwoIndex]->scEta << endl;
            }
        
        }
	//cout << "----- END --------\n" ;
    } // end elelg selection

    else if (params->selection == "tautaug") {
        if (muons.size() != 1) // avoid multi-muon events (DY contamination)
            return kTRUE;
        hTotalEvents->Fill(5);
        if (taus.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (muons[0]->pt <= 25.) 
            return kTRUE;
        hTotalEvents->Fill(8);

        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);

        unsigned int tau_index = 0;
        for (unsigned int i = 0; i < taus.size(); ++i) {
            if (taus[i]->q != muons[0]->q && taus[i]->pt > 20.) {
                tau_index = i;
                break;
            }
        }

        if (taus[tau_index]->pt <= 20.)
            return kTRUE;
        hTotalEvents->Fill(9);

        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(taus[tau_index]->pt, taus[tau_index]->eta, taus[tau_index]->phi, taus[tau_index]->m);

        if (photons[0]->calibPt <= 15.)
            return kTRUE;
        hTotalEvents->Fill(10);

        // muon transverse mass cut to reject ttbar and w+jets
     
        // event passed the selection; fill output variables
        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
        photonOnePt  = photonOneP4.Pt();
        photonOneEta = photonOneP4.Eta();
        photonOnePhi = photonOneP4.Phi();
        photonOneMVA = photons[0]->mvaFall17V2;
        passElectronVeto = photons[0]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[0]);
        else 
            photonOneR9 = photons[0]->r9;

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

        tauDecayMode    = taus[tau_index]->decaymode;
        tauMVA          = taus[tau_index]->rawIsoMVA3newDMwLT;

        if (muonP4.Pt() > tauP4.Pt()) {
            leptonOnePt     = muonP4.Pt();
            leptonOneEta    = muonP4.Eta();
            leptonOnePhi    = muonP4.Phi();
            leptonOneIso    = GetMuonIsolation(muons[0]);
            leptonOneFlavor = muons[0]->q*13;
            leptonOneDZ     = muons[0]->dz;
            leptonOneD0     = muons[0]->d0;

            leptonTwoPt     = tauP4.Pt();
            leptonTwoEta    = tauP4.Eta();
            leptonTwoPhi    = tauP4.Phi();
            leptonTwoIso    = 0.;
            leptonTwoFlavor = 15*taus[tau_index]->q;
            leptonTwoDZ     = taus[tau_index]->dzLeadChHad;
            leptonTwoD0     = taus[tau_index]->d0LeadChHad;
        }
        else {
            leptonOnePt     = tauP4.Pt();
            leptonOneEta    = tauP4.Eta();
            leptonOnePhi    = tauP4.Phi();
            leptonOneIso    = 0.;
            leptonOneFlavor = 15*taus[tau_index]->q;
            leptonOneDZ     = taus[tau_index]->dzLeadChHad;
            leptonOneD0     = taus[tau_index]->d0LeadChHad;

            leptonTwoPt     = muonP4.Pt();
            leptonTwoEta    = muonP4.Eta();
            leptonTwoPhi    = muonP4.Phi();
            leptonTwoIso    = GetMuonIsolation(muons[0]);
            leptonTwoFlavor = muons[0]->q*13;
            leptonTwoDZ     = muons[0]->dz;
            leptonTwoD0     = muons[0]->d0;
        }
    
        isDijetTag = false; // to ensure proper jet filling

        // MC event weights
        if (!isData) {

            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            //eventWeight *= weights->GetMuonISOEff(muonP4);
            eventWeight *= weights->GetPhotonMVAIdEff(*photons[0]);
            eventWeight *= 0.95; // flat tau id scale factor

            float eff_data = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4).first; 
            float eff_mc = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4).second; 
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;
            cout << "eventWeight = " << eventWeight << endl;
            
        }

    } // end tautaug selection


    ///////////////////
    // Fill jet info //
    ///////////////////
    

    if (!isDijetTag) {
        if (jets.size() > 0) {
            jetOnePt   = jets[0]->pt;
            jetOneEta  = jets[0]->eta;
            jetOnePhi  = jets[0]->phi;
            jetOneM    = jets[0]->mass;
            jetOneTag  = jets[0]->csv;
        } else {
            jetOnePt  = 0.;
            jetOneEta = 0.;
            jetOnePhi = 0.;
            jetOneM   = 0.;
            jetOneTag = 0.;
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
    Sgen = 0;
    SgenAccep = 0;

    if (sync_print_precut) {
        cout << "event should have been filled" << endl;
        }
    return kTRUE;
}

void zgAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zgAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zgAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<zgAnalyzer> selector(new zgAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zgAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float zgAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float zgAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
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

float zgAnalyzer::GetGenIsolation(const TGenParticle* pho)
{

	float genIso = 0;
	int isStable = 1;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
		TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

		float dR = sqrt(pow((pho->eta - particle->eta),2) + pow((pho->phi - particle->phi),2));
		if(particle->status == isStable && dR < 0.4)
			genIso += particle->pt;
	}
	return genIso;
}

bool zgAnalyzer::SignalRegionPass(const baconhep::TPhoton *photon){
        bool signalPass = false;
        if(photon->scEta <= 1.48){
                if(photon->chHadIso <= 2.0)
                        signalPass = true;
        }
        else {
                if(photon->chHadIso <= 1.5)
                        signalPass = true;
        }
        return signalPass;
}


float zgAnalyzer::GetWorstChIsolation(const TPhoton* pho)
{

        unsigned int nMuon = fMuonArr    ->GetEntries();
        unsigned int nElec = fElectronArr->GetEntries();
        unsigned int nTau  = fTauArr     ->GetEntries();
        unsigned int nJets = fAK4CHSArr  ->GetEntries();


        float worstChIso = 0;
        int isStable = 1;
        for(unsigned int i ; i < nMuon ; i++){
                TMuon *muon = (TMuon * ) fMuonArr->At(i);
                float dR = sqrt(pow((pho->scEta - muon->eta),2) + pow((pho->phi - muon->phi),2));
                if( dR < 0.4 && dR > 0.1)
                        worstChIso += muon->pt;
        }


        for(unsigned int i ; i < nElec ; i++){
                TElectron *electron = (TElectron *) fElectronArr->At(i);
                float dR = sqrt(pow((pho->scEta - electron->scEta),2) + pow((pho->phi - electron->phi),2));
                if( dR < 0.4 && dR > 0.1)
                        worstChIso += electron->pt;
        }

        for(unsigned int i ; i < nTau ; i++){
                TTau *tau = (TTau *) fTauArr->At(i);
                float dR = sqrt(pow((pho->scEta - tau->eta),2) + pow((pho->phi - tau->phi),2));
                if( dR < 0.4 && dR > 0.1)
                        worstChIso += tau->pt;
        }


        for(unsigned int i ; i < nJets ; i++){
                TJet *jets = (TJet *) fAK4CHSArr->At(i);
                float dR = sqrt(pow((pho->scEta - jets->eta),2) + pow((pho->phi - jets->phi),2));
                if( dR < 0.4 && dR > 0.1)
                        worstChIso += jets->pt;
        }


        return worstChIso;
}


 void zgAnalyzer::find_optimized(double* p, double &e1, double& e2)
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
