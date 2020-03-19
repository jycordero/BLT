#include "zAnalyzerTrig.h"
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

zAnalyzerTrig::zAnalyzerTrig() : BLTSelector()
{

}

zAnalyzerTrig::~zAnalyzerTrig()
{

}

void zAnalyzerTrig::Begin(TTree *tree)
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
    if(params->period == "2016" || params->period == "2016ReReco" || params->period == "2016Legacy"){
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
    else if(params->period == "2017" || params->period == "2017ReReco" || params->period == "2017Legacy"){
        if ( params->selection == "mumug") {
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
        }
	else if( params->selection == "mumu" ){
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_v*");
		triggerNames.push_back("HLT_Mu8_TrkIsoVVL_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
		triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
		
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
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if( params->period == "2016Legacy"){
	jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt";
    	muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/ReReco2016/roccor.Run2.v3/RoccoR2016.txt");
    }	
    else if(params->period == "2016ReReco"){
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"; 
    	muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/ReReco2016/roccor.Run2.v3/RoccoR2017.txt");
    }
    else if(params->period == "2017ReReco"){
	//jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt";
	jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
	//jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/" + PeriodFolder + "/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt";
    	muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/ReReco2016/roccor.Run2.v3/RoccoR2017.txt");
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

   
    outTree->Branch("genWeight"   , &genWeight);
    outTree->Branch("eventWeight" , &eventWeight);
    outTree->Branch("puWeight"    , &puWeight);

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

    // Gen objects
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

    outTree->Branch("ProbeTrigPass"         , &ProbeTrigPass);
    outTree->Branch("ProbeGlobalTrigMu17Pass"     , &ProbeGlobalTrigMu17Pass);
    outTree->Branch("ProbeGlobalTrigMu8Pass"     , &ProbeGlobalTrigMu8Pass);
    outTree->Branch("ProbeTrigMu17Pass"     , &ProbeTrigMu17Pass);
    outTree->Branch("ProbeTrigMu8Pass"      , &ProbeTrigMu8Pass);
    outTree->Branch("ProbeTrigMu17LegPass"  , &ProbeTrigMu17LegPass);
    outTree->Branch("ProbeTrigMu8LegPass"   , &ProbeTrigMu8LegPass);
    outTree->Branch("ProbeIDPass"           , &ProbeIDPass);
    outTree->Branch("ProbeIDPass"           , &ProbeIDPass);
    outTree->Branch("ProbeISOPass"          , &ProbeISOPass);
    outTree->Branch("ProbeWorstPass"        , &ProbeWorstPass);
    outTree->Branch("ProbeSigPass"          , &ProbeSigPass);
    // photons
    outTree->Branch("photonOnePt"     , &photonOnePt);
    outTree->Branch("photonOneEta"    , &photonOneEta);
    outTree->Branch("photonOnePhi"    , &photonOnePhi);
    outTree->Branch("photonOneR9"     , &photonOneR9);
    outTree->Branch("photonOneMVA"    , &photonOneMVA);
    outTree->Branch("photonOneERes"   , &photonOneERes);
    outTree->Branch("photonOneLeptonDR", &photonOneLeptonDR);
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



    // object counters
    outTree->Branch("nMuons"      , &nMuons);
    outTree->Branch("nElectrons"  , &nElectrons);
    outTree->Branch("nTaus"       , &nTaus);
    outTree->Branch("nPhotons"    , &nPhotons);
    outTree->Branch("nJets"       , &nJets);
    outTree->Branch("nFwdJets"    , &nFwdJets);
    outTree->Branch("nCentralJets", &nCentralJets);
    outTree->Branch("nBJets"      , &nBJets);
    
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
    

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    string outHistNameGen = params->get_output_treename("TotalEventsGen");
    hTotalEventsGen = new TH1D(outHistNameGen.c_str(),"TotalEventsGen",5,0.5,5.5);

    int ptNBins  = 6;
    int etaNBins = 10;

    double ptBins [ptNBins + 1] = {0,20,40,50,90,150,500};
    //double etaBins[etaNBins + 1] = {-2.5,-2,-1.566,-1.4442,-1.0,0,1.0,1.4442,1.566,2,2.5};
    double etaBins[etaNBins + 1] = {-2.5,-1.566,-1.4442,0,1.4442,1.566,2.5};

    PhotonProbe          = new TH2F("EGammaProbe"          , "EGammaProbe"          , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeIDPass    = new TH2F("EGammaProbeIDPass"    , "EGammaProbeIDPass"    , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeIDFail    = new TH2F("EGammaProbeIDFail"    , "EGammaProbeIDFail"    , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeISOPass   = new TH2F("EGammaProbeISOPass"   , "EGammaProbeISOPass"   , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeISOFail   = new TH2F("EGammaProbeISOFail"   , "EGammaProbeISOFail"   , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeSigPass   = new TH2F("EGammaProbeSigPass"   , "EGammaProbeSigPass"   , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeSigFail   = new TH2F("EGammaProbeSigFail"   , "EGammaProbeSigFail"   , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeWorstPass = new TH2F("EGammaProbeWorstPass" , "EGammaProbeWorstPass" , ptNBins, ptBins, etaNBins, etaBins);
    PhotonProbeWorstFail = new TH2F("EGammaProbeWorstFail" , "EGammaProbeWorstFail" , ptNBins, ptBins, etaNBins, etaBins);

    MuonProbe          = new TH2F("Probe"          , "Probe"          , ptNBins, ptBins, etaNBins, etaBins);
    MuonProbeIDPass    = new TH2F("ProbeIDPass"    , "ProbeIDPass"    , ptNBins, ptBins, etaNBins, etaBins);
    MuonProbeIDFail    = new TH2F("ProbeIDFail"    , "ProbeIDFail"    , ptNBins, ptBins, etaNBins, etaBins);
    MuonProbeIDISOPass = new TH2F("ProbeIDISOPass" , "ProbeIDISOPass" , ptNBins, ptBins, etaNBins, etaBins);
    MuonProbeIDISOFail = new TH2F("ProbeIDISOFail" , "ProbeIDISOFail" , ptNBins, ptBins, etaNBins, etaBins);

    ReportPostBegin();
}

Bool_t zAnalyzerTrig::Process(Long64_t entry)
{
    //bool debug = true;
    bool debug = false;
    //bool debugSelect = true;
    bool debugSelect = false;

    //bool useGlobalTrigger = true;
    bool useGlobalTrigger = false;
   
    //bool debugSF = true;
    bool debugSF = false;


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
                
            if (fabs(particle->pdgId) == 22) 
                genPhotons.push_back(particle);    

            }
            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

    } else {
        nPartons = 0;
    }
      

    if (debug)
	cout<< "Before Lumi mask passed\n";
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

    if (!passTrigger && useGlobalTrigger)// && isData)
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
    vector<TMuon*> muons_trig;

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

	if(
		   muonP4.Pt() > 10 
		&& fabs(muonP4.Eta()) < 2.4
	)
		muons_trig.push_back(muon);
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    sort(muons_trig.begin(), muons_trig.end(), sort_by_higher_pt<TMuon>);

    if(debug)
	cout << "----Muon collection";
    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i < fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        //assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->calibPt, electron->scEta, electron->phi, ELE_MASS);
        //if(debugSelect){
	//	cout    << "Electron info " << endl 
	//		<< "-ID:: " << particleSelector->PassElectronID(electron, cuts->tightElID) << endl 
	//		<< "-ISO:: " << particleSelector->PassElectronIso(electron,cuts->tightElIso) << endl; 
	//}
 
        if(sync_print_precut) {
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
                && fabs(electron->scEta) < 2.1
                //&& particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)
                && (fabs(electron->scEta) <= 1.4442 || fabs(electron->scEta) >= 1.566)
                //&& particleSelector->PassElectronID(electron, cuts->tightElID)
		//&& particleSelector->PassElectronIso(electron,cuts->tightElIso)
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
        photonP4.SetPtEtaPhiM(photon->calibPt, photon->scEta, photon->phi, 0.);
    	//cout << " ---------PhotonP4 PT :: " << photonP4.Pt() << endl;
    	if(debug){
		cout << " ---------Photon PT   :: " << photon->pt << " cal: " << photon->calibPt << endl;
		cout << "----------Photon Eta  :: " << photon->eta << " scEta: " << photon->scEta << endl;
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
                photon->calibPt > 15
                //photon->calibPt > 10
                && fabs(photon->scEta) < 2.5 
                && (fabs(photon->scEta) <= 1.4442 || fabs(photon->scEta) >= 1.566)
                //&& particleSelector->PassPhotonID(photon, cuts->preSelPhID)
                //&& particleSelector->PassPhotonIso(photon, cuts->preSelPhIso,EAPho)	
                //&& photon->passElectronVeto
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

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nTaus      = taus.size();
    nPhotons   = photons.size();

    if(debug)
	cout << "---- Start Selection" << endl;
    
    if (params->selection == "ee") {
	    /////////////////////////////////////////////////////////
	    if (nElectrons < 1)
	    	return kTRUE;
	    if (nPhotons < 1)
	    	return kTRUE;

	    float zMassDiff = 999;
	
	    int TagIndex   = 0;
	    int ProbeIndex = 0;


	    bool hasValidPair = false;

	    ProbeIDPass    = false;
	    ProbeISOPass   = false;
	    ProbeWorstPass = false;
	    ProbeSigPass   = false;

	    bool ElGmOverlap = false;
	    float dR;
	    //////////////////////////////////////////////////////////
	    if(debugSelect)
		cout << "---- Start Electron Tag" <<endl;
	    for (unsigned int i = 0; i < electrons.size(); ++i) {
	    	TElectron* electron = electrons[i];
		assert(electron);

		TLorentzVector tempElec;
		tempElec.SetPtEtaPhiM(electron->calibPt, electron->scEta, electron->phi, ELE_MASS);
		
		if( 	   tempElec.Pt() > 30  
			&& fabs(tempElec.Eta()) < 2.1 
			&& (fabs(tempElec.Eta()) <= 1.4442 || fabs(tempElec.Eta()) >= 1.566) 
			&& trigger->pass("HLT_Ele27_WPTight_Gsf_v*", fInfo->triggerBits)
			&& particleSelector->PassElectronID(electron, cuts->tightElID)
			//&& particleSelector->PassElectronIso(electron,cuts->tightElIso)
		){
	    		// Electron has been TAGGED
			if(debugSelect)
				cout << "---- Start Electron Probe" <<endl;
			for (unsigned int j = 0; j < photons.size(); j++) {
				TPhoton* photon = photons[j];
				assert(photon);
				
				TLorentzVector tempPhoton;
				tempPhoton.SetPtEtaPhiM(photon->calibPt, photon->scEta, photon->phi, ELE_MASS);
				
				ElGmOverlap = false;
				//cout << "----------------------------\n";
				for (unsigned int k = 0; k < electrons.size(); ++k) {
					if(k != i){
						TElectron* elec = electrons[k];
						assert(elec);

						dR = sqrt(pow(elec->scEta - photon->scEta,2) + pow( elec->phi - photon->phi,2));
						//cout << "    dR = " << dR << endl;
						if(dR < 0.1 ){
							ElGmOverlap = true;
							photonOneLeptonDR = dR;
						}
					}
				}
				//cout << " Elec Gm Overlap Flag :: "<< ElGmOverlap << endl;	
				TLorentzVector tempDilep;	
				tempDilep = tempElec + tempPhoton;
				if(		
					   tempPhoton.Pt() > 15 
					&& fabs(tempPhoton.Eta()) < 2.5  
					&& (fabs(tempPhoton.Eta()) <= 1.4442 || fabs(tempPhoton.Eta()) >= 1.566) 
					//SignalRegionPass(photon) && // need access to variables in the TPhoton object
					&& ElGmOverlap
				){
					if (tempDilep.M() > 60.0 && tempDilep.M() < 120) {
						if (hasValidPair) {
							if (fabs( tempDilep.M()- ZMASS) < zMassDiff) {
								zMassDiff = fabs(tempDilep.M() - ZMASS);

								TagIndex   = i;
								ProbeIndex = j;
							}
						}
						else {
							zMassDiff = fabs(tempDilep.M() - ZMASS);

							TagIndex   = i;
							ProbeIndex = j;
							hasValidPair = true;
						}
					}
				}
			}
		}
	    }
	    if(!hasValidPair)
	    	return kTRUE;


	    TLorentzVector Dilep, Photon, Electron;
	    Electron.SetPtEtaPhiM(electrons[TagIndex]->calibPt, electrons[TagIndex]->scEta, electrons[TagIndex]->phi, ELE_MASS);
	    Photon  .SetPtEtaPhiM(photons[ProbeIndex]->calibPt, photons[ProbeIndex]->scEta, photons[ProbeIndex]->phi, ELE_MASS);

	    Dilep = Electron + Photon;

	    PhotonProbe->Fill(Photon.Pt(),Photon.Eta(),1);



	    if(    particleSelector->PassPhotonID (photons[ProbeIndex], cuts->preSelPhID) ){
		    ProbeIDPass = true;
		    PhotonProbeIDPass->Fill(Photon.Pt(),Photon.Eta(),1);
	    }
	    else{
		    ProbeIDPass = false;
		    PhotonProbeIDFail->Fill(Photon.Pt(),Photon.Eta(),1);
	    }

	    if(    particleSelector->PassPhotonIso(photons[ProbeIndex], cuts->preSelPhIso,EAPho) ){
		    ProbeISOPass = true;
		    PhotonProbeISOPass->Fill(Photon.Pt(),Photon.Eta(),1);
	    }
	    else{
		    ProbeISOPass = false;
		    PhotonProbeISOFail->Fill(Photon.Pt(),Photon.Eta(),1);
	    }

	    if( SignalRegionPass(photons[ProbeIndex]) ){
		    ProbeSigPass = true;
		    PhotonProbeSigPass->Fill(Photon.Pt(),Photon.Eta(),1);
	    }
	    else{
		    ProbeSigPass = false;
		    PhotonProbeSigFail->Fill(Photon.Pt(),Photon.Eta(),1);
	    }

	    if( GetWorstChIsolation(photons[ProbeIndex]) < 15 ){
		    ProbeWorstPass = true;
		    PhotonProbeWorstPass->Fill(Photon.Pt(),Photon.Eta(),1);
	    }
	    else{
		    ProbeWorstPass = false;
		    PhotonProbeWorstPass->Fill(Photon.Pt(),Photon.Eta(),1);
	    }

	    leptonOnePt  = Electron.Pt();
	    leptonOneEta = Electron.Eta();
	    leptonOnePhi = Electron.Phi();

	    leptonTwoPt  = Photon.Pt();
	    leptonTwoEta = Photon.Eta();
	    leptonTwoPhi = Photon.Phi();	

	    dileptonPt   = Dilep.Pt();
	    dileptonEta  = Dilep.Eta();
	    dileptonPhi  = Dilep.Phi();
	    dileptonM    = Dilep.M();

	    if(debugSelect)
		cout << "---- Start Fill" <<endl;

    }
    else if (params->selection == "mumu") {
	    /////////////////////////////////////////////////////////
	    if (nMuons < 2)
	    	return kTRUE;
	

	    float zMassDiff = 999;
	
	    int TagIndex   = 0;
	    int ProbeIndex = 0;

	    bool hasValidPair = false;

	    ProbeTrigPass        = false; 
	    ProbeTrigMu17Pass    = false; 
	    ProbeTrigMu8Pass     = false; 
	    ProbeTrigMu17LegPass = false; 
	    ProbeTrigMu8LegPass  = false; 
	    ProbeIDPass          = false;
	    ProbeISOPass         = false;

	    float dR;
	    //debugSelect = true;
	    if(debugSelect)
		cout << "---- Start Muon Tag" <<endl;
	    for (unsigned int i = 0; i < muons_trig.size(); ++i) {
	    	TMuon* muon = muons_trig[i];
		assert(muon);

		TLorentzVector tempMuon;
		tempMuon.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
		
		if( 	   tempMuon.Pt() > 30  
			&& fabs(tempMuon.Eta()) < 2.1
			&& particleSelector->PassMuonID(muon, cuts->tightMuID)
			//&& particleSelector->PassElectronIso(electron,cuts->tightElIso)
		){
	    		// Electron has been TAGGED
			if(debugSelect)
				cout << "---- Start Muon Probe" <<endl;
			for (unsigned int j = 0; j < muons_trig.size(); j++) {
				if(j != i){
					TMuon* muonP = muons_trig[j];
					assert(muonP);
					
					TLorentzVector tempMuonProbe;
					tempMuonProbe.SetPtEtaPhiM(muonP->pt, muonP->eta, muonP->phi, MUON_MASS);
					
					//cout << " Elec Gm Overlap Flag :: "<< ElGmOverlap << endl;	
					TLorentzVector tempDilep;	
					tempDilep = tempMuon + tempMuonProbe;
					if(		
						   tempMuonProbe.Pt() > 15 
						&& fabs(tempMuonProbe.Eta()) < 2.4   
					){
						if (tempDilep.M() > 60.0 && tempDilep.M() < 120) {
							if (hasValidPair) {
								if (fabs( tempDilep.M()- ZMASS) < zMassDiff) {
									zMassDiff = fabs(tempDilep.M() - ZMASS);

									TagIndex   = i;
									ProbeIndex = j;
								}
							}
							else {
								zMassDiff = fabs(tempDilep.M() - ZMASS);

								TagIndex   = i;
								ProbeIndex = j;
								hasValidPair = true;
							}
						}
					}
				}
			}
		}
	    }
	    if(!hasValidPair)
	    	return kTRUE;
	    //cout << "Pair obtained!!" << endl;

	    TLorentzVector Dilep, Muon, MuonProbe;
	    Muon     .SetPtEtaPhiM(muons_trig[  TagIndex]->pt, muons_trig[  TagIndex]->eta, muons_trig[  TagIndex]->phi, MUON_MASS);
	    MuonProbe.SetPtEtaPhiM(muons_trig[ProbeIndex]->pt, muons_trig[ProbeIndex]->eta, muons_trig[ProbeIndex]->phi, MUON_MASS);

	    Dilep = Muon + MuonProbe;

	    PhotonProbe->Fill(MuonProbe.Pt(),MuonProbe.Eta(),1);

	    //string trigNAME = "HLT_Mu17_TrkIsoVVL_v*";
	    //if(trigger->pass( trigNAME, muons_trig[ProbeIndex]->hltMatchBits))
	    //    cout << " Global Trigger " << trigNAME << " PASS" << endl;
	    //else	
	    //	cout << " Global Trigger " << trigNAME << " FAILED" << endl;

	    string trigNAME;
	    //cout << " Global Trigger " << trigNAME << " PASS" << endl;
	    cout << " HLT Trig GLOBAL:: " <<  fInfo->triggerBits << endl;

	    trigNAME = "HLT_Mu17_TrkIsoVVL_v*";
	    if(trigger->pass( trigNAME, fInfo->triggerBits) ){
	        //cout << " Global Trigger " << trigNAME << " PASS" << endl;
		ProbeGlobalTrigMu17Pass = true;
	    }
	    else 
		ProbeGlobalTrigMu17Pass = false;
	    trigNAME = "HLT_Mu8_TrkIsoVVL_v*";
	    if(trigger->pass( trigNAME, fInfo->triggerBits) ){
	        //cout << " Global Trigger " << trigNAME << " PASS" << endl;
		ProbeGlobalTrigMu8Pass = true;
	    }
	    else
		ProbeGlobalTrigMu8Pass = false;
	   

	
	    ///////////////////////////////////////////////////
	    bool DoubleMuTrigLeg1 = false;
	    bool DoubleMuTrigLeg2 = false;
	    for(unsigned int iT = 0; iT < triggerNames.size()  ; iT++){
	    	DoubleMuTrigLeg1 |= trigger->passObj(triggerNames.at(iT), 1, muons_trig[ProbeIndex]->hltMatchBits);
	    	DoubleMuTrigLeg2 |= trigger->passObj(triggerNames.at(iT), 2, muons_trig[ProbeIndex]->hltMatchBits);
	    }



	    if( DoubleMuTrigLeg1 ){
		ProbeTrigMu17LegPass = true;
		cout << " Leg Mu17 Pass Obj" << endl;
	    }
	    else{
		ProbeTrigMu17LegPass = false;
	    }
	    if( DoubleMuTrigLeg2 ){
		ProbeTrigMu8LegPass = true;
		cout << " Leg Mu8 Pass Obj" << endl;
	    }
	    else{
		ProbeTrigMu8LegPass = false;
	    }
	    ///////////////////////////////////////////////////

	    //DoubleMuTrigLeg1 = false;
	    //DoubleMuTrigLeg2 = false;
	    //for(unsigned int iT = 0; iT < triggerNames.size()  ; iT++){
	    //	DoubleMuTrigLeg1 |= trigger->passObj(triggerNames.at(iT), 1, fInfo->triggerBits);
	    //	DoubleMuTrigLeg2 |= trigger->passObj(triggerNames.at(iT), 2, fInfo->triggerBits);
	    //}



	    //if( DoubleMuTrigLeg1 ){
	    //    ProbeTrigMu17LegPass = true;
	    //    cout << " Leg Mu17 Pass Obj" << endl;
	    //}
	    //else{
	    //    ProbeTrigMu17LegPass = false;
	    //}
	    //if( DoubleMuTrigLeg2 ){
	    //    ProbeTrigMu8LegPass = true;
	    //    cout << " Leg Mu8 Pass Obj" << endl;
	    //}
	    //else{
	    //    ProbeTrigMu8LegPass = false;
	    //}
	    //     
	    ///////////////////////////////////////////////////
	    
	    //cout << " HLT Trig:: " <<  muons_trig[ProbeIndex]->hltMatchBits << endl;

	    if( trigger->passObj( "HLT_Mu17_TrkIsoVVL_v*" , 1, muons_trig[ProbeIndex]->hltMatchBits) ){
	    	ProbeTrigMu17Pass = true;
	    }
	    else{
	    	ProbeTrigMu17Pass = false;
	    }

	    if( trigger->passObj( "HLT_Mu8_TrkIsoVVL_v*" , 1, muons_trig[ProbeIndex]->hltMatchBits) ){
	    	ProbeTrigMu8Pass = true;
	    }
	    else{
	    	ProbeTrigMu8Pass = false;
	    }
	
	    leptonOnePt  = Muon.Pt();
	    leptonOneEta = Muon.Eta();
	    leptonOnePhi = Muon.Phi();

	    leptonTwoPt  = MuonProbe.Pt();
	    leptonTwoEta = MuonProbe.Eta();
	    leptonTwoPhi = MuonProbe.Phi();	

	    dileptonPt   = Dilep.Pt();
	    dileptonEta  = Dilep.Eta();
	    dileptonPhi  = Dilep.Phi();
	    dileptonM    = Dilep.M();

    }

        

    outTree->Fill();
    this->passedEvents++;

    if (sync_print_precut) {
        cout << "event should have been filled" << endl;
        }
    return kTRUE;
}

void zAnalyzerTrig::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zAnalyzerTrig::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zAnalyzerTrig::ReportPostTerminate()
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
    std::unique_ptr<zAnalyzerTrig> selector(new zAnalyzerTrig());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zAnalyzerTrig::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float zAnalyzerTrig::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float zAnalyzerTrig::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
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

float zAnalyzerTrig::GetGenIsolation(const TGenParticle* pho)
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

float zAnalyzerTrig::GetWorstChIsolation(const TPhoton* pho)
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

void zAnalyzerTrig::find_optimized(double* p, double &e1, double& e2)
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
bool zAnalyzerTrig::SignalRegionPass(const baconhep::TPhoton *photon){
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