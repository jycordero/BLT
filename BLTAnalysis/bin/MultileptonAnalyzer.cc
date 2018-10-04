#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

TLorentzVector GetP4Sum(vector<TLorentzVector> p4)
{
    TLorentzVector sum(0, 0, 0, 0);
    for (auto p4_ = p4.begin(); p4_ != p4.end(); ++p4_)
        sum += *p4_;
    return sum;
}


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
{
    rng = new TRandom3();

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));


    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);


    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));


    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Trigger_requirements
    // https://indico.cern.ch/event/682891/contributions/2810364/attachments/1570825/2477991/20171206_CMSWeek_MuonHLTReport_KPLee_v3_1.pdf
    if (params->selection == "double")
    {
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        electronTriggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
    }
    else
        cout << "WARNING: No triggers selected!" << endl;


    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
//  string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/";
    size_t isDoubleMuon = params->datasetgroup.find("muon_2017");
    if (isDoubleMuon != string::npos)   // Remove period 2017B from muon data
        jsonFileName += "DoubleMuon_299368-306462_JSON.txt";
    else
        jsonFileName += "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
    cout << jsonFileName << endl;
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/RoccoR2017v1.txt");


    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");



    //--- BRANCHES ---//

    // Event
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("nPV", &nPV);

    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("PUWeight", &PUWeight);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    outTree->Branch("evtMuonTriggered", &evtMuonTriggered);
    outTree->Branch("evtElectronTriggered", &evtElectronTriggered);


    // Counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nLeptons", &nLeptons);

    outTree->Branch("nIsoMVAElectrons", &nIsoMVAElectrons);
    outTree->Branch("nNoIsoMVAElectrons", &nNoIsoMVAElectrons);
    outTree->Branch("nLooseMuons", &nLooseMuons);

    outTree->Branch("nHZZMuons", &nHZZMuons);
    outTree->Branch("nHZZElectrons", &nHZZElectrons);
    outTree->Branch("nHZZLeptons", &nHZZLeptons);

    outTree->Branch("nGenMuons", &nGenMuons);
    outTree->Branch("nGenElectrons", &nGenElectrons);
    outTree->Branch("nGenLeptons", &nGenLeptons);


    // Muons
    outTree->Branch("muonP4", &muonsP4, 32000, 1);
    outTree->Branch("muonQ", &muonsQ);
    outTree->Branch("muonFiredLeg1", &muonFiredLeg1);
    outTree->Branch("muonFiredLeg2", &muonFiredLeg2);

    outTree->Branch("muonIsGhost", &muonIsGhost);
    outTree->Branch("muonIsLoose", &muonIsLoose);
    outTree->Branch("muonIsHZZ", &muonIsHZZ);

    outTree->Branch("muonEnergySF", &muonEnergySF);
    outTree->Branch("muonHZZIDSF", &muonHZZIDSF);
    outTree->Branch("muonTrigEffLeg1Data", &muonTrigEffLeg1Data);
    outTree->Branch("muonTrigEffLeg1MC", &muonTrigEffLeg1MC);
    outTree->Branch("muonTrigEffLeg2Data", &muonTrigEffLeg2Data);
    outTree->Branch("muonTrigEffLeg2MC", &muonTrigEffLeg2MC);

    outTree->Branch("muonCombIso", &muonCombIso);
    outTree->Branch("muonTrkIso", &muonsTrkIso);
    outTree->Branch("muonD0", &muonD0);
    outTree->Branch("muonDz", &muonDz);
    outTree->Branch("muonSIP3d", &muonSIP3d);
    outTree->Branch("muonPtErr", &muonPtErr);
    outTree->Branch("muonNMatchStn", &muonNMatchStn);
    outTree->Branch("muonNPixHits", &muonNPixHits);
    outTree->Branch("muonNTkLayers", &muonNTkLayers);
    outTree->Branch("muonIsPF", &muonIsPF);
    outTree->Branch("muonIsGlobal", &muonIsGlobal);
    outTree->Branch("muonIsTracker", &muonIsTracker);
    outTree->Branch("muonBestTrackType", &muonBestTrackType);


    // Electrons
    outTree->Branch("electronP4", &electronsP4, 32000, 1);
    outTree->Branch("electronQ", &electronsQ);
    outTree->Branch("electronFiredLeg1", &electronFiredLeg1);
    outTree->Branch("electronFiredLeg2", &electronFiredLeg2);

    outTree->Branch("electronIsGhost", &electronIsGhost);
    outTree->Branch("electronPassIsoMVA", &electronPassIsoMVA);
    outTree->Branch("electronPassNoIsoMVA", &electronPassNoIsoMVA);
    outTree->Branch("electronIsHZZ", &electronIsHZZ);

    outTree->Branch("electronEnergySF", &electronEnergySF);
    outTree->Branch("electronHZZIDSF", &electronHZZIDSF);
    outTree->Branch("electronTrigEffLeg1Data", &electronTrigEffLeg1Data);
    outTree->Branch("electronTrigEffLeg1MC", &electronTrigEffLeg1MC);
    outTree->Branch("electronTrigEffLeg2Data", &electronTrigEffLeg2Data);
    outTree->Branch("electronTrigEffLeg2MC", &electronTrigEffLeg2MC);

    outTree->Branch("electronIsoMVA", &electronIsoMVA);
    outTree->Branch("electronNoIsoMVA", &electronNoIsoMVA);
    outTree->Branch("electronCombIso", &electronCombIso);
    outTree->Branch("electronTrkIso", &electronsTrkIso);
    outTree->Branch("electronD0", &electronD0);
    outTree->Branch("electronDz", &electronDz);
    outTree->Branch("electronSIP3d", &electronSIP3d);
    outTree->Branch("electronScEta", &electronScEta);
    outTree->Branch("electronNMissHits", &electronNMissHits);
    outTree->Branch("electronIsGap", &electronIsGap);


    // Gen particles
    outTree->Branch("genMuonP4", &genMuonsP4, 32000, 1);
    outTree->Branch("genMuonQ", &genMuonsQ);
    outTree->Branch("genMuonStatus", &genMuonStatus);

    outTree->Branch("genElectronP4", &genElectronsP4, 32000, 1);
    outTree->Branch("genElectronQ", &genElectronsQ);
    outTree->Branch("genElectronStatus", &genElectronStatus);


    // Histograms
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);
    outHistName = params->get_output_treename("PhaseSpaceEvents");
    hPhaseSpaceEvents = new TH1D(outHistName.c_str(), "PhaseSpaceEvents", 10, 0.5, 10.5);
    outHistName = params->get_output_treename("FiducialEvents");
    hFiducialEvents = new TH1D(outHistName.c_str(), "FiducialEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{


    //--- CLEAR CONTAINERS ---//

    nMuons = 0;                         nElectrons = 0;                     nLeptons = 0; 
    nLooseMuons = 0;                    nIsoMVAElectrons = 0;               nNoIsoMVAElectrons = 0;
    nHZZMuons = 0;                      nHZZElectrons = 0;                  nHZZLeptons = 0; 
    nGenMuons = 0;                      nGenElectrons = 0;                  nGenLeptons = 0; 

    muonsP4ptr.Delete();                muonsQ.clear();                     muonFiredLeg1.clear();              muonFiredLeg2.clear();
    muonIsGhost.clear();                muonIsLoose.clear();                muonIsHZZ.clear();
    muonEnergySF.clear();               muonHZZIDSF.clear();
    muonTrigEffLeg1Data.clear();        muonTrigEffLeg1MC.clear();          muonTrigEffLeg2Data.clear();        muonTrigEffLeg2MC.clear(); 
    muonCombIso.clear();                muonsTrkIso.clear();                muonD0.clear();                     muonDz.clear();                     
    muonSIP3d.clear();                  muonPtErr.clear();                  muonBestTrackType.clear();
    muonNMatchStn.clear();              muonNPixHits.clear();               muonNTkLayers.clear();
    muonIsPF.clear();                   muonIsGlobal.clear();               muonIsTracker.clear();
                                                                                                        
    electronsP4ptr.Delete();            electronsQ.clear();                 electronFiredLeg1.clear();          electronFiredLeg2.clear();
    electronIsGhost.clear();            electronPassIsoMVA.clear();         electronPassNoIsoMVA.clear();       electronIsHZZ.clear();
    electronEnergySF.clear();           electronHZZIDSF.clear();
    electronTrigEffLeg1Data.clear();    electronTrigEffLeg1MC.clear();      electronTrigEffLeg2Data.clear();    electronTrigEffLeg2MC.clear();
    electronCombIso.clear();            electronsTrkIso.clear();            electronD0.clear();                 electronDz.clear();
    electronSIP3d.clear();              electronScEta.clear();              electronIsoMVA.clear();             electronNoIsoMVA.clear();
    electronNMissHits.clear();          electronIsGap.clear();

    genMuonsP4ptr.Delete();             genMuonsQ.clear();                  genMuonStatus.clear();
    genElectronsP4ptr.Delete();         genElectronsQ.clear();              genElectronStatus.clear();




    /////////////////////
    //      START      //
    /////////////////////
    
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    if (sync_print)
    {
        if (fInfo->runNum == 275311 && fInfo->evtNum == 560373421)
        {
            cout << "START!" << endl;

            cout << "========================================" << endl;
            cout << "Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << std::endl;
            cout << "========================================\n" << endl;
        }
        else
            return kTRUE;
    }

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    weights->SetDataBit(isData);
   


    //--- EVENT WEIGHTS ---//

    genWeight   = 1.;
    PUWeight    = 1.;
    nPU         = 0;
    nPartons    = 0;
    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    nPV         = fPVArr->GetEntries();


    if (!isData)
    {
        // Save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (genWeight < 0)
            hTotalEvents->Fill(10);


        // Pileup reweighting
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(nPU);




        /////////////////////
        //  GEN PARTICLES  //
        /////////////////////


        hPhaseSpaceEvents->Fill(1, genWeight);

        vector<TLorentzVector> genLeps, genMuons, genElecs;     // Momenta
        vector<int> genLepsQ;                                   // Charges
        int elecTotalQ = 0, muonTotalQ = 0;



        //--- HARD PROCESS ---//

        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);


            // Parton counting for jet-binned Drell-Yan samples
            if (particle->status == 23 && particle->parent != -2 && (abs(particle->pdgId) < 6 || particle->pdgId == 21)) 
                nPartons++;


            // Find final state leptons
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13)
                && (particle->status == 23 || particle->status == 1 || particle->status == 2)
                && particle->parent != -2)
            {
                // Find if the lepton comes from a Z
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                int origin = abs(mother->pdgId);

                // (status 23 particle seems to be guaranteed as a hard scatter product
                //  and tends to be missing mother info)
                if (origin == 23 || particle->status == 23)
                {
                    TLorentzVector lep;
                    int q = copysign(1, particle->pdgId);

                    if (abs(particle->pdgId) == 13)
                    {   
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        new(genMuonsP4ptr[nGenMuons]) TLorentzVector(lep);

                        nGenMuons++;
                        genMuons.push_back(lep);
                        muonTotalQ += q;
                        genMuonsQ.push_back(q);
                        genMuonStatus.push_back(particle->status);
                    }
                    else if (abs(particle->pdgId) == 11)
                    {   
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        new(genElectronsP4ptr[nGenElectrons]) TLorentzVector(lep);

                        nGenElectrons++;
                        genElecs.push_back(lep);
                        elecTotalQ += q;
                        genElectronsQ.push_back(q);
                        genElectronStatus.push_back(particle->status);
                    }

                    nGenLeptons++;
                    genLeps.push_back(lep);
                    genLepsQ.push_back(q);
                }
            }
        }



        //--- PHASE SPACE ---//

        bool isPhaseSpace = kTRUE;
        TLorentzVector muonSum = GetP4Sum(genMuons), elecSum = GetP4Sum(genElecs), lepSum = GetP4Sum(genLeps);

        // Mass window
        if (lepSum.M() < M_MIN || lepSum.M() > M_MAX)
            isPhaseSpace = kFALSE;


        // Charge requirement
        if (elecTotalQ != 0 || muonTotalQ != 0)
            isPhaseSpace = kFALSE;


        // Sort events by decay channel
        unsigned idx = 0;

        if (nGenMuons == 2 && nGenElectrons == 0)           // mumu = 3
            idx = 3;

        else if (nGenMuons == 0 && nGenElectrons == 2)      // ee   = 4
            idx = 4;

        else if (nGenMuons == 4 && nGenElectrons == 0)      // 4m   = 6
            idx = 6;

        else if (nGenMuons == 2 && nGenElectrons == 2       // 2m2e = 7
                && muonSum.M() > elecSum.M())
            idx = 7;

        else if (nGenMuons == 2 && nGenElectrons == 2       // 2e2m = 8
                && muonSum.M() < elecSum.M())
            idx = 8;

        else if (nGenMuons == 0 && nGenElectrons == 4)      // 4e   = 9
            idx = 9;

        else
            isPhaseSpace = kFALSE;


        // Dilepton mass requirement
        if (nGenMuons == 2 && nGenElectrons == 2)       // Mixed flavor
        {
            if (elecSum.M() < MLL_MIN)
                isPhaseSpace = kFALSE;

            if (muonSum.M() < MLL_MIN)
                isPhaseSpace = kFALSE;
        }
        else if (nGenMuons == 4 || nGenElectrons == 4)  // Single flavor
        {
            for (unsigned j = 1; j < 4; j++)
            {
                for (unsigned i = 0; i < j; i++)
                {
                    if (genLepsQ[i] != genLepsQ[j])
                    {
                        TLorentzVector dilep = genLeps[i] + genLeps[j];
                        if (dilep.M() < MLL_MIN)
                            isPhaseSpace = kFALSE;
                    }
                }
            }
        }


        // Remaining events must be in phase space 
        if (isPhaseSpace)
        {
            hPhaseSpaceEvents->Fill(idx, genWeight);
            hFiducialEvents->Fill(1, genWeight);



            //--- FIDUCIAL REGION ---//

            bool isFiducial = kTRUE;
            sort(genLeps.begin(), genLeps.end(), P4SortCondition);


            // Eta
            for (unsigned i = 0; i < nGenLeptons; i++)
            {
                if (fabs(genLeps[i].Eta()) > ETA_MAX)
                    isFiducial = kFALSE;
            }


            // Lepton Pt
            if (genLeps[0].Pt() < PT1_MIN)
                isFiducial = kFALSE;

            if (genLeps[1].Pt() < PT2_MIN)
                isFiducial = kFALSE;

            if (nGenLeptons == 4)
            {
                if (genLeps[3].Pt() < PT_MIN || genLeps[2].Pt() < PT_MIN)
                    isFiducial = kFALSE;
            }


            // Remaining events must be in fiducial region
            if (isFiducial)
                hFiducialEvents->Fill(idx, genWeight);
        }
    }


    // Lumi mask
    if (isData)
    {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);



    //--- TRIGGER SELECTION ---//

    evtMuonTriggered = kFALSE;
    for (unsigned i = 0; i < muonTriggerNames.size(); i++)
    {
        if (trigger->pass(muonTriggerNames[i], fInfo->triggerBits))
            evtMuonTriggered = kTRUE;
    }

    evtElectronTriggered = kFALSE;
    for (unsigned i = 0; i < electronTriggerNames.size(); i++)
    {
        if (trigger->pass(electronTriggerNames[i], fInfo->triggerBits))
            evtElectronTriggered = kTRUE;
    }

    bool passTrigger = evtMuonTriggered || evtElectronTriggered;

    if (sync_print)
        cout << "trigger status: " << passTrigger << "\n" << endl;

    if (!passTrigger) // &isData)
        return kTRUE;

    if (passTrigger)    
        hTotalEvents->Fill(3);




    ////////////////////
    // SELECT OBJECTS //
    ////////////////////


    //--- VERTICES ---//

    if (fInfo->hasGoodPV)
    {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    }
    else
        return kTRUE;

    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);



    //--- MUONS ---//

    vector<TMuon*> muons;
    for (int i = 0; i < fMuonArr->GetEntries(); i++)
    {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);

        if (muon->pt > 4)
        {
            nMuons++;
            muons.push_back(muon);
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);


    // Second pass: store P4 & SF, ghost clean, identify loose & HZZ muons
    for (unsigned i = 0; i < nMuons; i++)
    {
        TMuon* muon = muons[i];
        TLorentzVector muonP4;
        copy_p4(muons[i], MUON_MASS, muonP4);

        if (sync_print)
            cout << "(" << muonP4.Pt() << ", " << muon->pt << ", " << muon->eta << ", " << muon->phi 
                 << ") , " << muon->q << ", " << muon->trkIso << endl;


        // Store P4 before corrections
        new(muonsP4ptr[i]) TLorentzVector(muonP4);


        // Apply Rochester correction
        muonEnergySF.push_back(GetRochesterCorrection(muon, muonCorr, rng, isData));
        muon->pt *= muonEnergySF.back();
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);


        // Ghost cleaning
        // Does not include track segment matching...
        float minDeltaR = 0.02;
        bool isGhost = kFALSE;
        for (unsigned j = 0; j < i; j++)
        {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR)
            {
                isGhost = kTRUE;
                break;
            }
        }
        muonIsGhost.push_back(isGhost);


        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i + 1; j < nMuons; j++)
        {
            TLorentzVector muon_j;
            copy_p4(muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3)
            {
                muon->trkIso = max(0., muon->trkIso - muon_j.Pt());
                muons[j]->trkIso = max(0., muons[j]->trkIso - muonP4.Pt());
            }
        }
        if (sync_print)
            cout << "(pt, pt_raw, eta, phi), q, trk_iso" << endl;


        // Check muon ID
        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
            nLooseMuons++;
        if (PassMuonHZZTightID(muon))
            nHZZMuons++;
    }
    if (sync_print) cout << endl;



    //--- ELECTRONS ---//

    vector<TElectron*> electrons;
    for (int i = 0; i < fElectronArr->GetEntries(); i++)
    {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (electron->pt > 4)
        {
            nElectrons++;
            electrons.push_back(electron);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);


    // Second pass: store P4 & SF; identify loose, HZZ, & tight electrons
    for (unsigned i = 0; i < nElectrons; i++)
    {
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);

        // Store P4 before corrections
        new(electronsP4ptr[i]) TLorentzVector(electronP4);

        // Apply electron energy correction
        electronEnergySF.push_back(GetElectronPtSF(electron));
        electronP4 *= electronEnergySF.back();


        // Cross cleaning
        float minDeltaR = 0.05;
        bool isGhost = kFALSE;
        for (unsigned j = 0; j < nMuons; j++)
        {
            TMuon* tmpMuon = muons[j];
            TLorentzVector tmpMuonP4;
            copy_p4(tmpMuon, MUON_MASS, tmpMuonP4);

            float dr = electronP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR && PassMuonHZZTightID(tmpMuon))
            {
                isGhost = kTRUE;
                break;
            }
        }
        electronIsGhost.push_back(isGhost);


        // Check electron ID
        if (particleSelector->PassElectronMVA(electron, cuts->hzzIsoV1))
            nIsoMVAElectrons++;
        if (particleSelector->PassElectronMVA(electron, cuts->hzzNoIsoV1))
            nNoIsoMVAElectrons++;
        if (PassElectronHZZTightID(electron, fInfo->rhoJet))
            nHZZElectrons++;
    }

    nLeptons        = nMuons + nElectrons;
    nHZZLeptons     = nHZZMuons + nHZZElectrons;



    //--- MET ---//

    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    if (!isData)
        met = 0.96*met; //met = MetKluge(met)*met;

    if (sync_print)
    {
        cout << "\npfmet, pfmet_type1" << endl;
        cout << fInfo->pfMET << ", " << fInfo->pfMETC << "\n" << endl;
    }




    ///////////////
    // SELECTION //
    ///////////////


    // Require at least two same-flavor leptons passing HZZ cuts
    if (nHZZMuons < 2 && nHZZElectrons < 2)
        return kTRUE;
    hTotalEvents->Fill(5);



    /////////////////////
    // FILL CONTAINERS //
    /////////////////////


    //--- MUONS ---//

    for (unsigned i = 0; i < nMuons; i++)
    {
        // Basic quantities
        TMuon* muon = muons[i];
        TLorentzVector muonP4;
        copy_p4(muons[i], MUON_MASS, muonP4);

        muonsQ.push_back(muon->q);


        // Boolean ID
        muonIsLoose.push_back(muon->pogIDBits & baconhep::kPOGLooseMuon);
        muonIsHZZ.push_back(PassMuonHZZTightID(muon));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZMuonIDEff(muon);
        muonHZZIDSF.push_back(effCont.GetSF());


        // Trigger bools and SFs
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE;
        for (unsigned i = 0; i < muonTriggerNames.size(); i++)
        {
            if (trigger->passObj(muonTriggerNames[i], 1, muon->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(muonTriggerNames[i], 2, muon->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        muonFiredLeg1.push_back(firedLeg1);
        muonFiredLeg2.push_back(firedLeg2);

        // Leptons will automatically fail leg 2 for single trigger selection
//      if (params->selection == "single")
//          effCont = weights->GetSingleMuonTriggerEff(muon);
//      else if (params->selection == "double")
            effCont = weights->GetDoubleMuonTriggerEff(muon, 1);
        trigEff = effCont.GetEff();

        if (!firedLeg1)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        muonTrigEffLeg1Data.push_back(trigEff.first);
        muonTrigEffLeg1MC.push_back(trigEff.second);

        effCont = weights->GetDoubleMuonTriggerEff(muon, 2);
        trigEff = effCont.GetEff();
        if (!firedLeg2)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        muonTrigEffLeg2Data.push_back(trigEff.first);
        muonTrigEffLeg2MC.push_back(trigEff.second);


        // ID criteria for downstream use
        muonCombIso.push_back(GetMuonIsolation(muon));
        muonsTrkIso.push_back(muon->trkIso);
        muonD0.push_back(muon->d0);
        muonDz.push_back(muon->dz);
        muonSIP3d.push_back(muon->sip3d);
        muonPtErr.push_back(muon->ptErr);
        muonNMatchStn.push_back(muon->nMatchStn);
        muonNPixHits.push_back(muon->nPixHits);
        muonNTkLayers.push_back(muon->nTkLayers);
        muonIsPF.push_back(muon->typeBits & baconhep::kPFMuon);
        muonIsGlobal.push_back(muon->typeBits & baconhep::kGlobal);
        muonIsTracker.push_back(muon->typeBits & baconhep::kTracker);
        muonBestTrackType.push_back(muon->btt);
    }



    //--- ELECTRONS ---//

    for (unsigned i = 0; i < nElectrons; i++)
    {
        // Basic quantities
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);

        electronsQ.push_back(electron->q);

        // Boolean ID
        electronPassIsoMVA.push_back(particleSelector->PassElectronMVA(electron, cuts->hzzIsoV1));
        electronPassNoIsoMVA.push_back(particleSelector->PassElectronMVA(electron, cuts->hzzNoIsoV1));
        electronIsHZZ.push_back(PassElectronHZZTightID(electron, fInfo->rhoJet));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZElectronIDRecoEff(electron);
        electronHZZIDSF.push_back(effCont.GetSF());


        // Trigger
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE;
        for (unsigned i = 0; i < electronTriggerNames.size(); i++)
        {
            if (trigger->passObj(electronTriggerNames[i], 1, electron->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(electronTriggerNames[i], 2, electron->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        electronFiredLeg1.push_back(firedLeg1);
        electronFiredLeg2.push_back(firedLeg2);

        effCont = weights->GetDoubleElectronTriggerEff(electron, 1);
        trigEff = effCont.GetEff();
        if (!firedLeg1)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        electronTrigEffLeg1Data.push_back(trigEff.first);
        electronTrigEffLeg1MC.push_back(trigEff.second);

        effCont = weights->GetDoubleElectronTriggerEff(electron, 2);
        trigEff = effCont.GetEff();
        if (!firedLeg2)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        electronTrigEffLeg2Data.push_back(trigEff.first);
        electronTrigEffLeg2MC.push_back(trigEff.second);


        // ID criteria
        electronIsoMVA.push_back(electron->mvaIso);
        electronNoIsoMVA.push_back(electron->mva);
        electronCombIso.push_back(GetElectronIsolation(electron, fInfo->rhoJet));
        electronsTrkIso.push_back(electron->trkIso);
        electronD0.push_back(electron->d0);
        electronDz.push_back(electron->dz);
        electronSIP3d.push_back(electron->sip3d);
        electronNMissHits.push_back(electron->nMissingHits);
        electronIsGap.push_back(electron->fiducialBits & kIsGap);
    }


    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void MultileptonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void MultileptonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void MultileptonAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("output") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<MultileptonAnalyzer> selector(new MultileptonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



float MultileptonAnalyzer::MetKluge(float met)
{
    if (met > 500) {
        return 1.;
    }

    float bins[]        = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 250., 500.};
    float corrections[] = {1., 0.99, 0.96, 0.95, 0.948, 0.947, 0.95, 0.957, 0.966, 0.97, 0.971, 0.98};

    int bin = 0;
    for (int i = 0; i < 12; ++i) {
        if (met > bins[i] && met <= bins[i+1]) {
            bin = i;
            break;
        }
    }
    
    return corrections[bin];
}



// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muons
float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}



// roccor.Run2.v1.tgz: https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
// From README:
// double mcSF = rc.kSpreadMC(Q, pt, eta, phi, genPt, s=0, m=0); //(recommended), MC scale and resolution correction when matched gen muon is available
// double mcSF = rc.kSmearMC(Q, pt, eta, phi, nl, u, s=0, m=0); //MC scale and extra smearing when matched gen muon is not available
float MultileptonAnalyzer::GetRochesterCorrection(const baconhep::TMuon* muon, RoccoR* muonCorr, TRandom3* rng, bool isData)
{   
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muon_scale_and_resolution_correc
    if (muon->pt < 200 && muon->btt == 1 && muon->nTkLayers > 5)
    {
        if (isData)
            return muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        else
            return muonCorr->kSmearMC(muon->q, muon->pt, muon->eta, muon->phi, muon->nTkLayers, rng->Rndm(), 0, 0);
    }
    else
        return 1;
}



// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muons
// No change since 2016
bool MultileptonAnalyzer::PassMuonHZZTightID(const baconhep::TMuon* muon)
{
    // AFTER Rochester correction!

    if (    muon->pt > 5
        &&  fabs(muon->eta) < 2.4
        &&  fabs(muon->d0) < 0.5        // These are, indeed, corrected: https://github.com/NWUHEP/BaconProd/blob/jbueghly_2017/Ntupler/src/FillerMuon.cc
        &&  fabs(muon->dz) < 1.0
        &&  ((muon->typeBits & baconhep::kGlobal) || ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) // Global muon or (arbitrated) tracker muon
        &&  muon->btt != 2
        &&  fabs(muon->sip3d) < 4.0                 
        &&  GetMuonIsolation(muon)/muon->pt < 0.35)
    {
        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
            return kTRUE;
        else if (   muon->pt > 200
                 && muon->nMatchStn > 1
                 && (muon->ptErr/muon->pt) < 0.3
                 && fabs(muon->d0) < 0.2
                 && fabs(muon->dz) < 0.5
                 && muon->nPixHits > 0
                 && muon->nTkLayers > 5)
            return kTRUE;
        else
            return kFALSE;
    }
    else
        return kFALSE;
}


// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Electrons
// seems to have typos?  So this is unchanged from 2016
//
// Effective area from
// https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
// which is NOT the page linked on the HZZ twiki!!!
float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 5.0};
    float effArea[7] = {0.1566, 0.1626, 0.1073, 0.0854, 0.1051, 0.1204, 0.1524};
    for (unsigned i = 0; i < 7; i++)
    {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1])
            iEta = i;
    }

    // These seem to be defined in CMSSW for 0.3 iso cone (see HZZ twiki) and they are pulled directly in the ntupler
    float combIso = el->chHadIso + std::max(0., (double)el->neuHadIso + el->gammaIso - rho*effArea[iEta]);
    return combIso;
}


// FIXME
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
float MultileptonAnalyzer::GetElectronPtSF(baconhep::TElectron* electron)
{
    TLorentzVector electronP4;
    electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);
    return electron->calibE / electronP4.E();
}


// FIXME
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Electrons
bool MultileptonAnalyzer::PassElectronHZZTightID(const baconhep::TElectron* electron, float rho)
{
    if (    electron->calibPt > 7
        &&  fabs(electron->scEta) < 2.5
//      &&  (particleSelector->PassElectronMVA(electron, cuts->hzzIsoV1) || particleSelector->PassElectronMVA(electron, cuts->hzzMVANoIsoV1))   // Needs to be checked separately!!
        &&  fabs(electron->d0) < 0.5
        &&  fabs(electron->dz) < 1
        &&  fabs(electron->sip3d) < 4.0
        &&  GetElectronIsolation(electron, rho)/electron->calibPt < 0.35)   // Not sure if this needs to be included for iso MVA?
        return kTRUE;
    else
        return kFALSE;
}
