#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

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
/*
    if (params->selection == "mumu" || params->selection == "emu")
    {
        //triggerNames.push_back("HLT_IsoMu22_v*");
        //triggerNames.push_back("HLT_IsoTkMu22_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    }
    else if (params->selection == "ee")
    {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }
*/

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
                                    // "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // electron scale corrections
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");


    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // Branches
    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("nPV", &nPV);

    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("PUWeight", &PUWeight);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    outTree->Branch("evtMuonTriggered", &evtMuonTriggered);
    outTree->Branch("evtElectronTriggered", &evtElectronTriggered);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nLeptons", &nLeptons);

    outTree->Branch("nVetoElectrons", &nVetoElectrons);

    outTree->Branch("nLooseMuons", &nLooseMuons);
    outTree->Branch("nLooseElectrons", &nLooseElectrons);
    outTree->Branch("nLooseLeptons", &nLooseLeptons);

    outTree->Branch("nMediumElectrons", &nMediumElectrons);

    outTree->Branch("nTightMuons", &nTightMuons);
    outTree->Branch("nTightElectrons", &nTightElectrons);
    outTree->Branch("nTightLeptons", &nTightLeptons);

    outTree->Branch("nHZZMuons", &nHZZMuons);
    outTree->Branch("nHZZElectrons", &nHZZElectrons);
    outTree->Branch("nHZZLeptons", &nHZZLeptons);

    outTree->Branch("nGenMuons", &nGenMuons);
    outTree->Branch("nGenElectrons", &nGenElectrons);
    outTree->Branch("nGenLeptons", &nGenLeptons);

    // muons
    outTree->Branch("muonP4", &muonsP4, 32000, 1);
    outTree->Branch("muonQ", &muonsQ);
    outTree->Branch("muonCombIso", &muonCombIso);
    outTree->Branch("muonTrkIso", &muonsTrkIso);

    outTree->Branch("muonIsPF", &muonIsPF);
    outTree->Branch("muonIsGlobal", &muonIsGlobal);
    outTree->Branch("muonIsTracker", &muonIsTracker);
    outTree->Branch("muonIsLoose", &muonIsLoose);
    outTree->Branch("muonIsTight", &muonIsTight);
    outTree->Branch("muonIsHZZ", &muonIsHZZ);
    outTree->Branch("muonTriggered", &muonTriggered);

    outTree->Branch("muonSF", &muonSF);
    outTree->Branch("muonLooseIDEff", &muonLooseIDEff);
    outTree->Branch("muonTightIDEff", &muonTightIDEff);
    outTree->Branch("muonLooseIsoEff", &muonLooseIsoEff);
    outTree->Branch("muonTightIsoEff", &muonTightIsoEff);
    outTree->Branch("muonHZZIDEff", &muonHZZIDEff);
    outTree->Branch("muonTriggerEffData", &muonTriggerEffData);
    outTree->Branch("muonTriggerEffMC", &muonTriggerEffMC);

    outTree->Branch("muonD0", &muonD0);
    outTree->Branch("muonDz", &muonDz);
    outTree->Branch("muonSIP3d", &muonSIP3d);
    outTree->Branch("muonMuNChi2", &muonMuNChi2);
    outTree->Branch("muonPtErr", &muonPtErr);
    outTree->Branch("muonNMatchStn", &muonNMatchStn);
    outTree->Branch("muonNPixHits", &muonNPixHits);
    outTree->Branch("muonNTkLayers", &muonNTkLayers);
    outTree->Branch("muonBestTrackType", &muonBestTrackType);
    outTree->Branch("muonNValidHits", &muonNValidHits);

    // electrons
    outTree->Branch("electronP4", &electronsP4, 32000, 1);
    outTree->Branch("electronQ", &electronsQ);
    outTree->Branch("electronCombIso", &electronCombIso);
    outTree->Branch("electronTrkIso", &electronsTrkIso);

    outTree->Branch("electronPassVetoIso", &electronPassVetoIso);
    outTree->Branch("electronPassLooseIso", &electronPassLooseIso);
    outTree->Branch("electronPassMediumIso", &electronPassMediumIso);
    outTree->Branch("electronPassTightIso", &electronPassTightIso);
    outTree->Branch("electronIsVeto", &electronIsVeto);
    outTree->Branch("electronIsLoose", &electronIsLoose);
    outTree->Branch("electronIsMedium", &electronIsMedium);
    outTree->Branch("electronIsTight", &electronIsTight);
    outTree->Branch("electronIsHZZ", &electronIsHZZ);
    outTree->Branch("electronTriggered", &electronTriggered);

    outTree->Branch("electronSF", &electronSF);
    outTree->Branch("electronHZZRecoEff", &electronHZZRecoEff);
    outTree->Branch("electronRecoEff", &electronRecoEff);
    outTree->Branch("electronTriggerEffData", &electronTriggerEffData);
    outTree->Branch("electronTriggerEffMC", &electronTriggerEffMC);

    outTree->Branch("electronD0", &electronD0);
    outTree->Branch("electronDz", &electronDz);
    outTree->Branch("electronSIP3d", &electronSIP3d);
    outTree->Branch("electronScEta", &electronScEta);
    outTree->Branch("electronSieie", &electronSieie);
    outTree->Branch("electronEnergyInv", &electronEnergyInv);
    outTree->Branch("electronHOverE", &electronHOverE);
    outTree->Branch("electronDEtaIn", &electronDEtaIn);
    outTree->Branch("electronDPhiIn", &electronDPhiIn);
    outTree->Branch("electronNMissHits", &electronNMissHits);
    outTree->Branch("electronIsConv", &electronIsConv);

    // gen-level particles
    outTree->Branch("genMuonP4", &genMuonsP4, 32000, 1);
    outTree->Branch("genMuonQ", &genMuonsQ);
    outTree->Branch("genMuonStatus", &genMuonStatus);

    outTree->Branch("genElectronP4", &genElectronsP4, 32000, 1);
    outTree->Branch("genElectronQ", &genElectronsQ);
    outTree->Branch("genElectronStatus", &genElectronStatus);

//  outTree->Branch("genIntermID", &genIntermID);
//  outTree->Branch("genIntermMass", &genIntermMass);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    string outHistName2 = params->get_output_treename("AcceptedEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);
    hAcceptedEvents = new TH1D(outHistName2.c_str(), "AcceptedEvents", 12, -1.5, 10.5);


    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{
    /* Clear counters & containers */
    nMuons = 0;                         nElectrons = 0;                     nLeptons = 0; 
    nLooseMuons = 0;                    nLooseElectrons = 0;                nLooseLeptons = 0; 
    nVetoElectrons = 0;                 nMediumElectrons = 0;
    nHZZMuons = 0;                      nHZZElectrons = 0;                  nHZZLeptons = 0; 
    nTightMuons = 0;                    nTightElectrons = 0;                nTightLeptons = 0; 
    nGenMuons = 0;                      nGenElectrons = 0;                  nGenLeptons = 0; 

    muonsP4ptr.Delete();                muonsQ.clear();                     muonCombIso.clear();                muonsTrkIso.clear();
    muonIsPF.clear();                   muonIsGlobal.clear();               muonIsTracker.clear();
    muonIsLoose.clear();                muonIsTight.clear();
    muonIsHZZ.clear();                  muonTriggered.clear(); 
    muonSF.clear();                     muonLooseIDEff.clear();             muonTightIDEff.clear();             muonHZZIDEff.clear();
    muonLooseIsoEff.clear();            muonTightIsoEff.clear();            muonTriggerEffData.clear();         muonTriggerEffMC.clear();
    muonD0.clear();                     muonDz.clear();                     muonSIP3d.clear();
    muonMuNChi2.clear();                muonPtErr.clear();
    muonNMatchStn.clear();              muonNPixHits.clear();               muonNTkLayers.clear();
    muonBestTrackType.clear();          muonNValidHits.clear();
                                                                                                        
    electronsP4ptr.Delete();            electronsQ.clear();                 electronCombIso.clear();            electronsTrkIso.clear();
    electronPassVetoIso.clear();        electronPassLooseIso.clear();       electronPassMediumIso.clear();      electronPassTightIso.clear();
    electronIsVeto.clear();             electronIsLoose.clear();            electronIsMedium.clear();           electronIsTight.clear();
    electronIsHZZ.clear();              electronTriggered.clear(); 
    electronSF.clear();                 electronRecoEff.clear();            electronHZZRecoEff.clear();     
    electronTriggerEffData.clear();     electronTriggerEffMC.clear();
    electronD0.clear();                 electronDz.clear();                 electronSIP3d.clear();
    electronScEta.clear();              electronSieie.clear();
    electronEnergyInv.clear();          electronHOverE.clear();             electronDEtaIn.clear();             electronDPhiIn.clear();
    electronNMissHits.clear();          electronIsConv.clear();

    genMuonsP4ptr.Delete();             genMuonsQ.clear();                  genMuonStatus.clear();
    genElectronsP4ptr.Delete();         genElectronsQ.clear();              genElectronStatus.clear();
//  genIntermID.clear();                genIntermMass.clear();




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
   


    /* Event weights */
    eventWeight   = 1.;
    PUWeight      = 1.;
    nPU           = 0;
    nPartons      = 0;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();

    if (!isData)
    {
        // Set data period for 2016 MC scale factors
        if (rng->Rndm() < 0.468)
            weights->SetDataPeriod("2016BtoF");    
        else
            weights->SetDataPeriod("2016GH");


        // Save gen weight for amc@nlo Drell-Yan sample
        eventWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (eventWeight < 0)
            hTotalEvents->Fill(10);


        unsigned lepCount = 0, idSum = 0;
        TLorentzVector lepSum;
        hAcceptedEvents->Fill(-1, eventWeight);

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
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

                if (origin == 23 || particle->status == 23)
                {
                    lepCount++;
                    idSum += particle->pdgId;

                    double pMass = 0;
                    if (abs(particle->pdgId) == 11)
                        pMass = ELE_MASS;
                    else if (abs(particle->pdgId) == 13)
                        pMass = MUON_MASS;

                    TLorentzVector lepVec;
                    lepVec.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, pMass);
                    lepSum += lepVec;
                }
            }
        }
        if (idSum == 0 && lepSum.M() < 100 && lepSum.M() > 80)
            hAcceptedEvents->Fill(lepCount, eventWeight);
/* 
        if (lepCount != 2)
        {
            cout << "Index \tStatus\tID    \tParent" << endl;
            for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
            {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
                cout << i << "\t" << particle->status << "\t" << particle->pdgId << "\t" << particle->parent << endl;
            }           
            cout << endl << endl;
        }
*/


        // Pileup reweighting
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(fInfo->nPUmean);
    }



    /* Lumi mask */
    if (isData)
    {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);



    /* Trigger selection */
    evtMuonTriggered     = trigger->pass("HLT_IsoMu24_v*", fInfo->triggerBits) || trigger->pass("HLT_IsoTkMu24_v*", fInfo->triggerBits);
    evtElectronTriggered = trigger->pass("HLT_Ele27_WPTight_Gsf_v*", fInfo->triggerBits);

    bool passTrigger = evtMuonTriggered || evtElectronTriggered;

    if (sync_print)
        cout << "trigger status: " << passTrigger << "\n" << endl;

    if (!passTrigger) // &isData)
        return kTRUE;

    if (passTrigger)    
        hTotalEvents->Fill(3);




    ///////////////////
    // Select objects//
    ///////////////////

    /* VERTICES */
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



    /* MUONS */
    vector<TMuon*> muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++)
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

    // Second pass: store P4 & SF; identify loose, HZZ, & tight muons
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
        muonSF.push_back(GetRochesterCorrection(muon, muonCorr, rng, isData));
        muon->pt *= muonSF.back();
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

        // Remove muons with very small deltaR
        float minDeltaR = 1e6;
        for (unsigned j = 0; j < nMuons; j++)
        {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1051);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR)
                minDeltaR = dr;
        }

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
        if (PassMuonTightID(muon))
            nTightMuons++;
        if (PassMuonHZZTightID(muon))
            nHZZMuons++;
    }
    if (sync_print) cout << endl;



    /* ELECTRONS */
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
        electronSF.push_back(GetElectronPtSF(electron, electronScaler, rng, runNumber));
        electron->pt *= electronSF.back();
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);

        // Check electron ID
        if (particleSelector->PassElectronID(electron, cuts->vetoElID))
            nVetoElectrons++;
        if (particleSelector->PassElectronID(electron, cuts->looseElID))
            nLooseElectrons++;
        if (particleSelector->PassElectronID(electron, cuts->mediumElID))
            nMediumElectrons++;
        if (particleSelector->PassElectronID(electron, cuts->tightElID))
            nTightElectrons++;
        if (PassElectronHZZTightID(electron, particleSelector, cuts, fInfo->rhoJet))
            nHZZElectrons++;
    }

    nLeptons        = nMuons + nElectrons;
    nLooseLeptons   = nLooseMuons + nLooseElectrons;
    nTightLeptons   = nTightMuons + nTightElectrons;
    nHZZLeptons     = nHZZMuons + nHZZElectrons;



    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    if (!isData)
        met = 0.96*met; //met = MetKluge(met)*met;

    if (sync_print)
    {
        cout << "\npfmet, pfmet_type1" << endl;
        cout << fInfo->pfMET << ", " << fInfo->pfMETC << "\n" << endl;
    }




    ///////////////////////////////
    // Apply analysis selections //
    ///////////////////////////////

    if (params->selection == "emu")
    {
        // Require at least two HZZ leptons of any flavor
        if (nHZZLeptons < 2)
            return kTRUE;
        hTotalEvents->Fill(5);



        /* MUONS */
        for (unsigned i = 0; i < nMuons; i++)
        {
            TMuon* muon = muons[i];
            TLorentzVector muonP4;
            copy_p4(muons[i], MUON_MASS, muonP4);

            muonsQ.push_back(muon->q);
            muonCombIso.push_back(GetMuonIsolation(muon));
            muonsTrkIso.push_back(muon->trkIso);

            muonIsPF.push_back(muon->typeBits & baconhep::kPFMuon);
            muonIsGlobal.push_back(muon->typeBits & baconhep::kGlobal);
            muonIsTracker.push_back(muon->typeBits & baconhep::kTracker);

            muonIsLoose.push_back(muon->pogIDBits & baconhep::kPOGLooseMuon);
            muonIsTight.push_back(PassMuonTightID(muon));
            muonIsHZZ.push_back(PassMuonHZZTightID(muon));
            muonTriggered.push_back(trigger->passObj("HLT_IsoMu24_v*", 1, muon->hltMatchBits) || trigger->passObj("HLT_IsoTkMu24_v*", 1, muon->hltMatchBits));

            muonLooseIDEff.push_back(weights->GetLooseMuonIDEff(muonP4));
            muonTightIDEff.push_back(weights->GetMuonIDEff(muonP4));
            muonLooseIsoEff.push_back(weights->GetMuonLooseISOEff(muonP4));
            muonTightIsoEff.push_back(weights->GetMuonTightISOEff(muonP4));
            muonHZZIDEff.push_back(weights->GetHZZMuonIDEff(muonP4));

            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
            muonTriggerEffData.push_back(trigEff.first);
            muonTriggerEffMC.push_back(trigEff.second);

            muonD0.push_back(muon->d0);
            muonDz.push_back(muon->dz);
            muonSIP3d.push_back(muon->sip3d);
            muonMuNChi2.push_back(muon->muNchi2);
            muonPtErr.push_back(muon->ptErr);
            muonNMatchStn.push_back(muon->nMatchStn);
            muonNPixHits.push_back(muon->nPixHits);
            muonNTkLayers.push_back(muon->nTkLayers);
            muonBestTrackType.push_back(muon->btt);
            muonNValidHits.push_back(muon->nValidHits);
        }



        /* ELECTRONS */
        for (unsigned i = 0; i < nElectrons; i++)
        {
            TElectron* electron = electrons[i];
            TLorentzVector electronP4;
            copy_p4(electrons[i], ELE_MASS, electronP4);

            electronsQ.push_back(electron->q);
            electronCombIso.push_back(GetElectronIsolation(electron, fInfo->rhoJet));
            electronsTrkIso.push_back(electron->trkIso);

            electronPassVetoIso.push_back(particleSelector->PassElectronIso(electron, cuts->vetoElIso, cuts->EAEl));
            electronPassLooseIso.push_back(particleSelector->PassElectronIso(electron, cuts->looseElIso, cuts->EAEl));
            electronPassMediumIso.push_back(particleSelector->PassElectronIso(electron, cuts->mediumElIso, cuts->EAEl));
            electronPassTightIso.push_back(particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl));

            electronIsVeto.push_back(particleSelector->PassElectronID(electron, cuts->vetoElID));
            electronIsLoose.push_back(particleSelector->PassElectronID(electron, cuts->looseElID));
            electronIsMedium.push_back(particleSelector->PassElectronID(electron, cuts->mediumElID));
            electronIsTight.push_back(particleSelector->PassElectronID(electron, cuts->tightElID));
            electronIsHZZ.push_back(PassElectronHZZTightID(electron, particleSelector, cuts, fInfo->rhoJet));
            electronTriggered.push_back(trigger->passObj("HLT_Ele27_WPTight_Gsf_v*", 1, electron->hltMatchBits));

            electronRecoEff.push_back(weights->GetElectronRecoEff(electronP4));
            electronHZZRecoEff.push_back(weights->GetHZZElectronRecoEff(*electron));

            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
            electronTriggerEffData.push_back(trigEff.first);
            electronTriggerEffMC.push_back(trigEff.second);

            electronD0.push_back(electron->d0);
            electronDz.push_back(electron->dz);
            electronSIP3d.push_back(electron->sip3d);
            electronScEta.push_back(electron->scEta);
            electronSieie.push_back(electron->sieie);
            electronEnergyInv.push_back(fabs(1. - electron->eoverp) / electron->ecalEnergy);
            electronHOverE.push_back(electron->hovere);
            electronDEtaIn.push_back(electron->dEtaIn);
            electronDPhiIn.push_back(electron->dPhiIn);
            electronNMissHits.push_back(electron->nMissingHits);
            electronIsConv.push_back(electron->isConv);
        }



        /* GEN PARTICLES */
        if (!isData)
        {
            unsigned muonIdx = 0, electronIdx = 0;

            for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
            {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

                if (kFALSE)
                {
                    TLorentzVector mom;
                    mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                    std::cout << i << ", " << particle->status << ", " << particle->pdgId  << ", "
                        << particle->parent << ", " << mom.Pt() << endl;
                }

                if (abs(particle->pdgId) == 13) // && particle->status == 23
                {
                    nGenMuons++;
                    TLorentzVector mom;
                    mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                    new(genMuonsP4ptr[muonIdx]) TLorentzVector(mom);
                    genMuonsQ.push_back(copysign(1, particle->pdgId));
                    genMuonStatus.push_back(particle->status);
                    muonIdx++;
                }
                else if (abs(particle->pdgId) == 11) // && particle->status == 23
                {
                    nGenElectrons++;
                    TLorentzVector mom;
                    mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                    new(genElectronsP4ptr[electronIdx]) TLorentzVector(mom);
                    genElectronsQ.push_back(copysign(1, particle->pdgId));
                    genElectronStatus.push_back(particle->status);
                    electronIdx++;
                }
//              else if (particle->status == 22)
//              {
//                  genIntermID.push_back(particle->pdgId);
//                  genIntermMass.push_back(particle->mass);
//              }

                nGenLeptons = nGenMuons + nGenElectrons;

            }
        }
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



float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}



float MultileptonAnalyzer::GetRochesterCorrection(const baconhep::TMuon* muon, RoccoR* muonCorr, TRandom3* rng, bool isData)
{
    if (isData)
        return muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
    else
        return muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                  muon->nTkLayers, rng->Rndm(), rng->Rndm(), 0, 0);
}



bool MultileptonAnalyzer::PassMuonTightID(const baconhep::TMuon* muon)
{
    TLorentzVector muonP4;
    muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

    if (    muonP4.Pt()     > 10.   &&  fabs(muonP4.Eta()) < 2.4
        &&  (muon->typeBits & baconhep::kPFMuon) && (muon->typeBits & baconhep::kGlobal)
        &&  muon->muNchi2    < 10.  &&  muon->nMatchStn  > 1    &&  muon->nPixHits   > 0
        &&  fabs(muon->d0)   < 0.2  &&  fabs(muon->dz)   < 0.5
        &&  muon->nTkLayers  > 5    &&  muon->nValidHits > 0
        &&  GetMuonIsolation(muon)/muonP4.Pt()  < 0.15)
        return kTRUE;
    else
        return kFALSE;
}



bool MultileptonAnalyzer::PassMuonHZZTightID(const baconhep::TMuon* muon)
{
    TLorentzVector muonP4;
    muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

    if (    muonP4.Pt()     > 5.    &&  fabs(muonP4.Eta()) < 2.4
        &&  (((muon->typeBits & baconhep::kGlobal) || 
             ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) &&
               (muon->btt != 2)) // Global muon or (arbitrated) tracker muon
        &&  fabs(muon->d0)  < 0.5   &&  fabs(muon->dz) < 1.0
        &&  fabs(muon->sip3d) < 4.0                 
        &&  GetMuonIsolation(muon)/muonP4.Pt()  < 0.35)
    {
        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
            return kTRUE;
        else if (   muonP4.Pt() > 200.0 && (muon->ptErr/muon->pt) < 0.3
                 && muon->nMatchStn > 1 && muon->nPixHits > 0
                 && fabs(muon->d0)  < 0.2   &&  fabs(muon->dz)  < 0.5
                 && muon->nTkLayers > 5)
            return kTRUE;
        else
            return kFALSE;
    }
    else
        return kFALSE;
}



float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, float rho)
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



float MultileptonAnalyzer::GetElectronPtSF(baconhep::TElectron* electron, EnergyScaleCorrection* electronScaler, TRandom3* rng, int runNumber)
{
    bool isData = (runNumber != 1);

    if (isData)
    {
        scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
        return sdata.scale;
    }
    else
    {
        float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
        return rng->Gaus(1, sFactor);
    }
}



bool MultileptonAnalyzer::PassElectronGoodID(const baconhep::TElectron* electron, std::unique_ptr<ParticleSelector>& particleSelector, std::unique_ptr<Cuts>& cuts)
{
    if (    electron->pt > 10.  &&  fabs(electron->scEta) < 2.5
        &&  particleSelector->PassElectronID(electron, cuts->tightElID)
        &&  particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl))
        return kTRUE;
    else
        return kFALSE;
}



bool MultileptonAnalyzer::PassElectronHZZTightID(const baconhep::TElectron* electron, std::unique_ptr<ParticleSelector>& particleSelector, std::unique_ptr<Cuts>& cuts, float rho)
{
    TLorentzVector electronP4;
    electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);

    if (    electron->pt > 5.   &&  fabs(electron->scEta) < 2.5
//      &&  particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)
        &&  fabs(electron->d0)  < 0.5   &&  fabs(electron->dz)  < 1.0
        &&  fabs(electron->sip3d) < 4.0 
        &&  GetElectronIsolation(electron, rho)/electronP4.Pt() < 0.35)
        return kTRUE;
    else
        return kFALSE;

}
