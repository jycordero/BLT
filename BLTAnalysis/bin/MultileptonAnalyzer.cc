#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
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


    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));


    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2";
    trigger.reset(new baconhep::TTrigger(trigfilename));


    // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2012
    // http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2013/108
    if (params->selection == "double")
    {
        muonTriggerNames.push_back("HLT_Mu17_Mu8_v*");
        muonTriggerNames.push_back("HLT_Mu17_TkMu8_v*");
        electronTriggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }
    else if (params->selection == "single")
    {
        muonTriggerNames.push_back("HLT_IsoMu24_eta2p1_v*");
        muonTriggerNames.push_back("HLT_IsoMu24_v*");
        electronTriggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
//      electronTriggerNames.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*");    // Not so sure about this one...
    }


    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false));


    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);


    // muon momentum corrections
    muonCorr = new rochcor2012();


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

    outTree->Branch("nMVAElectrons", &nMVAElectrons);
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
    outTree->Branch("electronPassMVA", &electronPassMVA);
    outTree->Branch("electronIsHZZ", &electronIsHZZ);

    outTree->Branch("electronEnergySF", &electronEnergySF);
    outTree->Branch("electronHZZIDSF", &electronHZZIDSF);
    outTree->Branch("electronTrigEffLeg1Data", &electronTrigEffLeg1Data);
    outTree->Branch("electronTrigEffLeg1MC", &electronTrigEffLeg1MC);
    outTree->Branch("electronTrigEffLeg2Data", &electronTrigEffLeg2Data);
    outTree->Branch("electronTrigEffLeg2MC", &electronTrigEffLeg2MC);

    outTree->Branch("electronMVA", &electronMVA);
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
    outHistName = params->get_output_treename("AcceptedEvents");
    hAcceptedEvents = new TH1D(outHistName.c_str(), "AcceptedEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{


    //--- CLEAR CONTAINERS ---//

    nMuons = 0;                         nElectrons = 0;                     nLeptons = 0; 
    nLooseMuons = 0;                    nMVAElectrons = 0;
    nHZZMuons = 0;                      nHZZElectrons = 0;                  nHZZLeptons = 0; 
    nGenMuons = 0;                      nGenElectrons = 0;                  nGenLeptons = 0; 

    muonsP4ptr.Delete();                muonsQ.clear();                     muonFiredLeg1.clear();              muonFiredLeg2.clear();
    muonIsGhost.clear();                muonIsLoose.clear();                muonIsHZZ.clear();
    muonEnergySF.clear();               muonHZZIDSF.clear();
    muonTrigEffLeg1Data.clear();        muonTrigEffLeg1MC.clear();          muonTrigEffLeg2Data.clear();        muonTrigEffLeg2MC.clear(); 
    muonCombIso.clear();                muonsTrkIso.clear();                muonD0.clear();                     muonDz.clear();                     
    muonSIP3d.clear();                  muonPtErr.clear();
    muonNMatchStn.clear();              muonNPixHits.clear();               muonNTkLayers.clear();              muonBestTrackType.clear();
    muonIsPF.clear();                   muonIsGlobal.clear();               muonIsTracker.clear();

    electronsP4ptr.Delete();            electronsQ.clear();                 electronFiredLeg1.clear();          electronFiredLeg2.clear();
    electronIsGhost.clear();            electronPassMVA.clear();            electronIsHZZ.clear();
    electronEnergySF.clear();           electronHZZIDSF.clear();
    electronTrigEffLeg1Data.clear();    electronTrigEffLeg1MC.clear();      electronTrigEffLeg2Data.clear();    electronTrigEffLeg2MC.clear();
    electronCombIso.clear();            electronsTrkIso.clear();            electronD0.clear();                 electronDz.clear();
    electronSIP3d.clear();              electronMVA.clear();
    electronNMissHits.clear();          electronIsGap.clear();

    genMuonsP4ptr.Delete();             genMuonsQ.clear();                  genMuonStatus.clear();
    genElectronsP4ptr.Delete();         genElectronsQ.clear();              genElectronStatus.clear();




    /////////////////////
    //      START      //
    /////////////////////


    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0) { 
        std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;
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
        unsigned muCount = 0, eCount = 0;
        int idSum = 0;
        TLorentzVector lepSum, muSum, eSum;
        hAcceptedEvents->Fill(1, genWeight);

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);


            // Parton counting for jet-binned Drell-Yan samples
            if (particle->status == 3 && particle->parent != -2 && (abs(particle->pdgId) < 6 || particle->pdgId == 21))
                nPartons++;


            // Find final state leptons
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13)
                    && particle->status == 3
                    && particle->parent != -2)
            {
                // Find if the lepton comes from a Z
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                int origin = abs(mother->pdgId);

                if (origin == 23)
                {
                    idSum += particle->pdgId;

                    TLorentzVector lepVec;
                    if (abs(particle->pdgId) == 11)
                    {
                        eCount++;
                        lepVec.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        eSum += lepVec;
                    }
                    else if (abs(particle->pdgId) == 13)
                    {
                        muCount++;
                        lepVec.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        muSum += lepVec;
                    }

                    lepSum += lepVec;
                }
            }
        }

        if (idSum == 0 && lepSum.M() < 100 && lepSum.M() > 80)
        {
            if (muCount == 2 && eCount == 0)            // mumu = 3
                hAcceptedEvents->Fill(3, genWeight);
            else if (muCount == 0 && eCount == 2)       // ee   = 4
                hAcceptedEvents->Fill(4, genWeight);
            else if (muCount == 4 && eCount == 0)       // 4m   = 6
                hAcceptedEvents->Fill(6, genWeight);
            else if (muCount == 2 && eCount == 2        // 2m2e = 7
                    && muSum.M() > eSum.M())
                hAcceptedEvents->Fill(7, genWeight);
            else if (muCount == 2 && eCount == 2        // 2e2m = 8
                    && muSum.M() < eSum.M())
                hAcceptedEvents->Fill(8, genWeight);
            else if (muCount == 0 && eCount == 4)       // 4e   = 9
                hAcceptedEvents->Fill(9, genWeight);
        }

        // Sanity check
        if (kFALSE && muCount + eCount != 2)
        {
            cout << "Index \tStatus\tID    \tParent" << endl;
            for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
            {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
                cout << i << "\t" << particle->status << "\t" << particle->pdgId << "\t" << particle->parent << endl;
            }           
            cout << endl << endl;
        }


        // Pileup reweighting
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(nPU);
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


        // Store P4 before corrections
        new(muonsP4ptr[i]) TLorentzVector(muonP4);


        // Apply Rochester correction
        muonEnergySF.push_back(GetRochesterCorrection(muon, muonCorr, isData));
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
                muon->trkIso03 = max(0., muon->trkIso03 - muon_j.Pt());
                muons[j]->trkIso03 = max(0., muons[j]->trkIso03 - muonP4.Pt());
            }
        }


        // Check muon ID
        if (IsPOGLooseMuon(muon))
            nLooseMuons++;
        if (PassMuonHZZTightID(muon, muonP4))
            nHZZMuons++;
    }



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


    // Second pass: store P4 & SF, cross clean, identify MVA & HZZ electrons
    for (unsigned i = 0; i < nElectrons; i++)
    {
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);


        // Store P4 before corrections
        new(electronsP4ptr[i]) TLorentzVector(electronP4);


        // Apply electron energy correction
        electronEnergySF.push_back(GetElectronPtSF(electron));
        electronP4.SetPtEtaPhiM(electron->ptHZZ4l, electron->eta, electron->phi, ELE_MASS);


        // Cross cleaning
        float minDeltaR = 0.05;
        bool isGhost = kFALSE;
        for (unsigned j = 0; j < nMuons; j++)
        {
            TMuon* tmpMuon = muons[j];
            TLorentzVector tmpMuonP4;
            copy_p4(tmpMuon, MUON_MASS, tmpMuonP4);

            float dr = electronP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR && PassMuonHZZTightID(tmpMuon, tmpMuonP4))
            {
                isGhost = kTRUE;
                break;
            }
        }
        electronIsGhost.push_back(isGhost);


        // Check electron ID
        if (particleSelector->PassElectronMVA(electron, cuts->hzzMVAID))
            nMVAElectrons++;
        if (PassElectronHZZTightID(electron, fInfo->rhoJet))
            nHZZElectrons++;
    }

    nLeptons        = nMuons + nElectrons;
    nHZZLeptons     = nHZZMuons + nHZZElectrons;



    //--- MET ---//

    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;




    ///////////////
    // SELECTION //
    ///////////////


    // Require at least two HZZ leptons of any flavor
    if (nHZZLeptons < 2)
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


        // ID bools and SF
        muonIsLoose.push_back(IsPOGLooseMuon(muon));
        muonIsHZZ.push_back(PassMuonHZZTightID(muon, muonP4));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZMuonIDEff(muonP4);
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
        if (params->selection == "single")
            effCont = weights->GetSingleMuonTriggerEff(muonP4);
        else if (params->selection == "double")
            effCont = weights->GetDoubleMuonTriggerEff(muonP4, 1);
        trigEff = effCont.GetEff();

        if (!firedLeg1)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        muonTrigEffLeg1Data.push_back(trigEff.first);
        muonTrigEffLeg1MC.push_back(trigEff.second);

        effCont = weights->GetDoubleMuonTriggerEff(muonP4, 2);
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
        muonsTrkIso.push_back(muon->trkIso03);
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

        int btt = ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker)) ? 0 : 2;
        muonBestTrackType.push_back(btt);
    }



    //--- ELECTRONS ---//

    for (unsigned i = 0; i < nElectrons; i++)
    {
        // Basic quantities
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);

        electronsQ.push_back(electron->q);


        // ID bools and SF
        electronPassMVA.push_back(particleSelector->PassElectronMVA(electron, cuts->hzzMVAID));
        electronIsHZZ.push_back(PassElectronHZZTightID(electron, fInfo->rhoJet));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZElectronIDRecoEff(electron);
        electronHZZIDSF.push_back(effCont.GetSF());


        // Trigger bools and SFs
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

        // Leptons will automatically fail leg 2 for single trigger selection
        if (params->selection == "single")
            effCont = weights->GetSingleElectronTriggerEff(electronP4);
        else if (params->selection == "double")
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


        // ID criteria for downstream use
        electronMVA.push_back(electron->mva);
        electronCombIso.push_back(GetElectronIsolation(electron, fInfo->rhoJet));
        electronsTrkIso.push_back(electron->trkIso03);
        electronD0.push_back(electron->d0);
        electronDz.push_back(electron->dz);
        electronSIP3d.push_back(electron->sip3d);
        electronScEta.push_back(electron->scEta);
        electronNMissHits.push_back(electron->nMissingHits);
        electronIsGap.push_back(electron->fiducialBits & kIsGap);
    }



    //--- GEN PARTICLES ---//

    if (!isData)
    {
        unsigned muonIdx = 0, electronIdx = 0;

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            // Sanity check
            if (kFALSE)
            {
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                std::cout << i << ", " << particle->status << ", " << particle->pdgId  << ", "
                    << particle->parent << ", " << mom.Pt() << endl;
            }

            // Returns all leptons regardless of status and parent
            if (abs(particle->pdgId) == 13)
            {
                nGenMuons++;
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                new(genMuonsP4ptr[muonIdx]) TLorentzVector(mom);
                genMuonsQ.push_back(copysign(1, particle->pdgId));
                genMuonStatus.push_back(particle->status);
                muonIdx++;
            }
            else if (abs(particle->pdgId) == 11)
            {
                nGenElectrons++;
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                new(genElectronsP4ptr[electronIdx]) TLorentzVector(mom);
                genElectronsQ.push_back(copysign(1, particle->pdgId));
                genElectronStatus.push_back(particle->status);
                electronIdx++;
            }

            nGenLeptons = nGenMuons + nGenElectrons;

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


// DIFFERENT than
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Muon_isolation
// which uses 0.4 isolation cone
float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}


// DIFFERENT than
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Momentum_scale_and_resolution
// which usese MuScle
//
// The Rochester correction function for 2012 returns void...
// but modifies the passed TLorentzVector by scaling each Px, Py, Pz, E by "pt_tune/pt"
float MultileptonAnalyzer::GetRochesterCorrection(const baconhep::TMuon* muon, rochcor2012* muonCorr, const bool isData)
{
    TLorentzVector muonP4;
    muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

    float qter = 1.;    // This is essentially a dummy variable
    if (isData)
        muonCorr->momcor_data(muonP4, muon->q, 0, qter);
    else
        muonCorr->momcor_mc(muonP4, muon->q, 0, qter);

    return muonP4.Pt() / muon->pt;
}


// DIFFERENT than
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2017#Muons
//
// DIFFERENT than 13TeV analysis
// Does not include "bestTrackType != 2", but this is apparently covered by "kGlobal || kTRacker":
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
bool MultileptonAnalyzer::PassMuonHZZTightID(const baconhep::TMuon* muon, const TLorentzVector &muonP4)
{
    // muonP4 should include Rochester correction! 

    if (    muonP4.Pt() > 5
            &&  fabs(muonP4.Eta()) < 2.4 
            &&  ((muon->typeBits & baconhep::kGlobal) || ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) // Global muon or (arbitrated) tracker mun
//          &&  muon->btt != 2                                          // NOT required for legacy
            &&  fabs(muon->d0) < 0.5                                    // NOT corrected for legacy; I am pretty sure they are corrected in the Bacon however
            &&  fabs(muon->dz) < 1.0
            &&  fabs(muon->sip3d) < 4.0                 
            &&  GetMuonIsolation(muon)/muonP4.Pt()  < 0.35              // 0.4 for legacy; isolation is also defined differently
       )
    {   
        if (muon->typeBits & baconhep::kPFMuon)
            return kTRUE;

        else if (   muonP4.Pt() > 200                                   // Pretty sure this part isn't in the legacy
                && (muon->ptErr/muon->pt) < 0.3
                && muon->nMatchStn > 1
                && muon->nPixHits > 0
                && fabs(muon->d0) < 0.2
                && fabs(muon->dz) < 0.5
                && muon->nTkLayers > 5)
            return kTRUE;
        else
            return kFALSE;
    }
    else
        return kFALSE;
}


bool MultileptonAnalyzer::IsPOGLooseMuon(const baconhep::TMuon* muon)
{
    if  (   (muon->typeBits & baconhep::kPFMuon)
            && ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker))
        )
        return kTRUE;
    else
        return kFALSE;
}


// DIFFERENT than
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Electron_isolation
// which uses 0.4 isolation cone
//
// Effective area from
// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Isolation_cone_R_0_3
float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
{
    unsigned iEta = 0;
    float etaBins[8] = {0,      1,      1.479,  2,      2.2,    2.3,    2.4,    2.5};
    float effArea[7] = {0.13,   0.14,   0.07,   0.09,   0.11,   0.11,   0.14};

    for (unsigned i = 0; i < 7; i++)
    {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1])
            iEta = i;
    }

    float combIso = el->chHadIso03 + std::max(0., (double)el->neuHadIso03 + el->gammaIso03 - rho*effArea[iEta]);
    return combIso;
}


// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#E_scale_corrections_Smearing
// Looks like this (regression and correction) is embedded in the Bacon:
// https://github.com/NWUHEP/BaconProd/blob/04/Ntupler/src/FillerElectron.cc
float MultileptonAnalyzer::GetElectronPtSF(baconhep::TElectron* electron)
{
    return electron->ptHZZ4l / electron->pt;
}


// DIFFERENT than
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4lPaper2013#Electrons 
bool MultileptonAnalyzer::PassElectronHZZTightID(const baconhep::TElectron* electron, const float rho)
{
    if (    electron->ptHZZ4l > 7
            &&  fabs(electron->scEta) < 2.5
    //      &&  particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)     // Check separately downstream!
            &&  fabs(electron->d0) < 0.5
            &&  fabs(electron->dz) < 1
            &&  fabs(electron->sip3d) < 4.0
            &&  electron->nMissingHits < 2                                      // Needed for 2012 (later included in MVA)
            &&  GetElectronIsolation(electron, rho)/electron->ptHZZ4l < 0.35    // 0.4 used in legacy
       )
        return kTRUE;
    else
        return kFALSE;
}
