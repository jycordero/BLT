#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_gen_pt (TGenParticle *i, TGenParticle *j) { return (i->pt > j->pt); }

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
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "mumu" || params->selection == "emu" || params->selection == "4mu") {
        //triggerNames.push_back("HLT_IsoMu22_v*");
        //triggerNames.push_back("HLT_IsoTkMu22_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + 
            "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
        //"/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");
    rng = new TRandom3();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);
    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);
    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeQ", &leptonThreeQ);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeTrigger", &leptonThreeTrigger);
    outTree->Branch("leptonFourP4", &leptonFourP4);
    outTree->Branch("leptonFourIso", &leptonFourIso);
    outTree->Branch("leptonFourQ", &leptonFourQ);
    outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
    outTree->Branch("leptonFourTrigger", &leptonFourTrigger);

    outTree->Branch("GENleptonOneP4", &muon1);
    outTree->Branch("GENleptonTwoP4", &muon2);
    outTree->Branch("GENleptonThreeP4", &muon3);
    outTree->Branch("GENleptonFourP4", &muon4);

    // jets
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetD0", &jetD0);
    outTree->Branch("jetTag", &jetTag);
    outTree->Branch("jetPUID", &jetPUID);
    outTree->Branch("jetFlavor", &jetFlavor);

    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("bjetTag", &bjetTag);
    outTree->Branch("bjetPUID", &bjetPUID);
    outTree->Branch("bjetFlavor", &bjetFlavor);

    outTree->Branch("genBJetP4", &genBJetP4);
    outTree->Branch("genBJetTag", &genBJetTag);
    outTree->Branch("genJetP4", &genJetP4);
    outTree->Branch("genJetTag", &genJetTag);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    // histograms
    hZEta = new TH1D("z_eta", "z_eta", 100, -10, 10);
    hZPt = new TH1D("z_pt", "z_pt", 100, 0, 100);
    hZMass = new TH1D("z_mass", "z_mass", 100, 80, 100);
    hMu1Eta = new TH1D("mu1_eta", "mu1_eta", 100, -10, 10);
    hMu1Pt = new TH1D("mu1_pt", "mu1_pt", 100, 0, 100);
    hMu2Eta = new TH1D("mu2_eta", "mu2_eta", 100, -10, 10);
    hMu2Pt = new TH1D("mu2_pt", "mu2_pt", 100, 0, 100);

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;

    if (sync_print) {
        if (
            (fInfo->runNum == 275311 && fInfo->evtNum == 560373421  )
           ) {
                    cout << "START!" << endl;

                    cout << "========================================" << endl;
                    cout << "Run: " << fInfo->runNum 
                        << " Lumi: " << fInfo->lumiSec 
                        << " Event: " << fInfo->evtNum 
                        << std::endl;
                    cout << "========================================\n" << endl;
                } else {
                    return kTRUE;
                }
    }

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    ////////////////////////
    // Gen-level analysis //
    ////////////////////////

    bool isProcess = false, isAccepted = false;
    if (!isData)
    {
        vector<TGenParticle*> leptons;
        vector<TGenParticle*> zbosons;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (/*particle->status == 2 &&*/ abs(particle->pdgId) == 23)
                zbosons.push_back(particle);
            else if (particle->status == 1 && abs(particle->pdgId) == 13)
                leptons.push_back(particle);
        }
        std::sort(leptons.begin(), leptons.end(), sort_gen_pt);
        if (zbosons.size() > 1)
        {
            float zmass_diff = 1000;
            for (unsigned int i = 0; i < zbosons.size(); i++)
            {
                if (abs(zbosons[i]->mass - 91.2) < zmass_diff)
                {
                    zmass_diff = abs(zbosons[i]->mass - 91.2);
//                    cout << zbosons[i]->status << endl;
                    the_z.SetPtEtaPhiM(zbosons[i]->pt, zbosons[i]->eta, zbosons[i]->phi, zbosons[i]->mass);
                }
            }
        }

        if (params->selection == "mumu")
        {
            if (leptons.size() >= 2)
            {
                muon1.SetPtEtaPhiM(leptons[0]->pt, leptons[0]->eta, leptons[0]->phi, leptons[0]->mass);
                for (unsigned i = 1; i < leptons.size(); ++i)
                {
                    if (leptons[0]->pdgId != leptons[i]->pdgId)
                    {
                        muon2.SetPtEtaPhiM(leptons[i]->pt, leptons[i]->eta, leptons[i]->phi, leptons[i]->mass);
                        break;
                    }
                }
                hZEta->Fill(the_z.Eta());
                hZPt->Fill(the_z.Pt());
                hZMass->Fill(the_z.M());
                hMu1Eta->Fill(muon1.Eta());
                hMu1Pt->Fill(muon1.Pt());
                hMu2Eta->Fill(muon2.Eta());
                hMu2Pt->Fill(muon2.Pt());
            }
        }
        else if (params->selection == "4mu")
        {
            if (params->datasetgroup == "DYJetsToLL" && leptons.size() > 1)
            {
                for (unsigned k = 0; k < zbosons.size(); k++)
                {
                    for (unsigned int j = 0; j < leptons.size(); j++)
                    {
                        for (unsigned int i = 0; i < j; i++)
                        {
                            float Mll = sqrt(2 * leptons[i]->pt * leptons[j]->pt * (cosh(leptons[i]->eta - leptons[j]->eta) - cos(leptons[i]->phi - leptons[j]->phi)));
                            if (abs(Mll - zbosons[k]->mass) <= 3 && leptons[i]->pdgId != leptons[j]->pdgId)
                            {
                                isProcess = true;
                                break;
                            }
                        }
                    }
                }
            }

            if (leptons.size() >= 4)
            {
                unsigned muon2Index = 1, muon3Index = 2, muon4Index = 3;
                muon1.SetPtEtaPhiM(leptons[0]->pt, leptons[0]->eta, leptons[0]->phi, leptons[0]->mass);
                if (leptons[0]->pdgId == leptons[1]->pdgId)
                {
                    muon2Index = 2;
                    muon3Index = 1;
                }
                muon2.SetPtEtaPhiM(leptons[muon2Index]->pt, leptons[muon2Index]->eta, leptons[muon2Index]->phi, leptons[muon2Index]->mass);
                muon3.SetPtEtaPhiM(leptons[muon3Index]->pt, leptons[muon3Index]->eta, leptons[muon3Index]->phi, leptons[muon3Index]->mass);

                for (unsigned i = muon4Index; i < leptons.size(); ++i)
                {
                    if (leptons[muon3Index]->pdgId != leptons[i]->pdgId)
                    {
                        muon4.SetPtEtaPhiM(leptons[i]->pt, leptons[i]->eta, leptons[i]->phi, leptons[i]->mass);
                        muon4Index = i;
                        if (params->datasetgroup == "Zto4mu")
                            isProcess = true;
                        break;
                    }
                }

            }

            if (isProcess)
            {
                this->genEvents++;
                float Mll = (muon1 + muon2).M(), M4l = (muon1 + muon2 + muon3 + muon4).M();
                if (muon1.Pt() > 25 && muon2.Pt() > 25 && muon3.Pt() > 10 && muon4.Pt() > 10
                        && abs(muon1.Eta()) < 2.1 && abs(muon2.Eta()) < 2.1 && abs(muon3.Eta()) < 2.4 && abs(muon4.Eta()) < 2.4
                        && Mll > 12 && M4l > 80 && M4l < 100)
                {
                    this->genAcceptedEvents++;
                    isAccepted = true;
                }

            }
        }


//      cout << fGenParticleArr->GetEntries() << ": ";
//      for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
//      {
//          TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
//          cout << particle->pdgId << " ";
//          int n = 0;

//          int j = particle->parent;
//          if (abs(particle->pdgId) == 13)
//          if (j >= 0)
//          {
//              TGenParticle* mom = (TGenParticle*) fGenParticleArr->At(j);
//              cout << mom->pdgId << " ";
//              if (mom->pdgId == 22)
//                  n++;
//          }
//          if (n > 0)
//              cout << n << endl;
//      }
//    cout << endl;
    }

    // Apply lumi mask
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
//        if(!lumiMask.HasRunLumi(rl))
//        { 
//            cout << "Lumi" << endl;
//            return kTRUE;
//        }
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }

    if (!passTrigger)// && isData)
        return kTRUE;

    if (sync_print) {
        cout << "trigger status: " << passTrigger << "\n" << endl;
    }

    hTotalEvents->Fill(3);

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
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << particle->status << ", "
            //     << particle->pdgId  << ", "
            //     << particle->parent
            //     << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;

            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
    }

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    /* Apply a preselection so we can make a collection of muons to clean against */
    vector<TMuon*> tmp_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 5 
                && fabs(muon->eta) < 2.4
                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
           ) {
            tmp_muons.push_back(muon);
        }
    }
    sort(tmp_muons.begin(), tmp_muons.end(), sort_by_higher_pt<TMuon>);

    // Second pass
    //int trigger_muon_index = -1;
    vector<TLorentzVector> muons;
    vector<TLorentzVector> veto_muons;
    vector<float> muons_iso;
    vector<float> muons_q;
    vector<bool> muons_trigger;

    if (sync_print) {
        cout << "(pt, pt_raw, eta, phi), q, trk_iso" << endl;
    }
    for (unsigned i = 0; i < tmp_muons.size(); i++) {
        TMuon* muon = tmp_muons[i];
        TLorentzVector muonP4;
        copy_p4(tmp_muons[i], MUON_MASS, muonP4);

        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i+1; j < tmp_muons.size(); j++) {
            TLorentzVector muon_j;
            copy_p4(tmp_muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3) {
                muon->trkIso = max(0., muon->trkIso - muon_j.Pt());
                tmp_muons[j]->trkIso = max(0., tmp_muons[j]->trkIso - muonP4.Pt());
            }
        }

        // Apply rochester muon momentum corrections
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                                                muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                                                0, 0);
        }
        muonP4.SetPtEtaPhiM(muonSF*muon->pt, muon->eta, muon->phi, MUON_MASS);

        if (sync_print) {
            cout << "(" << muonP4.Pt() << ", " << muon->pt 
                << ", " << muon->eta << ", " << muon->phi 
                << ") , " << muon->q 
                << ", " << muon->trkIso
                << endl;
        }

        // Fill containers
        if (muon->trkIso/muon->pt < 0.1) {

            if (muonP4.Pt() > 20)
                veto_muons.push_back(muonP4);

            if (muonP4.Pt() > 5) {
                muons.push_back(muonP4);
                muons_iso.push_back(muon->trkIso);
                muons_q.push_back(muon->q);

                // trigger matching
                bool triggered = false;
                for (unsigned i = 0; i < triggerNames.size(); ++i) {
                    triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
                }
                muons_trigger.push_back(triggered);
            }
        }
    }
    std::sort(muons.begin(), muons.end(), P4SortCondition);

    if (sync_print) cout << endl;

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<float> electrons_q;
    vector<bool> electrons_trigger;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 5 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
            electrons_iso.push_back(0.);
            electrons_q.push_back(electron->q);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);


    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    std::vector<TJet*> genbjets;
    std::vector<TJet*> genjets;
    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;

    if (sync_print) {
        cout << "(pt, pt_raw, eta, phi), id, overlap, csv, bmva, puid_bdt" << endl;
    }

    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
        }

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        //bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

        //cout << jet->genpt << ", " << jet->genm << ", " 
        //     << jet->geneta << ", " << jet->genphi << ", " 
        //     << jet->partonFlavor << ", " << jet->hadronFlavor << ", " 
        //     << endl;

        if (sync_print) {
            cout << "(" << jet->pt << ", " << jet->ptRaw << ", " 
                << jet->eta << ", " << jet->phi << "), " 
                << particleSelector->PassJetID(jet, cuts->looseJetID) << ", " 
                << muOverlap << ", " 
                << jet->csv << ", " << jet->bmva << ", "
                << jet->mva 
                << endl;
        }

        if (!isData) {
            if (abs(jet->hadronFlavor) == 5) {
                genbjets.push_back(jet);
            } else {
                genjets.push_back(jet);
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        && jet->mva > -0.89
                        && !muOverlap 
                        //&& !elOverlap
                   ) { 
                    if (isData) {
                        if (jet->bmva > 0.9432) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "MVAT")) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    }
                }
            } else {
                if ((fabs(jet->eta) < 2.5 && jet->mva > -0.89) || fabs(jet->eta) > 2.5) {
                    fwdjets.push_back(jet);
                    ++nFwdJets;
                }
            }
        }
    }

    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genjets.begin(), genjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genbjets.begin(), genbjets.end(), sort_by_higher_pt<TJet>);

    // Add additional b jets to the central jet collection
    if (bjets.size() > 1) {
        jets.insert(jets.end(), bjets.begin()+1, bjets.end());
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    if (!isData) {
        //met = MetKluge(met)*met;
        met = 0.96*met;
    }

    if (sync_print) {
        cout << "\npfmet, pfmet_type1" << endl;
        cout << fInfo->pfMET << ", "
            << fInfo->pfMETC
            << "\n" << endl;
    }

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    // Synchronization printout
    if (sync_print) {
        if (nBJets >= 1 && nJets >= 1 && met < 40 && muons.size() >= 2) {

            jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);

            TLorentzVector dijet = bjetP4 + jetP4;
            TLorentzVector dimuon = muons[0] + muons[1];
            cout << "phi_mumu, phi_jj, dphi_mumujj" << endl;
            cout << dimuon.Phi() << ", " << dijet.Phi() << ", " << fabs(dimuon.DeltaPhi(dijet)) << endl;
        }
        cout << "STOP!" << endl;
        return kTRUE;

    } else if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4;
        muonOneP4 = muons[0];
        unsigned muonTwoIndex = 1;
        for (unsigned i = 1; i < muons.size(); ++i) {
            if (muons_q[0] != muons_q[i]) {
                muonTwoP4 = muons[i];
                muonTwoIndex = i;
                break;
            }
        }
//        if (muons_q[0] * muons_q[muonTwoIndex] > 0)
//            return kTRUE;

        if (
                !(muonOneP4.Pt() > 25 && fabs(muonOneP4.Eta()) < 2.1) 
                || !(muonTwoP4.Pt() > 25 && fabs(muonTwoP4.Eta()) < 2.1)
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muonOneP4 + muonTwoP4;
        if (sync_print) {
            cout << dimuon.M() << endl;
        }

        if (dimuon.M() < 12.)
            return kTRUE;
        hTotalEvents->Fill(7);
/*
        if (nBJets == 0 || (nBJets + nJets < 2 && nFwdJets == 0)) 
            return kTRUE;
        hTotalEvents->Fill(8);
*/
        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 1;

        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso[muonTwoIndex];
        leptonTwoQ       = muons_q[muonTwoIndex];
        leptonTwoTrigger = muons_trigger[muonTwoIndex];
        leptonTwoFlavor  = 1;

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    } else if (params->selection == "4mu") {
        if (muons.size() < 4)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, muonThreeP4, muonFourP4;
        unsigned muonTwoIndex = 1, muonThreeIndex = 2, muonFourIndex = 3;
        muonOneP4 = muons[0];
        if (muons_q[0] == muons_q[1])
        {
            muonTwoIndex = 2;
            muonThreeIndex = 1;
        }
        muonTwoP4 = muons[muonTwoIndex];
        muonThreeP4 = muons[muonThreeIndex];
        for (unsigned i = muonFourIndex; i < muons.size(); ++i) {
            if (muons_q[muonThreeIndex] != muons_q[i]) {
                muonFourP4 = muons[i];
                muonFourIndex = i;
                break;
            }
        }
        if (muons_q[muonThreeIndex] == muons_q[muonFourIndex])
            return kTRUE;
/*        
        muonTwoIndex = 1; muonThreeIndex = 2; muonFourIndex = 3;
        if (muons_q[0] * muons_q[1] > 0 || muons_q[2] * muons_q[3] > 0)
            return kTRUE;
        muonOneP4 = muons[0];
        muonTwoP4 = muons[1];
        muonThreeP4 = muons[2];
        muonFourP4 = muons[3];
*/
        if (
                !(muonOneP4.Pt() > 25 && fabs(muonOneP4.Eta()) < 2.1) 
                || !(muonTwoP4.Pt() > 25 && fabs(muonTwoP4.Eta()) < 2.1)
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muonOneP4 + muonTwoP4;
        if (sync_print) {
            cout << dimuon.M() << endl;
        }

        if (dimuon.M() < 12)
            return kTRUE;
        hTotalEvents->Fill(7);
/*
        TLorentzVector fourmuon;
        fourmuon = muonOneP4 + muonTwoP4 + muonThreeP4 + muonFourP4;
        if (
                fourmuon.M() < 80 || fourmuon.M() > 100
                || muonThreeP4.Pt() < 10 || muonFourP4.Pt() < 10
           )
            return kTRUE;
*/
/*
        if (nBJets == 0 || (nBJets + nJets < 2 && nFwdJets == 0)) 
            return kTRUE;
        hTotalEvents->Fill(8);
*/
        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 1;

        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso[1];
        leptonTwoQ       = muons_q[1];
        leptonTwoTrigger = muons_trigger[1];
        leptonTwoFlavor  = 1;

        leptonThreeP4      = muonThreeP4;
        leptonThreeIso     = muons_iso[2];
        leptonThreeQ       = muons_q[2];
        leptonThreeTrigger = muons_trigger[2];
        leptonThreeFlavor  = 1;

        leptonFourP4      = muonFourP4;
        leptonFourIso     = muons_iso[3];
        leptonFourQ       = muons_q[3];
        leptonFourTrigger = muons_trigger[3];
        leptonFourFlavor  = 1;

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // not adjusted for four leptons
            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    } else if (params->selection == "ee") {

        if (electrons.size() != 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector dielectron;
        dielectron = electrons[0] + electrons[1];
        if (dielectron.M() < 12. || dielectron.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electrons[0];
        leptonOneIso     = electrons_iso[0];
        leptonOneQ       = electrons_q[0];
        leptonOneTrigger = electrons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[1];
        leptonTwoIso     = electrons_iso[1];
        leptonTwoQ       = electrons_q[1];
        leptonTwoTrigger = electrons_trigger[1];
        leptonTwoFlavor  = 0;
    } else if (params->selection == "emu") {

        if (muons.size() != 1 || electrons.size() != 1)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0].Pt() < 25 || fabs(muons[0].Eta()) > 2.1)
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dilepton;
        dilepton = muons[0] + electrons[0];
        if (dilepton.M() < 12 || dilepton.M() > 70)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4      = muons[0];
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[0];
        leptonTwoIso     = electrons_iso[0];
        leptonTwoQ       = electrons_q[0];
        leptonTwoTrigger = electrons_trigger[0];
        leptonTwoFlavor  = 1;

        if (!isData && false) {
            eventWeight *= weights->GetMuonIDEff(muons[0]);

            // trigger efficiency
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muons[0]);
            eventWeight *= trigEff.first;
        }
    } 
        


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0     = bjets[0]->d0;
        bjetTag    = bjets[0]->csv;
        bjetPUID   = bjets[0]->mva;
        bjetFlavor = bjets[0]->hadronFlavor;
    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0     = 0.;
        bjetTag    = 0.;
        bjetPUID   = 0.;
        bjetFlavor = 0.;
    }

    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0     = fwdjets[0]->d0;
        jetTag    = 0.;
        jetPUID   = fwdjets[0]->mva;
        jetFlavor = fwdjets[0]->hadronFlavor;
    } else if (jets.size() > 0) {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0     = jets[0]->d0;
        jetTag    = jets[0]->csv;
        jetPUID   = jets[0]->mva;
        jetFlavor = jets[0]->hadronFlavor;
    } else {
        jetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetD0     = 0.;
        jetTag    = 0.;
        jetPUID   = 0.;
        jetFlavor = 0.;
    } 

    if (genbjets.size() > 0 && !isData) {
        genBJetP4.SetPtEtaPhiM(genbjets[0]->genpt, genbjets[0]->geneta, genbjets[0]->genphi, genbjets[0]->genm);
        genBJetTag = genbjets[0]->csv;
    } else {
        genBJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genBJetTag = 0;
    }

    if (genjets.size() > 0 && !isData) {
        genJetP4.SetPtEtaPhiM(genjets[0]->genpt, genjets[0]->geneta, genjets[0]->genphi, genjets[0]->genm);
        genJetTag = genjets[0]->csv;
    } else {
        genJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetTag = 0;
    }

    outTree->Fill();
    this->passedEvents++;

    if (!isData & isAccepted)
        this->passIdEvents++;

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
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;

    if (fInfo->runNum == 1)
    {
        std::cout << "           : Accepted " << this->genAcceptedEvents << " / " << this->genEvents << " gen-level events." << std::endl;
        std::cout << "           : Found " << this->passIdEvents << " / " << this->genAcceptedEvents << " events passing ID requirements." << std::endl << std::endl;
        std::cout << "  Acceptance: " << (float)this->genAcceptedEvents / (float)this->genEvents << std::endl;
        std::cout << "  Efficiency: " << (float)this->passIdEvents / (float)this->genAcceptedEvents << std::endl;
    }

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
