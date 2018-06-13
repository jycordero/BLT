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

//  if (params->selection == "mumu" || params->selection == "emu") {
//      //triggerNames.push_back("HLT_IsoMu22_v*");
//      //triggerNames.push_back("HLT_IsoTkMu22_v*");
//      triggerNames.push_back("HLT_IsoMu24_v*");
//      triggerNames.push_back("HLT_IsoTkMu24_v*");
//      //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
//      //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
//      //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
//      //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

//  } else if (params->selection == "ee") {
//      triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
//  }

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

    // electron scale corrections
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");

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
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("PUWeight", &PUWeight);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("passTrigger", &passTrigger);
    outTree->Branch("passHLT_IsoMu24", &passHLT_IsoMu24);
    outTree->Branch("passHLT_IsoTkMu24", &passHLT_IsoTkMu24);
    outTree->Branch("passHLT_Ele27_WPTight_Gsf", &passHLT_Ele27_WPTight_Gsf);

    // muons
    outTree->Branch("muonP4", &muonsP4, 32000, 1);
    outTree->Branch("muonQ", &muonsQ);
    outTree->Branch("muonIDEff", &muonIDEff);
    outTree->Branch("muonTightIsoEff", &muonTightIsoEff);
    outTree->Branch("muonLooseIsoEff", &muonLooseIsoEff);
    outTree->Branch("muonTriggerEffData", &muonTriggerEffData);
    outTree->Branch("muonTriggerEffMC", &muonTriggerEffMC);
    outTree->Branch("muonSF", &muonSF);
    outTree->Branch("muonTrkIso", &muonsTrkIso);
    outTree->Branch("muonIsPFMuon", &muonIsPFMuon);
    outTree->Branch("muonIsGLB", &muonIsGLB);
    outTree->Branch("muonMuNChi2", &muonMuNChi2);
    outTree->Branch("muonNMatchStn", &muonNMatchStn);
    outTree->Branch("muonNPixHits", &muonNPixHits);
    outTree->Branch("muonD0", &muonD0);
    outTree->Branch("muonDz", &muonDz);
    outTree->Branch("muonNTkLayers", &muonNTkLayers);
    outTree->Branch("muonNValidHits", &muonNValidHits);
    outTree->Branch("muonPassStdCuts", &muonPassStdCuts);
    outTree->Branch("muonPassTrigger", &muonPassTrigger);

    // electrons
    outTree->Branch("electronP4", &electronsP4, 32000, 1);
    outTree->Branch("electronQ", &electronsQ);
    outTree->Branch("electronRecoEff", &electronRecoEff);
    outTree->Branch("electronTriggerEffData", &electronTriggerEffData);
    outTree->Branch("electronTriggerEffMC", &electronTriggerEffMC);
    outTree->Branch("electronSF", &electronSF);
    outTree->Branch("electronTrkIso", &electronsTrkIso);
    outTree->Branch("electronCombIso", &electronCombIso);
    outTree->Branch("electronEnergyInv", &electronEnergyInv);
    outTree->Branch("electronScEta", &electronScEta);
    outTree->Branch("electronD0", &electronD0);
    outTree->Branch("electronDz", &electronDz);
    outTree->Branch("electronSieie", &electronSieie);
    outTree->Branch("electronHOverE", &electronHOverE);
    outTree->Branch("electronDEtaIn", &electronDEtaIn);
    outTree->Branch("electronDPhiIn", &electronDPhiIn);
    outTree->Branch("electronNMissHits", &electronNMissHits);
    outTree->Branch("electronIsConv", &electronIsConv);
    outTree->Branch("electronPassID", &electronPassID);
    outTree->Branch("electronPassIso", &electronPassIso);
    outTree->Branch("electronPassStdCuts", &electronPassStdCuts);
    outTree->Branch("electronPassTrigger", &electronPassTrigger);

    // gen-level particles
    outTree->Branch("genMuonP4", &genMuonsP4, 32000, 1);
    outTree->Branch("genMuonQ", &genMuonsQ);
    outTree->Branch("genElectronP4", &genElectronsP4, 32000, 1);
    outTree->Branch("genElectronQ", &genElectronsQ);
    outTree->Branch("genIntermID", &genIntermID);
    outTree->Branch("genIntermMass", &genIntermMass);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nLeptons", &nLeptons);
    outTree->Branch("nStdMuons", &nStdMuons);
    outTree->Branch("nStdElectrons", &nStdElectrons);
    outTree->Branch("nStdLeptons", &nStdLeptons);
    outTree->Branch("nGenMuons", &nGenMuons);
    outTree->Branch("nGenElectrons", &nGenElectrons);
    outTree->Branch("nGenLeptons", &nGenLeptons);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{
    // Clear vectors and TClonesArrays
    muonsP4ptr.Clear(); electronsP4ptr.Clear();
    muonsTrkIso.clear(); electronsTrkIso.clear();
    muonsQ.clear(); electronsQ.clear();
    muonIsPFMuon.clear(); muonIsGLB.clear(); muonPassStdCuts.clear(); muonPassTrigger.clear();
    muonIDEff.clear(); muonTightIsoEff.clear(); muonLooseIsoEff.clear();
    muonTriggerEffData.clear(); muonTriggerEffMC.clear();
    muonSF.clear(); muonMuNChi2.clear(); muonD0.clear(); muonDz.clear();
    muonNMatchStn.clear(); muonNPixHits.clear(); muonNTkLayers.clear(); muonNValidHits.clear();
    electronIsConv.clear(); electronPassID.clear(); electronPassIso.clear();
    electronPassStdCuts.clear(); electronPassTrigger.clear();
    electronRecoEff.clear(); electronTriggerEffData.clear(); electronTriggerEffMC.clear();
    electronSF.clear(); electronCombIso.clear(); electronEnergyInv.clear();
    electronScEta.clear(); electronD0.clear(); electronDz.clear(); electronSieie.clear();
    electronHOverE.clear(); electronDEtaIn.clear(); electronDPhiIn.clear();
    electronNMissHits.clear();
    genMuonsP4ptr.Clear(); genElectronsP4ptr.Clear();
    genMuonsQ.clear(); genElectronsQ.clear();
    genIntermID.clear(); genIntermMass.clear();

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

    // Apply lumi mask
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    passHLT_IsoMu24 = trigger->pass("HLT_IsoMu24_v*", fInfo->triggerBits);
    passHLT_IsoTkMu24 = trigger->pass("HLT_IsoTkMu24_v*", fInfo->triggerBits);
    passHLT_Ele27_WPTight_Gsf = trigger->pass("HLT_Ele27_WPTight_Gsf_v*", fInfo->triggerBits);
    passTrigger = (passHLT_IsoMu24 || passHLT_IsoTkMu24) || passHLT_Ele27_WPTight_Gsf;

    if (sync_print) {
        cout << "trigger status: " << passTrigger << "\n" << endl;
    }

//  if (!passTrigger && isData)
//      return kTRUE;

    if (passTrigger)
        hTotalEvents->Fill(3);


    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    PUWeight      = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        // Pileup reweighting
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(fInfo->nPUmean);

        // Save gen weight for amc@nlo Drell-Yan sample
        eventWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (eventWeight < 0)
            hTotalEvents->Fill(10);

        // Set data period for 2016 MC scale factors
        if (rng->Rndm() < 0.468)
            weights->SetDataPeriod("2016BtoF");    
        else
            weights->SetDataPeriod("2016GH");
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    if (!isData) {
        unsigned count = 0,  muonIdx = 0, electronIdx = 0;

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

//          TLorentzVector mom;
//          mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
//          cout << i << ", "
//               << particle->status << ", "
//               << particle->pdgId  << ", "
//               << particle->parent << ", "
//               << mom.Pt()
//               << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;

            } else if (
                    particle->status == 23
                    && abs(particle->pdgId) == 13
               ) {
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                new(genMuonsP4ptr[muonIdx]) TLorentzVector(mom);
        		genMuonsQ.push_back(copysign(1, particle->pdgId));
                muonIdx++;

            } else if (
                    particle->status == 23
                    && abs(particle->pdgId) == 11
               ) {
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                new(genElectronsP4ptr[electronIdx]) TLorentzVector(mom);
	        	genElectronsQ.push_back(copysign(1, particle->pdgId));
                electronIdx++;

            } else if (particle->status == 22)
            {
                genIntermID.push_back(particle->pdgId);
                genIntermMass.push_back(particle->mass);
            }

            nGenMuons     = genMuonsQ.size();
            nGenElectrons = genElectronsQ.size();
            nGenLeptons   = nGenMuons + nGenElectrons;

        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons      = 0;
        nGenMuons     = 0;
        nGenElectrons = 0;
        nGenLeptons   = 0;
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
    vector<TMuon*> all_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (muon->pt > 4) {
            all_muons.push_back(muon);
            muonPassStdCuts.push_back(kFALSE);
        }
    }
    sort(all_muons.begin(), all_muons.end(), sort_by_higher_pt<TMuon>);

    /* Apply a preselection so we can make a collection of muons to clean against */
    vector<TMuon*> tmp_muons;
    vector<unsigned> tmp_muons_idx;
    for (unsigned i = 0; i < all_muons.size(); i++) {
        TMuon* muon = all_muons[i];
        assert(muon);

        if (
                muon->pt > 10 
                && fabs(muon->eta) < 2.4
                // tight muon ID
                && (muon->typeBits & baconhep::kPFMuon) 
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
            tmp_muons_idx.push_back(i);
        }
    }

    // Second pass
    //int trigger_muon_index = -1;
    nStdMuons = 0;
    vector<unsigned> std_muons_idx;

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
        double muScale = 1.;
        if (isData) {
            muScale = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muScale = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                                                muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                                                0, 0);
        }
        muonP4.SetPtEtaPhiM(muScale*muon->pt, muon->eta, muon->phi, MUON_MASS);

        if (sync_print) {
            cout << "(" << muonP4.Pt() << ", " << muon->pt 
                << ", " << muon->eta << ", " << muon->phi 
                << ") , " << muon->q 
                << ", " << muon->trkIso
                << endl;
        }

        // Perform selection
        if (muon->trkIso/muonP4.Pt() < 0.1) {
//          if (muonP4.Pt() > 25 && fabs(muonP4.Eta()) < 2.1) {
            if (muonP4.Pt() > 20) {
                nStdMuons++;
                std_muons_idx.push_back(tmp_muons_idx[i]);
            }
        }
    }
    for (unsigned i = 0; i < std_muons_idx.size(); i++)
        muonPassStdCuts[std_muons_idx[i]] = kTRUE;

    if (sync_print) cout << endl;


    /* ELECTRONS */
    vector<TElectron*> all_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (electron->pt > 4)
            all_electrons.push_back(electron);
    }
    sort(all_electrons.begin(), all_electrons.end(), sort_by_higher_pt<TElectron>);

    nStdElectrons = 0;
    for (unsigned i = 0; i < all_electrons.size(); i++) {
        TElectron* electron = all_electrons[i];
        assert(electron);

        // Apply electron energy correction
        float eScale = 1.;
        if (isData) {
            scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
            eScale = sdata.scale;
        } else {
            float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
            eScale = rng->Gaus(1, sFactor);
        }
        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(eScale*electron->pt, electron->eta, electron->phi, ELE_MASS);

        if (
                electronP4.Pt() > 20 
                && fabs(electronP4.Eta()) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            nStdElectrons++;
            electronPassStdCuts.push_back(kTRUE);
        }
        else
            electronPassStdCuts.push_back(kFALSE);
    }

    nStdLeptons = nStdMuons + nStdElectrons;

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

    nMuons     = all_muons.size();
    nElectrons = all_electrons.size();
    nLeptons   = nMuons + nElectrons;

    // Synchronization printout
    if (params->selection == "emu") {
        if (nStdLeptons < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Fill containers
        for (unsigned i = 0; i < all_muons.size(); ++i) {
            TMuon* muon = all_muons[i];
            TLorentzVector muonP4;
            copy_p4(all_muons[i], MUON_MASS, muonP4);

            new(muonsP4ptr[i]) TLorentzVector(muonP4);
            muonsQ.push_back(muon->q);
            muonsTrkIso.push_back(muon->trkIso);
            muonIsPFMuon.push_back(muon->typeBits & baconhep::kPFMuon);
            muonIsGLB.push_back(muon->typeBits & baconhep::kGlobal);
            muonMuNChi2.push_back(muon->muNchi2);
            muonNMatchStn.push_back(muon->nMatchStn);
            muonNPixHits.push_back(muon->nPixHits);
            muonD0.push_back(muon->d0);
            muonDz.push_back(muon->dz);
            muonNTkLayers.push_back(muon->nTkLayers);
            muonNValidHits.push_back(muon->nValidHits);
            muonPassTrigger.push_back(trigger->passObj("HLT_IsoMu24_v*", 1, muon->hltMatchBits)
                    || trigger->passObj("HLT_IsoTkMu24_v*", 1, muon->hltMatchBits));

            if (isData)
            {
                muonSF.push_back(muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0));
                muonIDEff.push_back(1.);
                muonTightIsoEff.push_back(1.);      muonLooseIsoEff.push_back(1.);
                muonTriggerEffData.push_back(1.);   muonTriggerEffMC.push_back(1.);
            }
            else
            {
                muonSF.push_back(muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi, muon->nTkLayers, rng->Rndm(), rng->Rndm(), 0, 0));
                muonIDEff.push_back(weights->GetMuonIDEff(muonP4));
                muonTightIsoEff.push_back(weights->GetMuonTightISOEff(muonP4));
                muonLooseIsoEff.push_back(weights->GetMuonLooseISOEff(muonP4));
                pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                muonTriggerEffData.push_back(trigEff.first);    muonTriggerEffMC.push_back(trigEff.second);
            }
        }

        for (unsigned i = 0; i < all_electrons.size(); ++i) {
            TElectron* electron = all_electrons[i];
            TLorentzVector electronP4;
            copy_p4(all_electrons[i], ELE_MASS, electronP4);

            new(electronsP4ptr[i]) TLorentzVector(electronP4);
            electronsQ.push_back(electron->q);
            electronsTrkIso.push_back(electron->trkIso);
            electronEnergyInv.push_back(fabs(1. - electron->eoverp)/electron->ecalEnergy);
            electronScEta.push_back(electron->scEta);
            electronD0.push_back(electron->d0);
            electronDz.push_back(electron->dz);
            electronSieie.push_back(electron->sieie);
            electronHOverE.push_back(electron->hovere);
            electronDEtaIn.push_back(electron->dEtaIn);
            electronDPhiIn.push_back(electron->dPhiIn);
            electronNMissHits.push_back(electron->nMissingHits);
            electronIsConv.push_back(electron->isConv);
            electronPassID.push_back(particleSelector->PassElectronID(electron, cuts->tightElID));
            electronPassIso.push_back(particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl));
            electronPassTrigger.push_back(trigger->passObj("HLT_Ele27_WPTight_Gsf_v*", 1, electron->hltMatchBits));

            int iEta = 0;
            float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
            for (unsigned i = 0; i < 8; ++i) {
                if (fabs(electron->scEta) > etaBins[i] && fabs(electron->scEta) < etaBins[i+1]) {
                    iEta = i;
                    break;
                }
            }
            electronCombIso.push_back(electron->chHadIso + std::max(0.,(double)electron->neuHadIso + electron->gammaIso - fInfo->rhoJet * cuts->EAEl[iEta]));

            float eScale;
            if (isData)
            {
                scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
                eScale = sdata.scale;
                electronRecoEff.push_back(1.);
                electronTriggerEffData.push_back(1.);   electronTriggerEffMC.push_back(1.);
            }
            else
            {
                float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
                eScale = rng->Gaus(1, sFactor);
                electronRecoEff.push_back(weights->GetElectronRecoEff(electronP4));
                pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                electronTriggerEffData.push_back(trigEff.first);    electronTriggerEffMC.push_back(trigEff.second);
            }
            electronSF.push_back(eScale);
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
