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

//  // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
//  std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
//  trigger.reset(new baconhep::TTrigger(trigfilename));

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
//  weights.reset(new WeightUtils(params->period, "mumu", false)); // Lumi mask
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
    string outFileName = params->get_output_filename("skimmed");
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
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // muons
    outTree->Branch("muonP4", &muonsP4, 256000, 0);
    outTree->Branch("muonQ", &muonsQ);
    outTree->Branch("muonTrkIso", &muonsTrkIso);
    outTree->Branch("muonIsGLB", &muonIsGLB);
    outTree->Branch("muonMuNChi2", &muonMuNChi2);
    outTree->Branch("muonNMatchStn", &muonNMatchStn);
    outTree->Branch("muonNPixHits", &muonNPixHits);
    outTree->Branch("muonD0", &muonD0);
    outTree->Branch("muonDz", &muonDz);
    outTree->Branch("muonNTkLayers", &muonNTkLayers);
    outTree->Branch("muonNValidHits", &muonNValidHits);
    outTree->Branch("muonPassStdCuts", &muonPassStdCuts);
//  outTree->Branch("muonMatchBits", &muonMatchBits);

    // electrons
    outTree->Branch("electronP4", &electronsP4, 256000, 0);
    outTree->Branch("electronQ", &electronsQ);
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
//  outTree->Branch("electronMatchBits", &electronMatchBits);

    // gen-level particles
    outTree->Branch("genMuonP4", &genMuonsP4, 256000, 0);
    outTree->Branch("genMuonQ", &genMuonsQ);
    outTree->Branch("genElectronP4", &genElectronsP4, 256000, 0);
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
    muonIsGLB.clear(); muonPassStdCuts.clear();
    muonMuNChi2.clear(); muonD0.clear(); muonDz.clear();
    muonNMatchStn.clear(); muonNPixHits.clear(); muonNTkLayers.clear(); muonNValidHits.clear();
    electronIsConv.clear(); electronPassID.clear(); electronPassIso.clear();
    electronPassStdCuts.clear();
    electronCombIso.clear(); electronEnergyInv.clear();
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
//  if (!isData)
//      hTotalEvents->Fill(2);

//  /* Trigger selection */
//  bool passTrigger = false;
//  for (unsigned i = 0; i < triggerNames.size(); ++i) {
//      passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
//  }

//  if (!passTrigger && isData)
//      return kTRUE;

//  if (sync_print) {
//      cout << "trigger status: " << passTrigger << "\n" << endl;
//  }

//  hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
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

        // Perform selection
        if (muon->trkIso/muon->pt < 0.1) {
            if (muonP4.Pt() > 25 && fabs(muonP4.Eta()) < 2.1) {
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

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
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
            muonIsGLB.push_back(muon->typeBits & baconhep::kGlobal);
            muonMuNChi2.push_back(muon->muNchi2);
            muonNMatchStn.push_back(muon->nMatchStn);
            muonNPixHits.push_back(muon->nPixHits);
            muonD0.push_back(muon->d0);
            muonDz.push_back(muon->dz);
            muonNTkLayers.push_back(muon->nTkLayers);
            muonNValidHits.push_back(muon->nValidHits);
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

            int iEta = 0;
            float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
            for (unsigned i = 0; i < 8; ++i) {
                if (fabs(electron->scEta) > etaBins[i] && fabs(electron->scEta) < etaBins[i+1]) {
                    iEta = i;
                    break;
                }
            }
            electronCombIso.push_back(electron->chHadIso + std::max(0.,(double)electron->neuHadIso + electron->gammaIso - fInfo->rhoJet * cuts->EAEl[iEta]));
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
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
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
