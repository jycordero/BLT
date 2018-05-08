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

//  muonsP4->SetClass("TLorentzVector");
//  electronsP4->SetClass("TLorentzVector");

    // event data
//  outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
//  outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("muonP4", &muonsP4, 256000, 0);
    outTree->Branch("muonIso", &muonsIso);
    outTree->Branch("muonQ", &muonsQ);
    outTree->Branch("electronP4", &electronsP4, 256000, 0);
    outTree->Branch("electronIso", &electronsIso);
    outTree->Branch("electronQ", &electronsQ);


    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{
    muonsP4ptr.Clear();
    electronsP4ptr.Clear();

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
//  triggerStatus = passTrigger;
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
//              muon->pt > 10 
//              && fabs(muon->eta) < 2.4
                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                   (muon->typeBits & baconhep::kGlobal) 
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
//  vector<TLorentzVector> veto_muons;
    vector<float> muons_iso;
    vector<int> muons_q;
//  vector<bool> muons_trigger;

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
//      if (muon->trkIso/muon->pt < 0.1) {

//          if (muonP4.Pt() > 20 && fabs(muonP4.Eta()) < 2.1) {
//              veto_muons.push_back(muonP4);

//              if (muonP4.Pt() > 25) {
                    muons.push_back(muonP4);
                    muons_iso.push_back(muon->trkIso);
                    muons_q.push_back(muon->q);

//                  // trigger matching
//                  bool triggered = false;
//                  for (unsigned i = 0; i < triggerNames.size(); ++i) {
//                      triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
//                  }
//                  muons_trigger.push_back(triggered);
//              }
//          }
//      }
    }
    std::sort(muons.begin(), muons.end(), P4SortCondition);

    if (sync_print) cout << endl;

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<int> electrons_q;
//  vector<bool> electrons_trigger;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
//              electron->pt > 20 
//              && fabs(electron->eta) < 2.5
                   particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
            electrons_iso.push_back(0.);
            electrons_q.push_back(electron->q);

//          // trigger matching
//          bool triggered = false;
//          for (unsigned i = 0; i < triggerNames.size(); ++i) {
//              triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
//          }
//          electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);


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
//      if (nBJets >= 1 && nJets >= 1 && met < 40 && muons.size() >= 2) {

//          jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
//          bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);

//          TLorentzVector dijet = bjetP4 + jetP4;
//          TLorentzVector dimuon = muons[0] + muons[1];
//          cout << "phi_mumu, phi_jj, dphi_mumujj" << endl;
//          cout << dimuon.Phi() << ", " << dijet.Phi() << ", " << fabs(dimuon.DeltaPhi(dijet)) << endl;
//      }
//      cout << "STOP!" << endl;
//      return kTRUE;

    } else if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4;
//      unsigned muonTwoIndex = 1;
//      muonOneP4 = muons[0];
//      for (unsigned i = 1; i < muons.size(); ++i) {
//          if (muons_q[0] != muons_q[i]) {
//              muonTwoP4 = muons[i];
//              muonTwoIndex = i;
//              break;
//          }
//      }

//      if (
//              !(muonOneP4.Pt() > 25 && fabs(muonOneP4.Eta()) < 2.1) 
//              || !(muonTwoP4.Pt() > 25 && fabs(muonTwoP4.Eta()) < 2.1)
//         )
//          return kTRUE;
//      hTotalEvents->Fill(6);

//      TLorentzVector dimuon;
//      dimuon = muonOneP4 + muonTwoP4;
//      if (sync_print) {
//          cout << dimuon.M() << endl;
//      }

//      if (dimuon.M() < 12. || dimuon.M() > 110.)
//          return kTRUE;
//      hTotalEvents->Fill(7);

//      if (nBJets == 0 || (nBJets + nJets < 2 && nFwdJets == 0)) 
//          return kTRUE;
//      hTotalEvents->Fill(8);

        for (unsigned i = 0; i < muons.size(); ++i)
            new(muonsP4ptr[i]) TLorentzVector(muons[i]);
        muonsIso = muons_iso;
        muonsQ = muons_q;

        for (unsigned i = 0; i < electrons.size(); ++i)
            new(electronsP4ptr[i]) TLorentzVector(electrons[i]);
        electronsIso = electrons_iso;
        electronsQ = electrons_q;

//      if (!isData) {
//          eventWeight *= weights->GetMuonIDEff(muonOneP4);
//          eventWeight *= weights->GetMuonISOEff(muonOneP4);
//          eventWeight *= weights->GetMuonIDEff(muonTwoP4);
//          eventWeight *= weights->GetMuonISOEff(muonTwoP4);

//          // trigger weight
//          pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
//          pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
//          eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
//      }


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
