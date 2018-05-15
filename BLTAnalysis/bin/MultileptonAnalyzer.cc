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
//  std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2";
//  trigger.reset(new baconhep::TTrigger(trigfilename));

//  if (params->selection == "mumu" || params->selection == "emu") {
//      triggerNames.push_back("HLT_IsoMu24_eta2p1_v*");
//  } else if (params->selection == "ee") {
//      triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
//  }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false));

    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

    // muon momentum corrections
    muonCorr = new rochcor2012();

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
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // muons
//  outTree->Branch("muonP4", &muonsP4, 256000, 0);
//  outTree->Branch("muonQ", &muonsQ);
//  outTree->Branch("muonTrkIso", &muonsTrkIso);
//  outTree->Branch("muonIsGLB", &muonIsGLB);
//  outTree->Branch("muonMuNChi2", &muonMuNChi2);
//  outTree->Branch("muonNMatchStn", &muonNMatchStn);
//  outTree->Branch("muonNPixHits", &muonNPixHits);
//  outTree->Branch("muonD0", &muonD0);
//  outTree->Branch("muonDz", &muonDz);
//  outTree->Branch("muonNTkLayers", &muonNTkLayers);
//  outTree->Branch("muonNValidHits", &muonNValidHits);
//  outTree->Branch("muonPassStdCuts", &muonPassStdCuts);
//  outTree->Branch("muonMatchBits", &muonMatchBits);

    // electrons
//  outTree->Branch("electronP4", &electronsP4, 256000, 0);
//  outTree->Branch("electronQ", &electronsQ);
//  outTree->Branch("electronTrkIso", &electronsTrkIso);
//  outTree->Branch("electronCombIso", &electronCombIso);
//  outTree->Branch("electronEnergyInv", &electronEnergyInv);
//  outTree->Branch("electronScEta", &electronScEta);
//  outTree->Branch("electronD0", &electronD0);
//  outTree->Branch("electronDz", &electronDz);
//  outTree->Branch("electronSieie", &electronSieie);
//  outTree->Branch("electronHOverE", &electronHOverE);
//  outTree->Branch("electronDEtaIn", &electronDEtaIn);
//  outTree->Branch("electronDPhiIn", &electronDPhiIn);
//  outTree->Branch("electronNMissHits", &electronNMissHits);
//  outTree->Branch("electronIsConv", &electronIsConv);
//  outTree->Branch("electronPassID", &electronPassID);
//  outTree->Branch("electronPassIso", &electronPassIso);
//  outTree->Branch("electronPassStdCuts", &electronPassStdCuts);
//  outTree->Branch("electronMatchBits", &electronMatchBits);

    // gen-level particles
//  outTree->Branch("genMuonP4", &genMuonsP4, 256000, 0);
//  outTree->Branch("genMuonQ", &genMuonsQ);
//  outTree->Branch("genElectronP4", &genElectronsP4, 256000, 0);
//  outTree->Branch("genElectronQ", &genElectronsQ);
//  outTree->Branch("genIntermID", &genIntermID);
//  outTree->Branch("genIntermMass", &genIntermMass);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
//  outTree->Branch("nLeptons", &nLeptons);
//  outTree->Branch("nStdMuons", &nStdMuons);
//  outTree->Branch("nStdElectrons", &nStdElectrons);
//  outTree->Branch("nStdLeptons", &nStdLeptons);
//  outTree->Branch("nGenMuons", &nGenMuons);
//  outTree->Branch("nGenElectrons", &nGenElectrons);
//  outTree->Branch("nGenLeptons", &nGenLeptons);

    // Event counter
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

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0) { 
        std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;
    }

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    // Apply lumi mask
    if (isRealData) {
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
//  if (!passTrigger)
//      return kTRUE;

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
    if (!isRealData) {
        nPU = fInfo->nPUmean;
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    if (!isRealData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            
            //cout << particle->status << ", "
            //     << particle->pdgId  << ", "
            //     << particle->parent
            //     << endl;

            if (particle->status == 3 && (abs(particle->pdgId) < 6 || particle->pdgId == 21)) {
                ++count;
            }
        }
        nPartons = count-4; // This is saved for reweighting inclusive DY and combining it with parton binned DY
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
    vector<TLorentzVector> veto_muons;
    vector<TMuon*> muons;
    for (unsigned i = 0; i < tmp_muons.size(); i++) {
        TMuon* muon = tmp_muons[i];

        TLorentzVector muonP4;
        copy_p4(tmp_muons[i], MUON_MASS, muonP4);

        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i+1; j < tmp_muons.size(); j++) {
            TLorentzVector muon_j;
            copy_p4(tmp_muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3) {
                muon->trkIso03 = max(0., muon->trkIso03 - muon_j.Pt());
                tmp_muons[j]->trkIso03 = max(0., tmp_muons[j]->trkIso03 - muonP4.Pt());
            }
        }

        // Apply rochester muon momentum corrections
        float qter = 1.;
        if (isRealData) {
            muonCorr->momcor_data(muonP4, muon->q, 0, qter);
        } else {
            muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
        }

        // Fill containers
        if (muon->trkIso03/muonP4.Pt() < 0.1) {
            if (muonP4.Pt() > 20)
                veto_muons.push_back(muonP4);

            if (muonP4.Pt() > 5) {
                muon->pt  = muonP4.Pt();
                muon->eta = muonP4.Eta();
                muon->phi = muonP4.Phi();
                muons.push_back(muon);
                ++nMuons;
            }
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    std::vector<TElectron*> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
            ++nElectrons;

            // trigger matching
            //bool triggered = false;
            //for (unsigned i = 0; i < triggerNames.size(); ++i) {
            //    triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            //}
        }
    }

    std::sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);


    // MET //
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    // Apply analysis selections //
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4;
//      unsigned muonTwoIndex = 1;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        if (muons.size() == 2) {
            muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        } else if (muons.size() > 2) {
            if (muons[0]->q != muons[1]->q) {
                muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
            } else if (muons[0]->q != muons[2]->q) {
//              muonTwoIndex = 2;
                muonTwoP4.SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[2]->phi, MUON_MASS);
            }
        }

        if (
            muonOneP4.Pt() < 25 
            || fabs(muonOneP4.Eta()) > 2.1
            || muonTwoP4.Pt() < 25 
            || fabs(muonTwoP4.Eta()) > 2.1
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muonOneP4 + muonTwoP4;
        if (dimuon.M() < 12 || dimuon.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(7);

//      leptonOneP4      = muonOneP4;
//      leptonOneIso     = muons[0]->trkIso03;
//      leptonOneQ       = muons[0]->q;
//      leptonOneFlavor  = 13;

//      leptonTwoP4      = muonTwoP4;
//      leptonTwoIso     = muons[muonTwoIndex]->trkIso03;
//      leptonTwoQ       = muons[muonTwoIndex]->q;
//      leptonTwoFlavor  = 13;

//      // trigger matching
//      for (unsigned i = 0; i < triggerNames.size(); ++i) {
//          leptonOneTrigger |= trigger->passObj(triggerNames[i], 1, muons[0]->hltMatchBits);
//          leptonTwoTrigger |= trigger->passObj(triggerNames[i], 1, muons[muonTwoIndex]->hltMatchBits);
//      }

//      if (!isRealData) {
//          eventWeight *= weights->GetMuonRecoEff(muonOneP4);
//          eventWeight *= weights->GetMuonRecoEff(muonTwoP4);

//          // trigger weight
//          pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonOneP4);
//          pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonTwoP4);
//          eventWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
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
