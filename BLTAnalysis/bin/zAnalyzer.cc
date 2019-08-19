#include "zAnalyzer.h"
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

zAnalyzer::zAnalyzer() : BLTSelector()
{

}

zAnalyzer::~zAnalyzer()
{

}

void zAnalyzer::Begin(TTree *tree)
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
        
    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");
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
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    // weights
    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("elIDWeightOne", &elIDWeightOne);
    outTree->Branch("elIDWeightTwo", &elIDWeightTwo);
    outTree->Branch("elTrigWeightOne", &elTrigWeightOne);
    outTree->Branch("elTrigWeightTwo", &elTrigWeightTwo);
    outTree->Branch("muonIDWeightOne", &muonIDWeightOne);
    outTree->Branch("muonIDWeightTwo", &muonIDWeightTwo);
    outTree->Branch("muonTrigWeightOne", &muonTrigWeightOne);
    outTree->Branch("muonTrigWeightTwo", &muonTrigWeightTwo);
    outTree->Branch("photonIDWeight", &photonIDWeight);
    
    outTree->Branch("Sgen", &Sgen);
    outTree->Branch("SgenAccep", &SgenAccep);

    // leptons
    outTree->Branch("leptonOnePt", &leptonOnePt);
    outTree->Branch("leptonOneEta", &leptonOneEta);
    outTree->Branch("leptonOnePhi", &leptonOnePhi);
    outTree->Branch("leptonOnePtKin", &leptonOnePtKin);
    outTree->Branch("leptonOnePtKinJames", &leptonOnePtKinJames);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
    outTree->Branch("leptonOneECALDriven", &leptonOneECALDriven);
    outTree->Branch("leptonOneCharge", &leptonOneCharge);
    outTree->Branch("leptonOneTag", &leptonOneTag);

    outTree->Branch("leptonTwoPt", &leptonTwoPt);
    outTree->Branch("leptonTwoEta", &leptonTwoEta);
    outTree->Branch("leptonTwoPhi", &leptonTwoPhi);
    outTree->Branch("leptonTwoPtKin", &leptonTwoPtKin);
    outTree->Branch("leptonTwoPtKinJames", &leptonTwoPtKinJames);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
    outTree->Branch("leptonTwoECALDriven", &leptonTwoECALDriven);
    outTree->Branch("leptonTwoCharge", &leptonTwoCharge);
    outTree->Branch("leptonTwoTag", &leptonTwoTag);


    outTree->Branch("lepTagR9"    , &lepTagR9);
    //outTree->Branch("lepTagS4"    , &lepTagS4);
    outTree->Branch("lepTagSieie" , &lepTagSieie); 
    outTree->Branch("lepTagHoverE", &lepTagHoverE);
    outTree->Branch("lepTagIneu"  , &lepTagIneu);
    outTree->Branch("lepTagIph"   , &lepTagIph);
    outTree->Branch("lepTagIch"   , &lepTagIch);

    outTree->Branch("lepProbeR9"    , &lepProbeR9);
    //outTree->Branch("lepProbeS4"    , &lepProbeS4);
    outTree->Branch("lepProbeSieie" , &lepProbeSieie); 
    outTree->Branch("lepProbeHoverE", &lepProbeHoverE);
    outTree->Branch("lepProbeIneu"  , &lepProbeIneu);
    outTree->Branch("lepProbeIph"   , &lepProbeIph);
    outTree->Branch("lepProbeIch"   , &lepProbeIch);
 
    outTree->Branch("tauDecayMode", &tauDecayMode);
    outTree->Branch("tauMVA", &tauMVA);

    outTree->Branch("isLeptonTag", &isLeptonTag);
    outTree->Branch("isDijetTag", &isDijetTag);
    outTree->Branch("isTightDijetTag", &isTightDijetTag);

    outTree->Branch("vetoDY", &vetoDY);

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
    outTree->Branch("dileptonPt", &dileptonPt);
    outTree->Branch("dileptonEta", &dileptonEta);
    outTree->Branch("dileptonPhi", &dileptonPhi);
    outTree->Branch("dileptonM", &dileptonM);
    outTree->Branch("dileptonDEta", &dileptonDEta);
    outTree->Branch("dileptonDPhi", &dileptonDPhi);
    outTree->Branch("dileptonDR", &dileptonDR);
    outTree->Branch("dileptonMKin", &dileptonMKin);
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

    ReportPostBegin();
}

Bool_t zAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);
   
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    
    genWeight = 1;
    if (!isData) {
        if (fGenEvtInfo->weight < 0) {
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
            hTotalEvents->Fill(maxBin);
        }
    }

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
            ) return kTRUE;

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

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

    } else {
        nPartons = 0;
    }


    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

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
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    if (sync_print_precut) {
        cout << "pvx, pvy, pvz, ndof" << endl;
        cout << thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << thePV->ndof << endl;
    }

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
    }
    
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        //assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->calibPt, electron->eta, electron->phi, ELE_MASS);
        
        if (sync_print_precut) {
           cout << "electron pt, calibpt, eta, sc_eta, phi, d0, dz, sip3d, iso, mva, pass_mva" << endl;
           //cout << electronP4.Pt() << "," << electronP4.Eta() << "," << 
           cout << electron->pt << "," << electron->calibPt << "," <<  electron->eta << ", " << 
                   electron->scEta << ", " << electron->phi << "," << electron->d0 << "," << 
                   electron->dz << "," << electron->sip3d << "," << 
                   GetElectronIsolation(electron, fInfo->rhoJet) << ", " << electron->mva << ", " << 
                   particleSelector->PassElectronMVA(electron, cuts->hzzMVAID) <<endl;
        } 

        if (
                electron->calibPt > 7
                && fabs(electron->scEta) < 2.5
                //&& particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && GetElectronIsolation(electron, fInfo->rhoJet)/electronP4.Pt() < 0.35
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
                && fabs(electron->sip3d) < 4.0 
           ) {
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

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

    /* PHOTONS */
    vector <TPhoton*> photons;
    vector<TLorentzVector> veto_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photon->calibPt, photon->eta, photon->phi, 0.);
    	cout << " ---------Photon PT :: " << photonP4.Pt() << endl;
        if (sync_print_precut) {
            cout << "photon_pt, photon_calibpt, photon_eta, photon_sc_eta, photon_phi, photon_mva, pass_electron_veto" << endl;
            cout << photon->pt << ", " << photon->calibPt << ", " << photon->eta << ", " << photon->scEta << ", " << photon->phi << ", " << photon->mva 
                 << ", " << photon->passElectronVeto << endl;
        }

        if (
                // ID conditions
                photon->calibPt > 10
                && fabs(photon->scEta) < 2.5 
                && (fabs(photon->scEta) <= 1.4442 || fabs(photon->scEta) >= 1.566)
                //&& particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
                && particleSelector->PassPhotonID(photon, cuts->mediumPhID)
                && photon->passElectronVeto
            ) {
            photons.push_back(photon);
            veto_photons.push_back(photonP4);
        }
    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);

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
                && fabs(jet->eta) < 4.7
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

    
   if (params->selection == "ee") {
        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

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

        if (leptonOneP4.Pt() <= 25.0)
            return kTRUE;

        if (leptonTwoP4.Pt() <= 15.0)
            return kTRUE;
        
	leptonOneP4.SetPtEtaPhiM(electrons[electronOneIndex]->calibPt, electrons[electronOneIndex]->eta, electrons[electronOneIndex]->phi, ELE_MASS);
	leptonTwoP4.SetPtEtaPhiM(electrons[electronTwoIndex]->calibPt, electrons[electronTwoIndex]->eta, electrons[electronTwoIndex]->phi, ELE_MASS);


        TLorentzVector dielectron = leptonOneP4 + leptonTwoP4;
        if (electrons[electronOneIndex]->q == electrons[electronTwoIndex]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        
        leptonOneIso    = GetElectronIsolation(electrons[electronOneIndex], fInfo->rhoJet);
        leptonOneFlavor = electrons[electronOneIndex]->q*11;
        leptonOneDZ     = electrons[electronOneIndex]->dz;
        leptonOneD0     = electrons[electronOneIndex]->d0;
            
        leptonTwoIso    = GetElectronIsolation(electrons[electronTwoIndex], fInfo->rhoJet);
        leptonTwoFlavor = electrons[electronTwoIndex]->q*11;
        leptonTwoDZ     = electrons[electronTwoIndex]->dz;
        leptonTwoD0     = electrons[electronTwoIndex]->d0;

	lepTagR9     = electrons[electronOneIndex]->r9;
	lepTagSieie  = electrons[electronOneIndex]->sieie;
	lepTagHoverE = electrons[electronOneIndex]->hovere;
	lepTagIneu   = electrons[electronOneIndex]->neuHadIso;
	lepTagIph    = electrons[electronOneIndex]->gammaIso; 	
	lepTagIch    = electrons[electronOneIndex]->chHadIso;       

	lepProbeR9     = electrons[electronTwoIndex]->r9;
	lepProbeSieie  = electrons[electronTwoIndex]->sieie;
	lepProbeHoverE = electrons[electronTwoIndex]->hovere;
	lepProbeIneu   = electrons[electronTwoIndex]->neuHadIso;
	lepProbeIph    = electrons[electronTwoIndex]->gammaIso; 	
	lepProbeIch    = electrons[electronTwoIndex]->chHadIso;       
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
	
    } 
    


    ///////////////////
    // Fill jet info //
    ///////////////////
    

    outTree->Fill();
    this->passedEvents++;
    Sgen = 0;
    SgenAccep = 0;

    if (sync_print_precut) {
        cout << "event should have been filled" << endl;
        }
    return kTRUE;
}

void zAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<zAnalyzer> selector(new zAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float zAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float zAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
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

void zAnalyzer::EvalMuonEnergyResolution(std::map<string, float> mva_input_floats, std::map<string, int> mva_input_ints, float &mean, float &sigma, float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;

    // one-time MVA initialization: get trainings
    if (!mvaMean) {
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");
         // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;
         RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_muons"));
        if (!ws) FATAL("TFile::Get() failed");
        // current working directory
        invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");
        mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");
     } // one-time initialization
     // load necessary tree branches
    float rho                      = mva_input_floats["rho"];
    float muEnergy                 = mva_input_floats["muEnergy"];
    float muEta                    = mva_input_floats["muEta"];
    float muChi2NDF                = mva_input_floats["muTkChi2"];
    Int_t muNumberOfValidTrkLayers = mva_input_ints["muNumberOfValidTrkLayers"];
    Int_t muNumberOfValidPixelHits = mva_input_ints["muNumberOfValidPixelHits"];
    Int_t muNumberOfValidMuonHits  = mva_input_ints["muNumberOfValidMuonHits"];
    Int_t muStations               = mva_input_ints["muStations"];
    float muPFIsoR04_CH            = mva_input_floats["muPFIsoR04_CH"];
    float muPFIsoR04_NH            = mva_input_floats["muPFIsoR04_NH"];
    float muPFIsoR04_Pho           = mva_input_floats["muPFIsoR04_Pho"];
    float muPFIsoR04_PU            = mva_input_floats["muPFIsoR04_PU"];
     // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = muEnergy;
    *invar[2] = muEta;
    *invar[3] = muChi2NDF;
    *invar[4] = muNumberOfValidTrkLayers;
    *invar[5] = muNumberOfValidPixelHits;
    *invar[6] = muNumberOfValidMuonHits;
    *invar[7] = muStations;
    *invar[8] = muPFIsoR04_CH;
    *invar[9] = muPFIsoR04_NH;
    *invar[10] = muPFIsoR04_Pho;
    *invar[11] = muPFIsoR04_PU;
    mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();
     // NOTE: negative = infinite; powers were fixed at the training level
    powerL = -1;
    powerR = -1;
    //f.Close();
 }
 void zAnalyzer::EvalElectronEnergyResolution(std::map<string, float> mva_inputs, float &mean, float &sigma, 
                                                    float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;
    static RooAbsReal* mvaPowerL;
    static RooAbsReal* mvaPowerR;
     // one-time MVA initialization: get trainings
    if (!mvaMean) {
        // current working directory
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");
         // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;
         RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_electrons"));
        if (!ws) FATAL("TFile::Get() failed");
         invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");
        invar[12] = ws->var("var13");
        invar[13] = ws->var("var14");
        invar[14] = ws->var("var15");
         mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");
        mvaPowerL = ws->function("limPowerL");
        mvaPowerR = ws->function("limPowerR");
     } // one-time initialization
     // load necessary tree branches
    float rho = mva_inputs["rho"];
    float elEnergy = mva_inputs["elEnergy"];
    float elScEta = mva_inputs["elScEta"];
    float elScPhi = mva_inputs["elScPhi"];
    float elR9 = mva_inputs["elR9"];
    float elE1x5OverE = mva_inputs["elE1x5OverE"];
    float elE2x5OverE = mva_inputs["elE2x5OverE"];
    float elE5x5OverE = mva_inputs["elE5x5OverE"];
    float elFBrem = mva_inputs["elFBrem"];
    float elHOverE = mva_inputs["elHOverE"];
    float elSigmaIEtaIEta = mva_inputs["elSigmaIEtaIEta"];
    float elPFIsoR04_CH = mva_inputs["elPFIsoR04_CH"];
    float elPFIsoR04_NH = mva_inputs["elPFIsoR04_NH"];
    float elPFIsoR04_Pho = mva_inputs["elPFIsoR04_Pho"];
    float elPFIsoR04_PU = mva_inputs["elPFIsoR04_PU"];
     // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = elEnergy;
    *invar[2] = elScEta;
    *invar[3] = elScPhi;
    *invar[4] = elR9;
    *invar[5] = elE1x5OverE;
    *invar[6] = elE2x5OverE;
    *invar[7] = elE5x5OverE;
    *invar[8] = elFBrem;
    *invar[9] = elHOverE;
    *invar[10] = elSigmaIEtaIEta;
    *invar[11] = elPFIsoR04_CH;
    *invar[12] = elPFIsoR04_NH;
    *invar[13] = elPFIsoR04_Pho;
    *invar[14] = elPFIsoR04_PU;
     mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();
    powerL = mvaPowerL->getVal();
    powerR = mvaPowerR->getVal();
}

 void zAnalyzer::find_optimized(double* p, double &e1, double& e2)
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
