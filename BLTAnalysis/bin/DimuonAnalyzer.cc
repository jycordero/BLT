#include "DimuonAnalyzer.hh"

//
// See header file for class documentation
//


using namespace baconhep;

DimuonAnalyzer::DimuonAnalyzer() : BLTSelector()
{

}

DimuonAnalyzer::~DimuonAnalyzer()
{

}

void DimuonAnalyzer::Begin(TTree *tree)
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
    triggerSelector.reset(new TriggerSelector());

    // Prepare the output tree
    outFileName = params->get_output_filename("demoFile");
    outTreeName = params->get_output_treename("demoTree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("muonOne", &muonOne);
    outTree->Branch("muonTwo", &muonTwo);
    outTree->Branch("dimuon", &dimuon);
    outTree->Branch("genMuonOne", &genMuonOne);
    outTree->Branch("genMuonTwo", &genMuonTwo);
    outTree->Branch("genZ", &genZ);
    outTree->Branch("met", &met);
    outTree->Branch("met_phi", &met_phi);

    ReportPostBegin();
}

Bool_t DimuonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    //if (entry%1==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);
    triggerSelector->SetRealData(isRealData);

    bool printEvent = false;

    //////////////////
    // GenParticles //
    //////////////////

    if (printEvent) {
        if (!isRealData) {
            for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
                const TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
                assert(particle);
                std::cout << "GenParticle " << i << ": " << particle << std::endl;
            }
        }
    }

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
        //particleSelector->SetPV(TVector3());
    }
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);


    ///////////////////
    // Select objects//
    ///////////////////

    /* MUONS */
    std::vector<TMuon*> muons;

    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
            muon->pt > 20 
            && std::abs(muon->eta) < 2.4
            && particleSelector->PassMuonID(muon, cuts->tightMuID)
            && particleSelector->PassMuonIso(muon, cuts->tightMuIso)
           )
            muons.push_back(muon);
    }

    std::sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* JETS */
    std::vector<TJet*> jets;
    std::vector<TJet*> bjets;
    for (int i=0; i<fAK4CHSArr->GetEntries(); i++) {
        TJet* jet = (TJet*) fAK4CHSArr->At(i);
        assert(jet);
        if (
                jet->pt > 30
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (particleSelector->PassJetID(jet, cuts->bJetID))
                bjets.push_back(jet);
            else 
                jets.push_back(jet);
        }
    }

    /* MET */
    TMET* pfMET = new TMET();
    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;

    if (printEvent) {
        std::cout << "MET " << "(PF)" << ": " << pfMET << std::endl;
    }

    /* Apply dimuon selection */
    TLorentzVector tmp_muonOne, tmp_muonTwo;
    int            idx_muonOne, idx_muonTwo;

    bool found_dimuon = particleSelector->FindGoodDiMuons(muons, tmp_muonOne, tmp_muonTwo, idx_muonOne, idx_muonTwo);

    if (!found_dimuon)
        return kFALSE;

    /* Gen matching */
    TLorentzVector tmp_genZ, tmp_genMuonOne, tmp_genMuonTwo;
    int            idx_genZ, idx_genMuonOne, idx_genMuonTwo;

    bool found_genZ = particleSelector->FindGenZToLL(fGenParticleArr, tmp_genZ, tmp_genMuonOne, tmp_genMuonTwo, idx_genZ, idx_genMuonOne, idx_genMuonTwo);

    if (!found_genZ)
        return kFALSE;

    if (std::abs(((TGenParticle *)fGenParticleArr->At(idx_genMuonOne))->pdgId) != MUON_PDGID ||
        std::abs(((TGenParticle *)fGenParticleArr->At(idx_genMuonTwo))->pdgId) != MUON_PDGID)
        return kFALSE;


    //////////
    // Fill //
    //////////

    muonOne = tmp_muonOne;
    muonTwo = tmp_muonTwo;
    dimuon  = muonOne + muonTwo;

    met = pfMET->pt;
    met_phi = pfMET->phi;

    if (tmp_genMuonOne.Pt() > tmp_genMuonTwo.Pt()) {
        genMuonOne = tmp_genMuonOne;
        genMuonTwo = tmp_genMuonTwo;
    } else {
        genMuonOne = tmp_genMuonTwo;
        genMuonTwo = tmp_genMuonOne;
    }
    genZ = genMuonOne + genMuonTwo;

    outTree->Fill();
    this->passedEvents++;

    delete pfMET;

    return kTRUE;
}

void DimuonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DimuonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DimuonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<DimuonAnalyzer> selector(new DimuonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
