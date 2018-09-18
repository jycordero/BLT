#include "AcceptanceAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

TLorentzVector GetP4Sum(vector<TLorentzVector> p4)
{
    TLorentzVector sum(0, 0, 0, 0);
    for (auto p4_ = p4.begin(); p4_ != p4.end(); ++p4_)
        sum += *p4_;
    return sum;
}

AcceptanceAnalyzer::AcceptanceAnalyzer() : BLTSelector()
{

}

AcceptanceAnalyzer::~AcceptanceAnalyzer()
{

}

void AcceptanceAnalyzer::Begin(TTree *tree)
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

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(), "RECREATE");
    outFile->cd();

    // Event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    // Acceptance counters
    outHistName = params->get_output_treename("GenEvents");
    hGenEvents = new TH1D(outHistName.c_str(), "GenEvents", nbins, xbins);
    hGenEvents->Sumw2();

    outHistName = params->get_output_treename("AccEvents");
    hAccEvents = new TH1D(outHistName.c_str(), "AccEvents", nbins, xbins);
    hAccEvents->Sumw2();

    ReportPostBegin();
}

void AcceptanceAnalyzer::Init(TTree *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrentFile = tree->GetCurrentFile();

    fInfo                    = 0;
    fGenEvtInfo              = 0;
    fGenParticleArr          = 0;
    fLHEWeightArr            = 0;

    fChain->SetBranchAddress("Info", &fInfo, &b_Info);
    fChain->SetBranchAddress("GenEvtInfo", &fGenEvtInfo, &b_GenEvtInfo);
    fChain->SetBranchAddress("GenParticle", &fGenParticleArr, &b_GenParticleArr);
    fChain->SetBranchAddress("LHEWeight", &fLHEWeightArr, &b_LHEWeightArr);
}

Bool_t AcceptanceAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches included above
    this->totalEvents++;
    hTotalEvents->Fill(1);

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0) { 
        std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;
    }

    const bool isData = (fInfo->runNum != 1);


    ///////////////////////
    // Generator objects //
    ///////////////////////
        
    if (!isData)
    {
        vector<TLorentzVector> elecs, muons;

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

//          if (abs(particle->pdgId) == 23 || abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 || abs(particle->pdgId) == 15)
//              cout << i << "\t" << particle->status << "\t" << particle->pdgId  << "\t" << particle->parent << "\t" << mom.Pt() << endl;


            // Find leptons from hard-scattering process
            if (particle->status == 3)      // seems safe at least for DY sample (pythia 6?)
            {
                TLorentzVector mom;
                mom.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);

                if (abs(particle->pdgId) == 11)
                    elecs.push_back(mom);
                else if (abs(particle->pdgId) == 13)
                    muons.push_back(mom);
            }
        }


        // Apply selection (flavor)
        vector<TLorentzVector> leptons;

        if (params->selection == "ee")
            leptons = elecs;
        else if (params->selection == "mumu")
            leptons = muons;
        sort(leptons.begin(), leptons.end(), P4SortCondition);

        unsigned nLeps = leptons.size();
        if (nLeps != 2)                     // specific to Drell-Yan...
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector lepSum = GetP4Sum(leptons);
        hGenEvents->Fill(lepSum.M());

//      if (lepSum.M() > 116)
//          cout << leptons[0].Pt() << ", " << leptons[1].Pt() << endl;


        // Apply cuts
        bool isAccepted = kTRUE;
        if (leptons[0].Pt() < PT1_MIN || fabs(leptons[0].Eta()) > ETA_MAX)
            isAccepted = kFALSE;
        if (leptons[1].Pt() < PT2_MIN || fabs(leptons[1].Eta()) > ETA_MAX)
            isAccepted = kFALSE;
//      for (unsigned i = 0; i < nLeps; i++)
//      {
//          if (leptons[i].Pt() < PT_MIN || fabs(leptons[i].Eta()) > ETA_MAX)
//          {
//              isAccepted = kFALSE;
//              break;
//          }
//      }

        if (isAccepted)
        {
            hAccEvents->Fill(lepSum.M());
            hTotalEvents->Fill(5);
        }
    }


    this->passedEvents++;
    return kTRUE;
}

void AcceptanceAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void AcceptanceAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void AcceptanceAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<AcceptanceAnalyzer> selector(new AcceptanceAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
