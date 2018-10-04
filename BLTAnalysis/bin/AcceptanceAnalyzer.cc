#include "AcceptanceAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


AcceptanceAnalyzer::AcceptanceAnalyzer() : BLTSelector()
{

}

AcceptanceAnalyzer::~AcceptanceAnalyzer()
{

}

void AcceptanceAnalyzer::Begin(TTree *tree)
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


    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();


    //--- HISTOGRAMS ---//

    // Event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    // Acceptance counters
    outHistName = params->get_output_treename("GenEvents");
    hGenEvents = new TH1D(outHistName.c_str(), "GenEvents", 10, 0.5, 10.5);
    hGenEvents->Sumw2();

    outHistName = params->get_output_treename("AccEvents");
    hAccEvents = new TH1D(outHistName.c_str(), "AccEvents", 10, 0.5, 10.5);
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

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    const bool isData = (fInfo->runNum != 1);
   



    /////////////////////
    //  GEN PARTICLES  //
    /////////////////////


    genWeight = 1;

    if (!isData)
    {
        // Save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (genWeight < 0)
            hTotalEvents->Fill(10);
        hAcceptedEvents->Fill(1, genWeight);



        //--- PARTICLE LOOP ---//

        vector<TLorentzVector> leptons, elecs, muons;   // Momenta
        vector<int> leptons_q, elecs_q, muons_q;        // Charges

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);


            // Find leptons from hard-scattering process
            if (
                    (abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13)
                    &&  (particle->status == 23 || particle->status == 1 || particle->status == 2)
                    &&  particle->parent != -2
                )
            {
                // Determine if lepton comes from a Z
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                int origin = abs(mother->pdgId);

                // (status 23 particle seems to be guaranteed as a hard scatter product
                //  and tends to be missing mother info)
                if (origin == 23 || particle->status == 23)
                {
                    int q = copysign(particle->pdgId);
                    TLorentzVector lep;

                    if (abs(particle->pdgId) == 11)
                    {
                        eCount++;
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        elecs.push_back(lep);       elecs_q.push_back(q);
                        leptons.push_back(lep);     leptons_q.push_back(q);

                    }
                    else if (abs(particle->pdgId) == 13)
                    {
                        muCount++;
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        muons.push_back(lep);       muons_q.push_back(q);
                        leptons.push_back(lep);     leptons_q.push_back(q);
                    }
                }
            }
        }



        //--- SELECTION ---//

        unsigned muCount = muons.size(), eCount = elecs.size(), lepCount = leptons.size();
        TLorentzVector muSum = GetP4Sum(muons), eSum = GetP4Sum(elecs), lepSum = GetP4Sum(leptons);

        // Total charge for each lepton type
        int eCharge = accumulate(elec_q.begin(), elec_q.end(), 0);
        int muCharge = accumulate(muon_q.begin(), muon_q.end(), 0);



        //--- PHASE SPACE ---//

        // Mass window
        if (lepSum.M() < M_MIN || lepSum.M() > M_MAX)
            return kTRUE;


        // Charge requirement
        if (eCharge != 0 || muCharge != 0)
            return kTRUE;


        // Sort events by decay channel
        unsigned idx = 0;

        if (muCount == 2 && eCount == 0)            // mumu = 3
            idx = 3;

        else if (muCount == 0 && eCount == 2)       // ee   = 4
            idx = 4;

        else if (muCount == 4 && eCount == 0)       // 4m   = 6
            idx = 6;

        else if (muCount == 2 && eCount == 2        // 2m2e = 7
                && muSum.M() > eSum.M())
            idx = 7;

        else if (muCount == 2 && eCount == 2        // 2e2m = 8
                && muSum.M() < eSum.M())
            idx = 8;

        else if (muCount == 0 && eCount == 4)       // 4e   = 9
            idx = 9;

        else
            return kTRUE;


        // Dilepton mass requirement
        if (muCount == 2 && eCount == 2)        // Mixed flavor
        {
            if (eSum.M() < MLL_MIN)
                return kTRUE;

            if (muSum.M() < MLL_MIN)
                return kTRUE;
        }
        else if (muCount == 4 || eCount == 4)   // Single flavor
        {
            for (unsigned j = 1; j < 4; j++)
            {
                for (unsigned i = 0; i < j; i++)
                {
                    if (lep_q[i] != lep_q[j])
                    {
                        TLorentzVector dilep = leptons[i] + leptons[j];
                        if (dilep.M() < MLL_MIN)
                            return kTRUE;
                    }
                }
            }
        }


        // Remaining events must be in phase space region
        hGenEvents->Fill(idx, genWeight);



        //--- FIDUCIAL REGION ---//


        sort(leptons.begin(), leptons.end(), P4SortCondition); 

        // Eta
 
        for (unsigned i = 0; i < lepCount; i++)
        {
            if (fabs(leptons[i].Eta()) > ETA_MAX)
                return kTRUE;
        }


        // Lepton Pt

        if (leptons[0].Pt() < PT1_MIN)
            return kTRUE;

        if (leptons[1].Pt() < PT2_MIN)
            return kTRUE;
        
        if (lepCount == 4)
        {
            if (leptons[3].Pt() < PT_MIN || leptons[2].Pt() < PT_MIN)
                return kTRUE;
        }


        // Remaining events must be in fiducial region
        hAccEvents->Fill(idx, genWeight);


        if (kTRUE)
        {
            cout << "Index \tStatus\tID    \tParent" << endl;
            for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
            {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
                cout << i << "\t" << particle->status << "\t" << particle->pdgId << "\t" << particle->parent << endl;
            }           
            cout << endl << endl;
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

TLorentzVector GetP4Sum(vector<TLorentzVector> p4)
{
    TLorentzVector sum(0, 0, 0, 0);
    for (auto p4_ = p4.begin(); p4_ != p4.end(); ++p4_)
        sum += *p4_;
    return sum;
}
