#include "AcceptanceAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;


bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
TLorentzVector GetP4Sum(const vector<TLorentzVector> &p4);


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
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");



    //--- BRANCHES ---//

    // Event
    outTree->Branch(    "runNumber",            &runNumber);
    outTree->Branch(    "evtNumber",            &evtNumber,             "eventNumber/l");
    outTree->Branch(    "lumiSection",          &lumiSection);
    outTree->Branch(    "genWeight",            &genWeight);
    outTree->Branch(    "isFiducial",           &isFiducial);

    // Counters
    outTree->Branch(    "nGenMuons",            &nGenMuons);
    outTree->Branch(    "nGenElectrons",        &nGenElectrons);
    outTree->Branch(    "nGenLeptons",          &nGenLeptons);

    // Gen muons
    outTree->Branch(    "genMuonP4",            &genMuonsP4,            32000,  1);
    outTree->Branch(    "genMuonQ",             &genMuonsQ);
    outTree->Branch(    "genMuonStatus",        &genMuonStatus);

    // Gen electrons
    outTree->Branch(    "genElectronP4",        &genElectronsP4,        32000,  1);
    outTree->Branch(    "genElectronQ",         &genElectronsQ);
    outTree->Branch(    "genElectronStatus",    &genElectronStatus);



    //--- HISTOGRAMS ---//

    // Total event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    // Acceptance counters
    outHistName = params->get_output_treename("PhaseSpaceEvents");
    hPhaseSpaceEvents = new TH1D(outHistName.c_str(), "PhaseSpaceEvents", 10, 0.5, 10.5);

    outHistName = params->get_output_treename("FiducialEvents");
    hFiducialEvents = new TH1D(outHistName.c_str(), "FiducialEvents", 10, 0.5, 10.5);


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


    //--- CLEAR CONTAINERS ---//
    
    isFiducial = kTRUE;     // innocent until proven guilty...

    nGenMuons = 0;                      nGenElectrons = 0;                  nGenLeptons = 0; 
    genMuonsP4ptr.Delete();             genMuonsQ.clear();                  genMuonStatus.clear();
    genElectronsP4ptr.Delete();         genElectronsQ.clear();              genElectronStatus.clear();



    /////////////////
    //    START    //
    /////////////////


    GetEntry(entry, 1);  // load all branches included above
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    const bool isData = (fInfo->runNum != 1);
    if (isData)
        return kTRUE;



    //--- PARTICLE LOOP ---//

    genWeight   = 1;
    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;


    if (!isData)
    {
        // Save gen weight for amc@nlo samples
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (genWeight < 0)
            hTotalEvents->Fill(10);




        /////////////////////
        //  GEN PARTICLES  //
        /////////////////////


        hPhaseSpaceEvents->Fill(1, genWeight);

        vector<TLorentzVector> genLeps, genMuons, genElecs;     // Momenta
        vector<int> genLepsQ;                                   // Charges
        int elecTotalQ = 0, muonTotalQ = 0;


  
        //--- HARD PROCESS ---//

        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);


            // Find leptons from hard-scattering process
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13)
                &&  (particle->status == 23 || particle->status == 1 || particle->status == 2)
                &&  particle->parent != -2)
            {
                // Determine if lepton comes from a Z
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                int origin = abs(mother->pdgId);

                // (status 23 particle seems to be guaranteed as a hard scatter product
                //  and tends to be missing mother info)
                if (origin == 23 || particle->status == 23)
                {
                    int q = copysign(1, particle->pdgId);
                    TLorentzVector lep;

                    if      (abs(particle->pdgId) == 13)
                    {
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        new(genMuonsP4ptr[nGenMuons]) TLorentzVector(lep);

                        nGenMuons++;
                        genMuons.push_back(lep);
                        muonTotalQ += q;
                        genMuonsQ.push_back(q);
                        genMuonStatus.push_back(particle->status);
                    }
                    else if (abs(particle->pdgId) == 11)
                    {
                        lep.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        new(genElectronsP4ptr[nGenElectrons]) TLorentzVector(lep);

                        nGenElectrons++;
                        genElecs.push_back(lep);
                        elecTotalQ += q;
                        genElectronsQ.push_back(q);
                        genElectronStatus.push_back(particle->status);
                    }

                    nGenLeptons++;
                    genLeps.push_back(lep);
                    genLepsQ.push_back(q);
                }
            }
        }



        //--- PHASE SPACE ---//

        TLorentzVector muonSum = GetP4Sum(genMuons), elecSum = GetP4Sum(genElecs), lepSum = GetP4Sum(genLeps);


        // Mass window
        if (lepSum.M() < M_MIN || lepSum.M() > M_MAX)
            return kTRUE;


        // Charge requirement
        if (elecTotalQ != 0 || muonTotalQ != 0)
            return kTRUE;


        // Sort events by decay channel
        unsigned idx = 0;                                   // Index

        if      (nGenMuons == 2 && nGenElectrons == 0)      // mumu = 3
            idx = 3;

        else if (nGenMuons == 0 && nGenElectrons == 2)      // ee   = 4
            idx = 4;

        else if (nGenMuons == 4 && nGenElectrons == 0)      // 4m   = 6
            idx = 6;

        else if (nGenMuons == 2 && nGenElectrons == 2       // 2m2e = 7
                && muonSum.M() > elecSum.M())
            idx = 7;

        else if (nGenMuons == 2 && nGenElectrons == 2       // 2e2m = 8
                && muonSum.M() < elecSum.M())
            idx = 8;

        else if (nGenMuons == 0 && nGenElectrons == 4)      // 4e   = 9
            idx = 9;

        else
            return kTRUE;


        // Dilepton mass requirement for 4l events
        if (nGenMuons == 2 && nGenElectrons == 2)           // Mixed flavor
        {
            if (elecSum.M() < MLL_MIN)
                return kTRUE;

            if (muonSum.M() < MLL_MIN)
                return kTRUE;
        }
        else if (nGenMuons == 4 || nGenElectrons == 4)      // Single flavor
        {
            for (unsigned j = 1; j < 4; j++)
            {
                for (unsigned i = 0; i < j; i++)
                {
                    if (genLepsQ[i] != genLepsQ[j])
                    {
                        TLorentzVector dilep = genLeps[i] + genLeps[j];
                        if (dilep.M() < MLL_MIN)
                            return kTRUE;
                    }
                }
            }
        }


        // Remaining events must be in phase space region
        hPhaseSpaceEvents->Fill(idx, genWeight);
        hFiducialEvents->Fill(1, genWeight);



        //--- FIDUCIAL REGION ---//


        sort(genLeps.begin(), genLeps.end(), P4SortCondition); 


        // Eta
        for (unsigned i = 0; i < nGenLeptons; i++)
        {
            if (fabs(genLeps[i].Eta()) > ETA_MAX)
                isFiducial = kFALSE;
        }


        // Lepton Pt
        if (genLeps[0].Pt() < PT1_MIN)
            isFiducial = kFALSE;

        if (genLeps[1].Pt() < PT2_MIN)
            isFiducial = kFALSE;
        
        if (nGenLeptons == 4)
        {
            if (genLeps[3].Pt() < PT_MIN || genLeps[2].Pt() < PT_MIN)
                isFiducial = kFALSE;
        }


        // Remaining events must be in fiducial region
        if (isFiducial)
            hFiducialEvents->Fill(idx, genWeight);
    }




    /////////////////////
    //    FILL TREE    //
    /////////////////////


    outTree->Fill();
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




///////////////////
//    HELPERS    //
///////////////////


TLorentzVector GetP4Sum(const vector<TLorentzVector> &p4)
{
    TLorentzVector p4sum;

    for (unsigned i = 0; i < p4.size(); i++)
        p4sum += p4[i];

    return p4sum;
}
