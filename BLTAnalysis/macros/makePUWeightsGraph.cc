#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"

using namespace std;

void makePUWeightsGraph()
{
    TString outFileName = "pileup_2017_69200_100bins.root";



    // Get histograms from files
 
    TString path = "../data/pileup/";
    TString dataFileName = "DataPileupHistogram2017_69200_100bins.root", dataHistName = "pileup";
    TString mcFileName = "histProbFunction.root", mcHistName = "RunIIFall17";

    TH1 *dataHist, *mcHist;

    TFile *dataFile = new TFile(path + dataFileName);
    dataFile->GetObject(dataHistName, dataHist);
    dataHist->SetDirectory(0);
    dataFile->Close();

    TFile *mcFile = new TFile(path + mcFileName);
    mcFile->GetObject(mcHistName, mcHist);
    mcHist->SetDirectory(0);
    mcFile->Close();



    // Create ratio histogram and graph

    TH1D *normHist = (TH1D*) dataHist->Clone();
    normHist->Divide(dataHist, mcHist, 1./dataHist->Integral(), 1./mcHist->Integral());
    TGraph *normGraph = new TGraph(normHist);
    normGraph->SetName("pileup_sf");



    // A couple of checks 

    normHist->Draw("HIST");
    cout << normGraph->GetMean(2) << endl;



    // Write everything to file

    TFile *outFile = new TFile(path + outFileName, "RECREATE");
    normGraph->Write();
    dataHist->Write();
    mcHist->Write();
    outFile->Close();
}
