#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TRandom3.h"

using namespace std;

void makeSmearMuonHZZ()
{
    // Get histograms
 
    TString inPath = "/uscms/home/jrainbol/nobackup/CMSSW_7_4_14/src/BLT/BLTAnalysis/data/muon_id/";
    TFile *inFile = TFile::Open(inPath + "hzz_muon_id_sf.root");
    TH2 *h_err;
    inFile->GetObject("ERROR", h_err);

//  h_err->Draw("COL");



    // Retrieve pt, eta binning

    Int_t x_bins, y_bins;
    x_bins = h_err->GetNbinsX();
    y_bins = h_err->GetNbinsY();

    Double_t x_edge[20], y_edge[20];
    h_err->GetXaxis()->GetLowEdge(x_edge);
    h_err->GetYaxis()->GetLowEdge(y_edge);

    // Have to do the overflow bin manually because ROOT sucks
    x_edge[x_bins] = h_err->GetXaxis()->GetBinLowEdge(x_bins+1);
    y_edge[y_bins] = h_err->GetYaxis()->GetBinLowEdge(y_bins+1);

/*
    // Print binning

    cout << "eta: ";
    for (unsigned i = 0; i < x_bins+1; i++)
        cout << x_edge[i] << ", ";
    cout << endl;

    cout << "pt: ";
    for (unsigned j = 0; j < y_bins+1; j++)
        cout << y_edge[j] << ", ";
    cout << endl;
*/


    // Seed RNG (don't change!!)

    TRandom3 rng(13);



    // Create and fill histogram

    TH2D *h_smr = new TH2D("SMEAR", "SMEAR", x_bins, x_edge, y_bins, y_edge);
    h_smr->SetDirectory(0);

    for (unsigned i = 1; i < x_bins+1; i++)
    {
        for (unsigned j = 1; j < y_bins; j++)
        {
            Double_t err = h_err->GetBinContent(i, j);
            Double_t smr = rng.Gaus(0, err);

            h_smr->SetBinContent(i, j, smr);
        }
    }

//  h_smr->Draw("COLZ");



    // Write to file
    
    inFile->Close();

    TFile *outFile = new TFile(inPath + "hzz_muon_id_smear.root", "RECREATE");
    h_smr->Write();
    outFile->Close();
}
