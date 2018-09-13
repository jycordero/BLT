#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

void sumEvents(const unsigned nFiles, const unsigned nStart, const TString filePath, const TString suffix)
{
    Long64_t nEvents = 0, nTotal = 0;
    TString filePath, num;

    for (unsigned i = nStart; i < nStart + nFiles; i++)
    {
        num.Form("%i", i);
        nEvents = 0;
        TString fullName = filePath + "_" + num + ".root";
        
        TFile *file = TFile::Open(fullName);
        TH1D *hist;

        file->GetObject("TotalEvents_" + suffix, hist);
        nEvents = hist->GetBinContent(1);

        cout << fullName << ": \t" << nEvents << " events" << endl;

        nTotal += nEvents;

        file->Close();
    }

    cout << "All " << nFiles << " files:\t" << nTotal << " events" << endl;
}
