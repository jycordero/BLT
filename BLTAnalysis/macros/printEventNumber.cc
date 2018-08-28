#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

using namespace baconhep;
using namespace std;



void printEventNumber(const UInt_t minIdx, const UInt_t nFileMax, const TString dirName, const TString rootFilePref)
{
    // Output text file
    TString textFileName = rootFilePref + "_EventNumbers.txt";
    fstream buff(textFileName, fstream::out);


    // ROOT file info
    TString eosPath = "root://cmsxrootd.fnal.gov/";


    // Loop over ROOT files
    UInt_t nFileMin = 0, nEvents = 0;
    for (unsigned i = nFileMin; i <= nFileMax; i++)
    {
        TString rootFileName;
        rootFileName.Form(rootFilePref + "_%i.root", i);
        cout << "Processing " << rootFileName << "..." << flush;

        TFile *rootFile = TFile::Open(eosPath + "/" + rootFileName);
        TTreeReader reader("Events", rootFile);
        TTreeReaderValue<TEventInfo> fInfo(reader, "Info");


        // Loop over events
        UInt_t eventCount = 0;
        while (reader.Next())
        {
            eventCount++;
            UInt_t eventNumber = fInfo->evtNum;
            buff << eventNumber << endl;
        }
        cout << eventCount << " events" << endl;


        // Cleanup
        nEvents += eventCount;
        rootFile->Close();
        delete rootFile;
    }
    buff.close();


    cout << endl;
    cout << "Ran over " << nEvents << " total events" << endl;
    cout << "Wrote output to " << textFileName << endl;
}
