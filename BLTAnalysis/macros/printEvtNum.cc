#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

using namespace std;



void printEvtNum(const TString rootFilePref, const UInt_t startIdx, const UInt_t nFiles, const TString inDir, const TString suffix)
{
    // Output text file
    TString textFileName = rootFilePref + "_EvtNum.txt";
    fstream buff(textFileName, fstream::out);


    // Loop over ROOT files
    TString eosPath = "root://cmsxrootd.fnal.gov/";
    UInt_t nEvents = 0;
    for (unsigned i = startIdx; i < nFiles + startIdx; i++)
    {
        TString rootFileName;
        rootFileName.Form(rootFilePref + "_%i.root", i);

        TFile *rootFile = TFile::Open(eosPath + inDir + "/" + rootFileName);
        TTreeReader reader("tree_" + suffix, rootFile);
        TTreeReaderValue<ULong64_t> evtNum(reader, "evtNumber.eventNumber");
        TTreeReaderValue<Int_t> lumiSec(reader, "lumiSection");
        TTreeReaderValue<Int_t> runNum(reader, "runNumber");


        // Loop over events
        // Print event index (for current tree), file index (for current file), and event number (global)
        cout << "Reading " << rootFileName << "..." << flush;

        UInt_t eventCount = 0;
        while (reader.Next())
        {
            eventCount++;
//          buff << reader.GetCurrentEntry() << '\t' << i << '\t';
            buff << *evtNum << '\t' << *lumiSec << '\t' << *runNum << endl;
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
    cout << "Wrote output to " << textFileName << endl << endl << endl;
}
