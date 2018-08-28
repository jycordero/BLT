#include <fstream>

#include "TString.h"

using namespace std;



void findDoubledEvents(const TString sampleName, const UInt_t nEvents)
{
    // Input text file
    TString inFileName = sampleName + ".txt", outFileName = sampleName + "_dupes.txt";
    ifstream file(inFileName);
    fstream buff(outFileName, fstream::out);


    // Double loop to find duplicates
    UInt_t *eventNumber = new UInt_t[nEvents];
    UInt_t *runNumber   = new UInt_t[nEvents];
    UInt_t *lumiSection = new UInt_t[nEvents];
    UInt_t nDupes = 0;

    cout << "Looping over " << nEvents << " events..." << endl;
    for (unsigned i = 0; i < nEvents; i++)
    {
        if (i % 10000 == 0)
            cout << "Processed " << i << " events" << endl;

        file >> eventNumber[i] >> runNumber[i] >> lumiSection[i];
//      cout << eventNumber[i] << "\t" << runNumber[i] << "\t" << lumiSection[i] << endl;
        for (unsigned j = 0; j < i; j++)
        {
            if (eventNumber[i] == eventNumber[j])
            {
                nDupes++;
                buff << eventNumber[i] << "\t" << runNumber[i] << "\t" << lumiSection[i] << endl;
            }
        }
    }
    file.close();
    buff.close();
    delete[] eventNumber;
    delete[] runNumber;
    delete[] lumiSection;

    cout << "Found " << nDupes << " duplicate events" << endl;
}
