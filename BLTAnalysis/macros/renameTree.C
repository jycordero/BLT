#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void renameTree(const TString filename, const TString oldname, const TString newname)
{
    TFile file(filename, "UPDATE");
    TTree *tree;
    file.GetObject(oldname, tree);
    tree->Write(newname);
    cout << "Copied " + oldname + " to " + newname + "." << endl;
    file.Close();

//  delete tree;
}
