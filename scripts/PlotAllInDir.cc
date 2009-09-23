////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Tue Aug 25 14:30:21 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"


int PlotAllInDir (TString const DirName, TString const FileName)
{
  // Open File
  TFile File(FileName, "read");
  if (!File.IsOpen()) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    return 1;
  }

  // Get the directory of interest
  TDirectory* Dir = (TDirectory*) File.Get(DirName);
  if (Dir == 0x0) {
    std::cerr << "ERROR: cannot get directory" << std::endl;
    return 1;
  }

  // Set some basic style options
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(11111111);

  // Get contents and loop over them
  TIter MyIter(Dir->GetListOfKeys());
  TKey* Key;
  while ( (Key = (TKey*) MyIter()) ) {
    TObject *Object = (TObject*) Key;
    if (Object->ClassName() != "TH1D") {
      continue;
    }
    TH1D* Hist = (TH1D*) Object;
    TString const Name = Hist->GetName();

    TCanvas Can;
    Can.cd();
    Hist->Draw();
    Can.SaveAs(Name+".gif");


  }
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [DirName] [FileName]" << std::endl;
    return 1;
  }

  PlotAllInDir(argv[1], argv[2]);

  return 0;
}
