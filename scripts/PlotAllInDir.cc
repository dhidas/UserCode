////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Tue Aug 25 14:30:21 CEST 2009
//
// To compile:
//   g++ -Wall `root-config --cflags` `root-config --libs` -lMinuit PlotAllInDir.cc -o PlotAllInDir
//
// To Run:
//   ./PlotAllInDir eps PlotElectrons/ TestInFile.root OutDir
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"



int PlotAllInDir (TString const Suffix, TString const FileName, TString const DirName = "", TString const OutDir = "")
{
  // Open File
  TFile File(FileName, "read");
  if (!File.IsOpen()) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    return 1;
  }

  if (DirName == "") {
    File.ls();
    return 0;
  }

  // Get the directory of interest
  TDirectory* Dir;
  if (DirName.BeginsWith(".")) {
    Dir = (TDirectory*) File.GetDirectory("");
  } else {
    Dir = (TDirectory*) File.Get(DirName);
  }
  if (Dir == 0x0) {
    std::cerr << "ERROR: cannot get directory" << std::endl;
    return 1;
  }

  Dir->ls();
  if (OutDir != "") {
    gSystem->mkdir(OutDir, kTRUE);
  }

  // Set some basic style options
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);

  // Get contents and loop over them
  TIter MyIter(Dir->GetListOfKeys());
  TKey* Key;
  while ( (Key = (TKey*) MyIter()) ) {
    TObject *Object = (TObject*) Key->ReadObj();
    std::cout << "ClassName: " << Object->ClassName() << std::endl;

    TCanvas Can;
    TString Name;
    TString const ClassName(Object->ClassName());
    if (ClassName == "TH1D" || ClassName == "TH1F") {
      TH1D* Hist = (TH1D*) Object;
      Name = Hist->GetName();
      Can.cd();
      Hist->SetFillColor(9);
      //Hist->Draw("hist");
      Hist->Draw("ep");
      Can.SetLogy(0);
    } else if (ClassName == "TH2D" || ClassName == "TH2F") {
      TH2D* Hist = (TH2D*) Object;
      Name = Hist->GetName();
      Can.cd();
      Hist->Draw();
    } else {
      continue;
    }
    std::cout << "Name: " << Name << std::endl;
    if (OutDir == "") {
      Can.SaveAs(Name+"."+Suffix);
    } else {
      Can.SaveAs(OutDir+"/"+Name+"."+Suffix);
    }


  }
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3 && argc != 4 && argc != 5) {
    std::cerr << "Usage: " << argv[0] << " [eps|gif..] [DirName] [FileName] ([OutDir])" << std::endl;
    return 1;
  }

  if (argc == 3) {
    PlotAllInDir(argv[1], argv[2]);
  } else if (argc == 4) {
    PlotAllInDir(argv[1], argv[3], argv[2]);
  } else if (argc == 5) {
    PlotAllInDir(argv[1], argv[3], argv[2], argv[4]);
  }

  return 0;
}
