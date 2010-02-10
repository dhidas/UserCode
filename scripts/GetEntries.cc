////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Feb 10 18:57:59 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"


int GetEntries (TString const HistName, TString const FileName)
{
  TFile MyFile(FileName, "read");
  if (!MyFile.IsOpen()) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    exit(1);
  }

  TH1D* Hist = (TH1D*) MyFile.Get(HistName);
  printf("%12i  Int:%12.4f  u:%12.4f  o:%12.4f\n", (int)
         Hist->GetEntries(),
         Hist->Integral(),
         Hist->GetBinContent(0),
         Hist->GetBinContent( Hist->GetNbinsX() + 1));

  delete Hist;
  MyFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [HistName] [FileName]" << std::endl;
    return 1;
  }

  GetEntries(argv[1], argv[2]);

  return 0;
}
