////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  8 14:38:24 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"


int GetHistFromCan (TString const InFileName, TString const OutFileName)
{
  TFile fin(InFileName, "read");
  TCanvas* c = (TCanvas*) fin.Get("c");
  TH1F* hist = (TH1F*) c->FindObject("Mjjj_70_20_160")->Clone("Mjjj_70_20_160_6jet");

  TFile fout(OutFileName, "recreate");
  TH1F* newhist = (TH1F*) hist->Clone();
  newhist->SetDirectory(&fout);
  fout.Write();
  fout.Close();
  fin.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFileName] [OutFileName]" << std::endl;
    return 1;
  }

  GetHistFromCan(argv[1], argv[2]);

  return 0;
}
