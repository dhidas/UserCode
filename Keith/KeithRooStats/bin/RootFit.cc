////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Aug 11 11:58:17 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"


int RootFit (TString const InFileName)
{
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file" << std::endl;
    throw;
  }

  TH1F* Hist = (TH1F*) InFile.Get("Extinction6_Ratio");

  TF1 Func("Func", "0.6 / (1.0 + TMath::Exp((x - [0]) / [1])) + 0.4", 507., 1684.);
  Func.SetParameter(0, 500);
  Func.SetParameter(1, 50);

  Hist->Fit("Func");
  printf("Fit Parameters:\n");
  printf("  CutOffPt:    %12.3E\n", Func.GetParameter(0));
  printf("  CutOffSpeed: %12.3E\n", Func.GetParameter(1));

  TCanvas Can;
  Can.cd();
  Hist->Draw();
  Can.SaveAs("2TeVFit.gif");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  RootFit(argv[1]);

  return 0;
}
