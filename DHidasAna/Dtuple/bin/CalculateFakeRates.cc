////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hida@cern.ch>
//
// Created on: Tue Apr 27 14:39:24 CEST 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"


int CalculateFakeRates ()
{
  // Calculate the fake rate parameterized however you want...

  // Open the relavent files
  TFile QCDFileA("SimpleAna_QCD.root", "read");
  TFile EWKWFile("SimpleAna_Wenu.root", "read");
  TFile EWKZFile("SimpleAna_Zee.root", "read");

  // The output file for fake rates
  TFile OutFile("FakeRates.root", "recreate");

  if ( !(QCDFileA.IsOpen() && EWKWFile.IsOpen() && EWKZFile.IsOpen()) ) {
    std::cerr << "ERROR: cannot open one of the files" << std::endl;
    return 1;
  }

  // Hist Names to use
  TString const NumerHistName = "PlotFakes/EleFakeNumeratorEtaPhi";
  TString const DenomHistName = "PlotFakes/EleFakeDenomJetEtaPt";

  // Histograms from qcd samples
  TH2D* hNumer = (TH2D*) QCDFileA.Get(NumerHistName);
  TH2D* hDenom = (TH2D*) QCDFileA.Get(DenomHistName);

  // Histograms from EWK contributions
  TH2D* hNumerW = (TH2D*) EWKWFile.Get(NumerHistName);
  TH2D* hDenomW = (TH2D*) EWKWFile.Get(DenomHistName);

  // Subtract the EWK contribtions from numerator and denom
  //hNumer->Add(hNumerW, -1);
  //hDenom->Add(hDenomW, -1);

  // The rate hist
  TH2D* hRate = (TH2D*) hNumer->Clone(0);
  //hRate->Clear();
  hRate->Divide(hDenom);

  OutFile.cd();
  hRate->SetTitle("FakeRate");
  hRate->SetName("FakeRate");
  hRate->SetDirectory(&OutFile);

  OutFile.Write();
  OutFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  CalculateFakeRates();

  return 0;
}
