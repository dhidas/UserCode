////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Jul 26 13:13:58 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"

int MakeFakeData ()
{
  // Number of entries and xmax
  int   const NEntries = 100000;
  float const XMax     =   2000;
  int   const NBins    =    100;
  // Open output file
  TFile OutFile("FakeData.root", "recreate");

  // Just some signal plus BG model for the ratio
  TF1 SignalModel  ("SignalModel",   "TMath::Exp(-0.006*x) * 1. / (1. + TMath::Exp((x - 1500.)/100.))", 0, XMax);
  TF1 StandardModel("StandardModel", "TMath::Exp(-0.006*x)", 0, XMax);

  TH1F HistB("SM",     "SM",     NBins, 0, XMax);
  TH1F HistB2("SM2",   "SM2",    NBins, 0, XMax);
  TH1F HistS("Signal", "Signal", NBins, 0, XMax);
  HistB.Sumw2();
  HistS.Sumw2();

  HistB.FillRandom("StandardModel", NEntries);
  HistB2.FillRandom("StandardModel", NEntries);

  float Rand;
  for (int i = 0; i < NEntries; ++i) {
    Rand = gRandom->Rndm() * XMax;
    if (gRandom->Rndm() < SignalModel.Eval(Rand)) {
      HistS.Fill(Rand);
    }
  }


  TH1F* HistR = (TH1F*) HistS.Clone("Extinction2");
  HistR->Divide(&HistB);
  HistR->SetDirectory(&OutFile);

  TCanvas C;
  C.cd();
  HistR->Draw();
  C.SaveAs("FakeSignalData.eps");


  TH1F* HistR2 = (TH1F*) HistB.Clone("StandardModel2");
  HistR2->Divide(&HistB2);
  HistR2->SetDirectory(&OutFile);
  C.cd();
  HistR2->Draw();
  C.SaveAs("FakeSMData.eps");


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

  MakeFakeData();

  return 0;
}
