////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Oct 27 15:21:28 EDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"


void MakePlotsFor (TString const HistName, TFile& fWpLO, float const xsWpLO,
                                           TFile& fWmLO, float const xsWmLO,
                                           TFile& fWpNLO, float const xsWpNLO,
                                           TFile& fWmNLO, float const xsWmNLO)
{
  TH1F* hWpLO = (TH1F*) fWpLO.Get(HistName);
  TH1F* hWmLO = (TH1F*) fWmLO.Get(HistName);
  TH1F* hWpNLO = (TH1F*) fWpNLO.Get(HistName);
  TH1F* hWmNLO = (TH1F*) fWmNLO.Get(HistName);

  hWpLO->Scale(xsWpLO/hWpLO->Integral());
  hWmLO->Scale(xsWmLO/hWmLO->Integral());
  hWpNLO->Scale(xsWpNLO/hWpNLO->Integral());
  hWmNLO->Scale(xsWmNLO/hWmNLO->Integral());

  TH1F* hLO = (TH1F*) hWpLO->Clone();
  hLO->Add(hWmLO);
  TH1F* hNLO = (TH1F*) hWpNLO->Clone();
  hNLO->Add(hWmNLO);

  TH1F* kFactor = (TH1F*) hNLO->Clone();
  kFactor->Divide(hLO);
  kFactor->SetLineColor(4);
  kFactor->SetFillColor(0);
  kFactor->SetMaximum(2);

  hLO->SetLineColor(1);
  hNLO->SetLineColor(2);
  hLO->SetFillColor(0);
  hNLO->SetFillColor(0);

  TH1F* hLO_Norm  = (TH1F*) hLO->Clone();
  TH1F* hNLO_Norm = (TH1F*) hNLO->Clone();
  hLO_Norm->SetNormFactor(1);
  hNLO_Norm->SetNormFactor(1);

  TCanvas Can(HistName, HistName, 300, 600);
  Can.Divide(1,3);

  Can.cd(1);
  hLO->SetYTitle("d#sigma/d[]");
  hLO->Draw("hist");
  hNLO->Draw("histsame");

  Can.cd(2);
  hLO_Norm->SetYTitle("");
  hLO_Norm->Draw("hist");
  hNLO_Norm->Draw("histsame");

  Can.cd(3);
  kFactor->SetYTitle("NLO/LO");
  kFactor->Draw("hist");

  TString MyName = HistName;
  MyName = MyName(MyName.Last('/')+1, MyName.Length()-MyName.Last('/')-1);
  Can.SaveAs(MyName+".eps");

  std::cout << "k-factor for " << HistName << " is " << hNLO->Integral() / hLO->Integral() << std::endl;

  delete hLO;
  delete hNLO;
  delete kFactor;

  return;
}


int MakeWgPlots (TString const InDir)
{
  TFile fWpLO( InDir+"SM_LO_Wp.root", "read");
  TFile fWmLO( InDir+"SM_LO_Wm.root", "read");
  TFile fWpNLO(InDir+"SM_NLO_Wp.root", "read");
  TFile fWmNLO(InDir+"SM_NLO_Wm.root", "read");

  float const xsWpLO  =  8.78;
  float const xsWmLO  = 11.8;
  float const xsWpNLO = 13.9;
  float const xsWmNLO =  9.46;

  MakePlotsFor("LHEPlots/PhotonEt", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/PhotonEta", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/LeptonPt", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/LeptonEta", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);




  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InDir]" << std::endl;
    return 1;
  }

  TString const InDir = argv[1];

  MakeWgPlots(InDir);

  return 0;
}
