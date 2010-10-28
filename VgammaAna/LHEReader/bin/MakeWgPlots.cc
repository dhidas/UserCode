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
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"


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

  hWpLO->SetLineColor(1);
  hWmLO->SetLineColor(1);
  hWpNLO->SetLineColor(2);
  hWmNLO->SetLineColor(2);
  hWpLO->SetFillColor(0);
  hWmLO->SetFillColor(0);
  hWpNLO->SetFillColor(0);
  hWmNLO->SetFillColor(0);

  TH1F* hLO = (TH1F*) hWpLO->Clone();
  hLO->Add(hWmLO);
  TH1F* hNLO = (TH1F*) hWpNLO->Clone();
  hNLO->Add(hWmNLO);

  TH1F* kFactor = (TH1F*) hNLO->Clone();
  kFactor->Divide(hLO);
  kFactor->SetLineColor(4);
  kFactor->SetFillColor(0);
  kFactor->SetMaximum(4);

  TH1F* kFactorP = (TH1F*) hWpNLO->Clone();
  kFactorP->Divide(hWpLO);
  kFactorP->SetLineStyle(2);
  kFactorP->SetLineColor(6);
  TH1F* kFactorM = (TH1F*) hWmNLO->Clone();
  kFactorM->Divide(hWmLO);
  kFactorM->SetLineStyle(3);
  kFactorM->SetLineColor(8);

  TH1F* hLO_Norm  = (TH1F*) hLO->Clone();
  TH1F* hNLO_Norm = (TH1F*) hNLO->Clone();
  hLO_Norm->SetNormFactor(1);
  hNLO_Norm->SetNormFactor(1);


  TCanvas Can(HistName, HistName, 300, 600);
  Can.Divide(1,3);

  TLegend Leg4(0.5, 0.7, 0.9, 0.9);
  Leg4.SetNColumns(2);
  Leg4.SetFillColor(0);
  Leg4.SetBorderSize(0);
  Leg4.AddEntry(hLO, "LO", "l");
  Leg4.AddEntry(hNLO, "NLO", "l");
  Leg4.AddEntry(hWpLO, "LO W^{+}", "l");
  Leg4.AddEntry(hWpNLO, "NLO W^{+}", "l");

  TLegend Leg2(0.5, 0.8, 0.9, 0.9);
  Leg2.SetNColumns(2);
  Leg2.SetFillColor(0);
  Leg2.SetBorderSize(0);
  Leg2.AddEntry(hLO, "LO", "l");
  Leg2.AddEntry(hNLO, "NLO", "l");

  TLegend LegK(0.5, 0.7, 0.9, 0.9);
  LegK.SetNColumns(2);
  LegK.SetFillColor(0);
  LegK.SetBorderSize(0);
  LegK.AddEntry(kFactor, "NLO/LO", "l");
  LegK.AddEntry(kFactorP, "NLO/LO W^{+}", "l");
  LegK.AddEntry(kFactorM, "NLO/LO W^{-}", "l");

  Can.cd(1);
  hLO->SetYTitle("d#sigma");
  if (hLO->GetMaximum() < hNLO->GetMaximum()) {
    hLO->SetMaximum( hNLO->GetMaximum()*1.1 );
  }
  hLO->Draw("hist");
  hWpLO->SetLineStyle(2);
  hWpLO->Draw("samehist");
  hWpNLO->SetLineStyle(2);
  hWpNLO->Draw("samehist");
  hNLO->Draw("histsame");
  Leg4.Draw("same");

  Can.cd(2);
  if (hLO_Norm->GetMaximum() < hNLO_Norm->GetMaximum()) {
    hLO_Norm->SetMaximum( hNLO_Norm->GetMaximum()*1.05 );
  }
  hLO_Norm->SetTitle("Normalized to unity");
  hLO_Norm->Draw("hist");
  hNLO_Norm->Draw("histsame");
  Leg2.Draw("same");

  Can.cd(3);
  kFactor->SetTitle("k-factor");
  kFactor->SetYTitle("NLO/LO");
  kFactor->Draw("hist");
  kFactorP->Draw("samehist");
  kFactorM->Draw("samehist");
  LegK.Draw("same");

  TString MyName = HistName;
  MyName = MyName(MyName.Last('/')+1, MyName.Length()-MyName.Last('/')-1);
  Can.SaveAs(MyName+".eps");

  printf("k-factor for All/W+/W- for %s is %8.4f / %8.4f / %8.4f\n", HistName.Data(),
      hNLO->Integral() / hLO->Integral(),
      hWpNLO->Integral() / hWpLO->Integral(), 
      hWmNLO->Integral() / hWmLO->Integral());

  delete hLO;
  delete hNLO;
  delete kFactor;
  delete kFactorP;
  delete kFactorM;

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

  // Set some style bullshit
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.09);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);


  MakePlotsFor("LHEPlots/PhotonEt", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/PhotonEta", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/LeptonPt", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/LeptonEta", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);
  MakePlotsFor("LHEPlots/LeptonPhotonDeltaR", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);




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
