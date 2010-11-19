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


void MakePlotsFor (TString const HistName, TFile& fLO,  float const xsLO,
                                           TFile& fNLO, float const xsNLO)
{
  bool const Together = false;

  TH1F* hLO  = (TH1F*) fLO.Get(HistName);
  TH1F* hNLO = (TH1F*) fNLO.Get(HistName);

  hLO->Scale(xsLO/hLO->Integral());
  hNLO->Scale(xsNLO/hNLO->Integral());

  hLO->SetLineColor(1);
  hNLO->SetLineColor(2);
  hLO->SetFillColor(0);
  hNLO->SetFillColor(0);

  TH1F* kFactor = (TH1F*) hNLO->Clone();
  kFactor->Divide(hLO);
  kFactor->SetLineColor(4);
  kFactor->SetFillColor(0);
  kFactor->SetMaximum(6);
  kFactor->SetMinimum(0);


  TH1F* hLO_Norm  = (TH1F*) hLO->Clone();
  TH1F* hNLO_Norm = (TH1F*) hNLO->Clone();
  hLO_Norm->Scale(1/hLO_Norm->Integral());
  hNLO_Norm->Scale(1/hNLO_Norm->Integral());


  TCanvas Can;
  TCanvas Can1;
  Can.Divide(2,2);

  TLegend Leg2(0.5, 0.7, 0.89, 0.89);
  Leg2.SetNColumns(2);
  Leg2.SetFillColor(0);
  Leg2.SetBorderSize(0);
  Leg2.AddEntry(hLO, "LO", "l");
  Leg2.AddEntry(hNLO, "NLO", "l");


  TLegend LegK(0.5, 0.7, 0.89, 0.89);
  LegK.SetNColumns(2);
  LegK.SetFillColor(0);
  LegK.SetBorderSize(0);
  LegK.AddEntry(kFactor, "NLO/LO", "l");

  TString MyName = HistName;
  MyName = MyName(MyName.Last('/')+1, MyName.Length()-MyName.Last('/')-1);

  printf("k-factor for %s is %8.4f\n", HistName.Data(),
      hNLO->Integral() / hLO->Integral());

  if (Together) Can.cd(1)->SetLogy(1); else {Can1.Clear(); Can1.SetLogy(1);}
  hLO->SetYTitle("d#sigma/dP_{T}^{#gamma}");
  if (hLO->GetMaximum() < hNLO->GetMaximum()) {
    hLO->SetMaximum( hNLO->GetMaximum()*1.1 );
  }
  hLO->Draw("hist");
  hNLO->Draw("samehist");
  //Leg2.Draw("same");
  if (!Together)  Can1.SaveAs(MyName+"_XS.eps");

  if (Together) Can.cd(2)->SetLogy(0); else {Can1.Clear(); Can1.SetLogy(0);}
  if (hLO_Norm->GetMaximum() < hNLO_Norm->GetMaximum()) {
    hLO_Norm->SetMaximum( hNLO_Norm->GetMaximum()*1.05 );
  }
  hLO_Norm->SetTitle("Baur Z#gamma");
  hLO_Norm->SetYTitle("Normalized to Unity");
  hLO_Norm->Draw("hist");
  hNLO_Norm->Draw("histsame");
  Leg2.Draw("same");
  if (!Together) Can1.SaveAs(MyName+"_Norm.eps");

  if (Together) Can.cd(3); else Can1.Clear();
  kFactor->SetTitle("k-factor");
  kFactor->SetYTitle("NLO/LO");
  kFactor->Draw("hist");
  LegK.Draw("same");
  if (!Together) Can1.SaveAs(MyName+"_kFactor.eps");

  if (Together) Can.cd(4)->SetLogy(1); else {Can1.Clear(); Can1.SetLogy(1);}
  hLO_Norm->Draw("hist");
  hNLO_Norm->Draw("histsame");
  Leg2.Draw("same");
  if (!Together) Can1.SaveAs(MyName+"_NormLog.eps");



  if (Together) Can.SaveAs(MyName+".eps");

  delete hLO;
  delete hNLO;
  delete kFactor;

  return;
}


int MakeZgPlots (TString const InDir)
{
  TFile fLO(  InDir+"/"+"SM_LO_Zg.root", "read");
  TFile fNLO( InDir+"/"+"SM_NLO_Zg.root", "read");

  // with Vg cuts. all from NLO prog (born wgts used for LO)
  float const xsLO  = 3.61;
  float const xsNLO = 5.12;


  // Set some style bullshit
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.09);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);


  MakePlotsFor("LHEPlots/PhotonEt", fLO, xsLO, fNLO, xsNLO);
  MakePlotsFor("LHEPlots/PhotonEta", fLO, xsLO, fNLO, xsNLO);
  MakePlotsFor("LHEPlots/LeptonEta", fLO, xsLO, fNLO, xsNLO);
  MakePlotsFor("LHEPlots/LeptonPt", fLO, xsLO, fNLO, xsNLO);
  //MakePlotsFor("LHEPlots/LeptonPhotonDeltaR", fWpLO, xsWpLO, fWmLO, xsWmLO, fWpNLO, xsWpNLO, fWmNLO, xsWmNLO);




  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InDir]" << std::endl;
    return 1;
  }

  TString const InDir = argv[1];

  MakeZgPlots(InDir);

  return 0;
}
