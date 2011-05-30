////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue May  3 16:08:52 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"


void SetStyle ()
{
  TStyle *DStyle = new TStyle("DStyle", "Dean Style Punk");

  DStyle->SetCanvasBorderMode(0);
  DStyle->SetPadBorderMode(0);
  DStyle->SetPadColor(0);
  DStyle->SetLineColor(0);
  //DStyle->SetFillColor(0);
  DStyle->SetFillStyle(1);
  DStyle->SetCanvasColor(0);
  //DStyle->SetTitleColor(0);
  DStyle->SetTitleFillColor(0);
  DStyle->SetStatColor(0);
  DStyle->SetTitleX(0.1);
  DStyle->SetTitleY(1.0);
  DStyle->SetTitleH(0.09);
  DStyle->SetTitleW(0.8);
  DStyle->SetTitleBorderSize(0);
  DStyle->SetStatY(0.88);
  DStyle->SetStatH(0.25);
  DStyle->SetOptStat(1);
  DStyle->SetCanvasBorderMode(0);

  gROOT->SetStyle("DStyle");
  return;
}



void Check (void* a)
{
  if (!a) {
    std::cerr << "I don't see that!!" << std::endl;
    throw;
  }
  return;
}

TCanvas* PlotForDiag (TFile& InFile, int const Diag, TString const Name = "", TString const Title = "")
{
  TString const PN = Diag >= 0 ? "P" : "N";
  char DD[10];
  sprintf(DD, "d%s%03i", PN.Data(), (int) Diag);
  TString const D = DD;

  TCanvas* c = new TCanvas(Title, Title);
  c->Divide(3,3);

  c->cd(1);
  TH2F* h2 = (TH2F*) InFile.Get("LeptonPlusJets/TriJetSumPt_vs_Mass");
  Check(h2);
  c->GetPad(1)->SetLogz(1);
  h2->SetStats(0);
  //h->SetOptStat(0);
  h2->Draw("colz");
  TLine* line = new TLine(Diag, 0, h2->GetXaxis()->GetXmax(), h2->GetXaxis()->GetXmax()-Diag);
  line->SetLineColor(2);
  line->SetLineWidth(1);
  line->Draw("same");

  c->cd(2);
  h2 = (TH2F*) InFile.Get("LeptonPlusJets/TriJetMas_vs_TripletPt");
  Check(h2);
  h2->SetStats(0);
  c->GetPad(2)->SetLogz(1);
  h2->Draw("colz");

  c->cd(3);
  TH1F* h = (TH1F*) InFile.Get("LeptonPlusJets/TriJetMass");
  Check(h);
  h->Draw("hist");

  c->cd(4);
  h = (TH1F*) InFile.Get("LeptonPlusJets/TriJetMass_"+D);
  Check(h);
  h->SetStats(0);
  h->Draw("hist");

  c->cd(5);
  h = (TH1F*) InFile.Get("LeptonPlusJets/LeptonPt_"+D);
  Check(h);
  h->SetStats(0);
  h->Draw("hist");

  c->cd(6);
  h = (TH1F*) InFile.Get("LeptonPlusJets/METMag_"+D);
  Check(h);
  h->Draw("hist");

  char NAME[200];
  for (int i = 0; i <= 2; ++i) {
    sprintf(NAME, "LeptonPlusJets/JetPt_TripletJet%i_%s", i, D.Data());
    c->cd(7+i);
    h = (TH1F*) InFile.Get(NAME);
    Check(h);
    h->SetStats(0);
    h->Draw("hist");
  }

  return c;
}





void TryPlots (TString const InFileName)
{
  TFile InFile(InFileName, "read");

  int diag;
  TCanvas* c = 0x0;
  std::cout << "Endter Diag cut (-100, 200) (by 10s only): ";
  while (std::cin >> diag) {
    //if (c) delete c; // comment out if you want to see all of them
    c = PlotForDiag(InFile, diag, "Test", "Test");
    //c->SaveAs("PLOT.eps");
    c->Update();

    std::cout << "Endter Diag cut (-100, 200) (by 10s only): ";
  }
  if (c) delete c;

  InFile.Close();

  return;
}


int main (int argc, char* argv[])
{
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " [InFileName] [Diag] [Name] [Title]" << std::endl;
    return 1;
  }

  TString const InFileName =argv[1];
  int const Diag = atoi(argv[2]);
  TString const Name = argv[3];
  TString const Title = argv[4];


  SetStyle();

  TFile InFile(InFileName, "read");

  PlotForDiag(InFile, Diag, Name, Title);

  return 0;
}
