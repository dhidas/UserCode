////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed May 25 16:03:50 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <stdlib.h>

#include "TFile.h"
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




int Plot2DWithCut (TString const InFileName, int const Diag, TString const Title = "")
{
  SetStyle();
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    throw;
  }

  TH2F* h2 = (TH2F*) InFile.Get("LeptonPlusJets/TriJetSumPt_vs_Mass_NJetGE04");
  if (h2 == 0x0) {
    std::cerr << "Cannot find hist" << std::endl;
    throw;
  }

  TCanvas c;
  gStyle->SetPalette(0);
  c.cd();
  h2->SetXTitle("Tri-jet SumPt");
  h2->SetYTitle("Tri-jet Mass");
  h2->SetTitle(Title);
  h2->Rebin2D(8);
  c.SetLogz(1);
  h2->SetStats(0);
  //h->SetOptStat(0);
  h2->Draw("colz");
  TLine* line = new TLine(Diag, 0, h2->GetXaxis()->GetXmax(), h2->GetXaxis()->GetXmax()-Diag);
  line->SetLineColor(2);
  line->SetLineWidth(1);
  line->Draw("same");
  c.SaveAs("TriJetMassVsSumPt.gif");




  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " [InFileName] [Diag]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Diag = atoi(argv[2]);
  TString const Title = argc == 4 ? argv[3] : "";

  Plot2DWithCut(InFileName, Diag, Title);

  return 0;
}
