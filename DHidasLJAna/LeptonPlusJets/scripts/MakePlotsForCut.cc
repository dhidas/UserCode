////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed May 25 14:58:44 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
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
  DStyle->SetOptStat(11111111);
  DStyle->SetCanvasBorderMode(0);

  gROOT->SetStyle("DStyle");
  return;
}


int MakePlotsForCut (TString const InFileName, TString const BaseName)
{
  SetStyle();
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open file: " << InFileName << std::endl;
    throw;
  }

  for (int cut = 0; cut <= 200; cut += 10) {
    TString const CutStr = TString::Format("%03i", cut);
    TString PlotName = BaseName + "_d" + (cut < 0 ? "N" : "P") + CutStr;


    TCanvas c;
    c.cd();
    InFile.Get(PlotName)->Draw();
    TString SaveName = PlotName( PlotName.Last('/')+1, PlotName.Length() - PlotName.Last('/') - 1) + ".gif";
    c.SaveAs(SaveName);
    
  }
  TString PlotName = BaseName;
  TCanvas c;
  c.cd();
  InFile.Get(PlotName)->Draw();
  TString SaveName = PlotName( PlotName.Last('/')+1, PlotName.Length() - PlotName.Last('/') - 1) + ".gif";
  c.SaveAs(SaveName);


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  TString const BaseName = argv[2];

  MakePlotsForCut(InFileName, BaseName);

  return 0;
}
