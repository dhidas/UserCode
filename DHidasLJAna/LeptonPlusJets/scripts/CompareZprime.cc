////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu May 19 06:03:58 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

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

int CompareZprime (TString const PlotName)
{
  // Set the style you want
  SetStyle();

  // Directory where files are
  TString const Dir = "/cms/data21/dhidas/TestPlots/";

  // Filenames, masses, root files
  std::vector<TString> FileName;
  FileName.push_back("ZprimeM500PAT.root");
  FileName.push_back("ZprimeM750PAT.root");
  FileName.push_back("ZprimeM1000PAT.root");
  FileName.push_back("ZprimeM1250PAT.root");
  FileName.push_back("ZprimeM1500PAT.root");
  FileName.push_back("ZprimeM2000PAT.root");
  FileName.push_back("ZprimeM3000PAT.root");
  FileName.push_back("ZprimeM4000PAT.root");
  FileName.push_back("TTJets.root");
  std::vector<TString> Masses;
  Masses.push_back("500");
  Masses.push_back("750");
  Masses.push_back("1000");
  Masses.push_back("1250");
  Masses.push_back("1500");
  Masses.push_back("2000");
  Masses.push_back("3000");
  Masses.push_back("4000");
  Masses.push_back("TTJets");
  std::vector<TFile*> Files;
  //Files.resize(8);
  for (size_t i = 0; i != FileName.size(); ++i) {
    Files.push_back(new TFile(Dir + FileName[i], "read"));
  }


  // Get a canvas and split it
  TCanvas Can(PlotName, PlotName);
  Can.Divide(3,3);

  // Loop over all fies
  for (size_t i = 0; i != Files.size(); ++i) {

    // Grab pad for log or not, done later
    TVirtualPad* Pad = Can.cd(i+1);

    // Get object and check if it's there...
    TObject* h = (TObject*) Files[i]->Get(PlotName);
    if (!h) {
      std::cerr << "ERROR: cannot get plot for file: " << FileName[i] << std::endl;
      continue;
    }

    // Set the title to the mass
    ((TH1F*) h)->SetTitle(Masses[i]);

    // For the 2D let's make it look nicer..
    if (TString(h->ClassName()).BeginsWith("TH2")) {
      // do some stuffs
      //gStyle->SetPalette(1);
      Pad->SetLogz(1);
      ((TH2F*) h)->Rebin2D(4);
      h->Draw("colz");

    } else {
      h->Draw();
    }
  }

  TString SaveName = PlotName( PlotName.Last('/')+1, PlotName.Length() - PlotName.Last('/') - 1) + ".gif";
  // Save the image
  Can.SaveAs(SaveName);


  // Clean up a little..
  for (size_t i = 0; i != Files.size(); ++i) {
    delete Files[i];
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [PlotName]" << std::endl;
    return 1;
  }

  CompareZprime(argv[1]);

  return 0;
}
