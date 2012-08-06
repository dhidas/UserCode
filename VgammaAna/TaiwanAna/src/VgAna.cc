#include "VgammaAna/TaiwanAna/interface/VgAna.h"

#include "VgammaAna/TaiwanAna/interface/TAnaHist.h"

#include <iostream>

VgAna::VgAna ()
{
}


VgAna::~VgAna ()
{
  OutFile_->Write();
  OutFile_->Close();
}


void VgAna::SetOutFile (TString const& OutFileName)
{
  OutFile_ = new TFile(OutFileName, "recreate");
  OutFile_->cd();
  return;
}


void VgAna::Loop ()
{

  for (unsigned int ientry = 0; tree_->GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 1000 == 0) {
      printf("Processing entry: %15i\n", ientry);
    }

    PlotElectrons();
    PlotMETs();
  }

  return;
}



void VgAna::PlotElectrons()
{
  static TAnaHist Hist(OutFile_, "PlotElectrons");

  for (int i = 0; i != nEle_; ++i) {
    Hist.FillTH1D("elePt", 20, 0, 200, elePt_[i]);
  }

  return;
}



void VgAna::PlotMETs ()
{
  static TAnaHist Hist(OutFile_, "PlotMETs");

  Hist.FillTH1D("MET", 20, 0, 200, MET_);

  return;
}
