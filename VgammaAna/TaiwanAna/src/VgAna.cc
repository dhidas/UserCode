#include "VgammaAna/TaiwanAna/interface/VgAna.h"

#include "VgammaAna/TaiwanAna/interface/TAnaHist.h"

#include "TLorentzVector.h"

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
    PlotJets();
    PlotElectronPhoton();
  }

  return;
}



void VgAna::PlotElectrons()
{
  // Just make a few electron plots

  // For histogram object
  static TAnaHist Hist(OutFile_, "PlotElectrons");

  // Loop over all electrons and plot whatever you want
  for (int i = 0; i != nEle_; ++i) {
    Hist.FillTH1D("elePt", 20, 0, 200, elePt_[i]);
  }

  return;
}



void VgAna::PlotMETs ()
{
  // Some plots of whatever MET variables you want

  // Histogram object
  static TAnaHist Hist(OutFile_, "PlotMETs");

  // Plot whatever blah you want here
  Hist.FillTH1D("MET", 20, 0, 200, MET_);

  return;
}




void VgAna::PlotJets ()
{
  // Just some basic plots of jets to see how we're doinghere..

  // Histogram object
  static TAnaHist Hist(OutFile_, "PlotJets");

  // Fill number of jets for each event
  Hist.FillTH1D("nJet", 15, 0, 15, nJet_);

  // Buffer for histogram names
  static char BUFF[150];

  // Loop over all jets and plot whatever you like
  for (int ijet = 0; ijet < nJet_; ++ijet) {
    Hist.FillTH1D("jetEta", 20, -4, 4, jetEta_[ijet]);

    // This plots for individual jets, assuming they are ordered in some sensible way
    sprintf(BUFF, "jetEt_%i", ijet);
    Hist.FillTH1D(BUFF, 20, -4, 4, jetEt_[ijet]);

  }


  return;
}




void VgAna::PlotElectronPhoton ()
{
  // Plot events with an electron and a photon

  // Histogram object
  static TAnaHist Hist(OutFile_, "PlotElectronPhoton");

  // Loop over electrons and pick out the ones we want to look at
  std::vector<int> Ele;
  for (int ielectron = 0; ielectron < nEle_; ++ielectron) {
    if (elePt_[ielectron] > 20.0 &&
        eleIsoEcalDR03_[ielectron] + eleIsoHcalDR03_[ielectron] + eleIsoTrkDR03_[ielectron] < 4.0) {
      Ele.push_back(ielectron);
    }
  }

  // Loop over photons and pick out the ones we want to keep
  std::vector<int> Pho;
  for (int iphoton = 0; iphoton < nPho_; ++iphoton) {
    if (phoEt_[iphoton] > 10.0 &&
        phoPos_[iphoton] == 0) {
      Pho.push_back(iphoton);
    }
  }

  // We're looking for one electron and one photon
  if (Ele.size() != 1 || Pho.size() != 1) {
    return;
  }

  // Require the MET be above a certain value
  if (tcMET_ < 20.0) {
    return;
  }

  // Easier index for electron and photon
  int const ie = Ele[0];
  int const ip = Pho[0];

  // Just because I want some 4-vectors
  TLorentzVector ve;
  ve.SetPtEtaPhiE(elePt_[ie], eleEta_[ie], elePhi_[ie], eleEn_[ie]);
  TLorentzVector vp;
  vp.SetPtEtaPhiE(phoEt_[ip], phoEta_[ip], phoPhi_[ip], phoE_[ip]);

  // Simple print of this event
  printf("elePt: %6.2f  phoEt: %6.2f  tcMET: %6.2f\n", elePt_[ie], phoEt_[ip], tcMET_);

  // Fill some histograms for these nicely selected events.
  Hist.FillTH1D("tcMET", 20, 0, 200, tcMET_);
  Hist.FillTH1D("elePt", 20, 0, 200, elePt_[ie]);
  Hist.FillTH1D("phoEt", 20, 0, 200, phoEt_[ip]);
  Hist.FillTH1D("ElePhoDeltaR", 20, 0, 8, ve.DeltaR(vp));
  Hist.FillTH1D("nJet", 10, 0, 10, nJet_);

  return;
}
