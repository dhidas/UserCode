////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 16:33:37 PDT 2010
//
////////////////////////////////////////////////////////////////////

#include "VgammaAna/LHEReader/interface/LHEPlotter.h"

#include "VgammaAna/LHEReader/interface/TAnaHist.h"

#include <iostream>







int LHEPlotter::Loop ()
{
  TAnaHist Hist(OutFile, "LHEPlots");

  unsigned int NEvents = 0;
  long double SumOfWeights = 0.0;
  while (NextEvent()) {
    ++NEvents;
    SumOfWeights += Weight;
    if (NEvents % 100 == 0) {
      printf("Processing event: %12i\n", NEvents);
    }



    int   const PtBins = 40;
    float const PtMin = 0;
    float const PtMax = 200;

    int   const EtaBins = 40;
    float const EtaMax = 5;

    // Plot anything you want down here!!
    std::vector<LHEParticle> Leptons;
    for (size_t ip = 0; ip != Particles.size(); ++ip) {
      LHEParticle P = Particles[ip];

      if (TMath::Abs( P.Id ) == 11 || TMath::Abs( P.Id ) == 13) {
        Leptons.push_back(P);
        Hist.FillTH1D("LeptonPt", "Lepton p_{T}", "Lepton p_{T}(GeV/c)", "# of Events", PtBins, PtMin, PtMax, P.Perp());
        Hist.FillTH1D("LeptonEta", "Lepton #eta", "Lepton #eta", "# of Events", EtaBins, -EtaMax, EtaMax, P.Eta(), Weight);
        Hist.FillTH1D("LeptonPhi", "Lepton #phi", "Lepton #phi", "# of Events", EtaBins, -TMath::Pi(), TMath::Pi(), P.Phi(), Weight);
        Hist.FillTH2D("LeptonEtaPhi", "Lepton #eta-#phi", "Lepton #eta", "Lepton #phi", EtaBins, -EtaMax, EtaMax, EtaBins, -TMath::Pi(), TMath::Pi(), P.Eta(), P.Phi(), Weight);
      }

      if (P.Id == 22) {
        Hist.FillTH1D("PhotonEt", "Photon E_{T}", "Photon E_{T}(GeV)", "# of Events", PtBins, PtMin, PtMax, P.Perp(), Weight);
        Hist.FillTH1D("PhotonEta", "Photon #eta", "Photon #eta", "# of Events", EtaBins, -EtaMax, EtaMax, P.Eta(), Weight);
        Hist.FillTH1D("PhotonPhi", "Photon #phi", "Photon #phi", "# of Events", EtaBins, -TMath::Pi(), TMath::Pi(), P.Phi(), Weight);
      }

      if (abs(P.Id) == 12 || abs(P.Id) == 14) {
        Hist.FillTH1D("NeutrinoEt", "Neutrino E_{T}", "Neutrino E_{T}(GeV)", "# of Events", PtBins, PtMin, PtMax, P.Perp(), Weight);
      }

    }

    if (Leptons.size() == 2) {
      Hist.FillTH1D("DileptonMass", "Dilepton Mass", "Dilepton Mass (GeV/c^{2})", "# of Events", PtBins, PtMin, PtMax, (Leptons[0] + Leptons[1]).M(), Weight);
    }



  }

  return NEvents;
}





LHEPlotter::LHEPlotter (TString const LHEFileName, TString const OutFileName) :
LHEEvent(LHEFileName)
{
  if (OutFileName.EndsWith("lhe") || OutFileName == LHEFileName) {
    std::cerr << "ERROR: this is just here to save you..." << std::endl;
    exit(1);
  }

  OutFile = new TFile(OutFileName, "recreate");
  if (!OutFile) {
    std::cerr << "ERROR: cannot open output file " << OutFileName << std::endl;
    exit(1);
  }
}



LHEPlotter::LHEPlotter (std::vector<TString> const& LHEFileNames, TString const OutFileName) :
LHEEvent(LHEFileNames)
{
  if (OutFileName.EndsWith("lhe")) {
    std::cerr << "ERROR: this is just here to save you..." << std::endl;
    exit(1);
  }

  OutFile = new TFile(OutFileName, "recreate");
  if (!OutFile) {
    std::cerr << "ERROR: cannot open output file " << OutFileName << std::endl;
    exit(1);
  }
}



LHEPlotter::~LHEPlotter ()
{
  if (OutFile) {
    OutFile->Write();
    OutFile->Close();
  }
}

