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


    // Plot anything you want down here!!

    std::vector<LHEParticle> Leptons;
    for (size_t ip = 0; ip != Particles.size(); ++ip) {
      LHEParticle P = Particles[ip];

      if (TMath::Abs( P.Id ) == 11 || TMath::Abs( P.Id ) == 13) {
        Leptons.push_back(P);
        Hist.FillTH1D("LeptonPt", "Lepton p_{T}", "Lepton p_{T}(GeV/c)", "# of Events", 20, 0, 200, P.Perp());
      }

      if (P.Id == 22) {
        Hist.FillTH1D("PhotonEt", "Photon E_{T}", "Photon E_{T}(GeV)", "# of Events", 20, 0, 200, P.Perp(), Weight);
        Hist.FillTH1D("PhotonEta", "Photon #Eta", "Photon #Eta", "# of Events", 20, -5, 5, P.Eta(), Weight);
      }

    }

    if (Leptons.size() == 2) {
      Hist.FillTH1D("DileptonMass", "Dilepton Mass", "Dilepton Mass (GeV/c^{2})", "# of Events", 20, 0, 200, (Leptons[0] + Leptons[1]).M(), Weight);
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



LHEPlotter::~LHEPlotter ()
{
  if (OutFile) {
    OutFile->Write();
    OutFile->Close();
  }
}

