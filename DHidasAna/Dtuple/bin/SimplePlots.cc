////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Mon Nov  9 20:05:50 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "DHidasAna/Dtuple/interface/DtupleReader.h"

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

int SimplePlots (TString const FileName)
{
  DtupleReader DR(FileName);

  DtupleReader::Event_Struct& Ev = DR.fEvent;

  TFile OutFile("SampleOutFile.root", "recreate");
  TH1D DeltaREGamma("DeltaREGamma", "DeltaREGamma", 50, 0, 4);
  DeltaREGamma.SetDirectory(&OutFile);

  for (unsigned long ientry = 0; DR.GetEntry(ientry) != 0; ++ientry) {
    if (ientry % 1000 == 0) {
      printf("Processing event: %15lu\n", ientry);
    }

    if (Ev.NLeptons == 1 && Ev.NPhotons == 1) {
      DeltaREGamma.Fill( TMath::Sqrt( TMath::Power(Ev.Lepton_Phi[0] - Ev.Photon_Phi[0], 2) + TMath::Power(Ev.Lepton_Eta[0] - Ev.Photon_Eta[0], 2) ) );
    }
  }

  OutFile.Write();
  OutFile.Close();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InputFile.root]" << std::endl;
    return 1;
  }

  SimplePlots(argv[1]);

  return 0;
}
