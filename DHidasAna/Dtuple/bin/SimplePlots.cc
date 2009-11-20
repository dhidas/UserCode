////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Mon Nov  9 20:05:50 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <map>

#include "DHidasAna/Dtuple/interface/DtupleReader.h"

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

int SimplePlots (TString const FileName)
{
  // Dtuple reader and event structure
  DtupleReader DR(FileName);
  DtupleReader::Event_Struct& Ev = DR.fEvent;

  // Output file
  TFile OutFile("SampleOutFile.root", "recreate");

  // Define histograms and set them to outfile
  std::map<TString, TH1D*> Hist;
  Hist["DeltaREGamma"] = new TH1D ("DeltaREGamma", "DeltaREGamma;#Delta R (e,#gamma);Events / 0.08", 50, 0, 4);
  Hist["PhotonEt"] = new TH1D("PhotonEt", "PhotonEt;Photon E_T (GeV);Events / (2 GeV)", 50, 0, 100);
  Hist["PhotonTrkIso"] = new TH1D("PhotonTrkIso", "PhotonTrkIso;Photon Track Isolation (GeV);Number of Photons", 50, 0, 20);
  Hist["PhotonCalIso"] = new TH1D("PhotonCalIso", "PhotonCalIso;Photon Calorimeter Isolation (GeV);Number of Photons", 50, 0, 20);
  Hist["PhotonHCalOverECal"] = new TH1D("PhotonHCalOverECal", "PhotonHCalOverECal;Photon E_{HCal} / E_{ECal};Number of Photons", 50, 0, 1);
  for (std::map<TString, TH1D*>::iterator It = Hist.begin(); It != Hist.end(); ++It) {
    It->second->SetDirectory(&OutFile);
    It->second->Sumw2();
  }


  // Loop over all events
  for (unsigned long ientry = 0; DR.GetEntry(ientry) != 0; ++ientry) {
    if (ientry % 1000 == 0) {
      printf("Processing event: %15lu\n", ientry);
    }

    // If it's interesting plot some stuff
    if (Ev.NLeptons == 1 && Ev.NPhotons == 1) {
      Hist["DeltaREGamma"]->Fill( TMath::Sqrt( TMath::Power(Ev.Lepton_Phi[0] - Ev.Photon_Phi[0], 2) + TMath::Power(Ev.Lepton_Eta[0] - Ev.Photon_Eta[0], 2) ) );
      Hist["PhotonEt"]->Fill( Ev.Photon_Pt[0] );
      Hist["PhotonTrkIso"]->Fill( Ev.Photon_TrkIso[0] );
      Hist["PhotonCalIso"]->Fill( Ev.Photon_CalIso[0] );
      Hist["PhotonHCalOverECal"]->Fill( Ev.Photon_HCalOverECal[0] );
    }
  }

  // Write and close
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
