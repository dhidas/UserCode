#include "DHidasAna/Dtuple/interface/SimpleSUSYPlots.h"

#include <iostream>

#include "TMath.h"

SimpleSUSYPlots::SimpleSUSYPlots (TString const& InFile) :
  DtupleReader(InFile)
{
  InitOutFile();
  InitializeHists();
}


SimpleSUSYPlots::SimpleSUSYPlots (std::vector<TString> const& InFiles) :
  DtupleReader(InFiles)
{
  InitOutFile();
  InitializeHists();
}


SimpleSUSYPlots::~SimpleSUSYPlots ()
{
  fOutFile->Write();
  fOutFile->Close();
}


void SimpleSUSYPlots::InitOutFile (TString const& OutFileName)
{
  // Initialize the output root file
  fOutFile = new TFile(OutFileName, "recreate");
  if (!fOutFile) {
    std::cerr << "ERROR: cannot open output root file" << std::endl;
    exit(1);
  }

  return;
}


void SimpleSUSYPlots::InitializeHists ()
{
  std::cout << "Initializing histograms" << std::endl;
  Hist1D["Electron_ECalIso"] = new TH1D("Electron_ECalIso", "Electron ECal Iso", 100, -1, 4);
  Hist1D["Electron_ECalIsoEE"] = new TH1D("Electron_ECalIsoEE", "Electron ECal Iso EE", 100, -1, 4);
  Hist1D["Electron_ECalIsoEB"] = new TH1D("Electron_ECalIsoEB", "Electron ECal Iso EB", 100, -1, 4);
  Hist1D["Electron_HCalIso"] = new TH1D("Electron_HCalIso", "Electron HCal Iso", 100, -1, 4);
  Hist1D["Electron_TrkIso"] = new TH1D("Electron_TrkIso", "Electron Trk Iso", 100, -1, 4);
  Hist1D["Electron_IsConvertedPhoton"] = new TH1D("Electron_IsConvertedPhoton", "Electron IsConvertedPhoton", 2, 0, 2);
  Hist1D["Muon_ECalIso"] = new TH1D("Muon_ECalIso", "Muon ECal Iso", 100, -1, 4);
  Hist1D["Muon_HCalIso"] = new TH1D("Muon_HCalIso", "Muon HCal Iso", 100, -1, 4);
  Hist1D["Muon_TrkIso"] = new TH1D("Muon_TrkIso", "Muon Trk Iso", 100, -1, 4);
  Hist1D["Muon_IsConvertedPhoton"] = new TH1D("Muon_IsConvertedPhoton", "Muon IsConvertedPhoton", 2, 0, 2);
  for (std::map<TString, TH1D*>::iterator It = Hist1D.begin(); It != Hist1D.end(); ++It) {
    It->second->SetDirectory(fOutFile);
    It->second->Sumw2();
  }

  return;
}

void SimpleSUSYPlots::Loop ()
{
  // Loop over all events in the Dtuple

  // Shorthand
  DtupleReader::Event_Struct& Ev = fEvent;


  for (long unsigned int ientry = 0; GetEntry(ientry) != 0; ++ientry) {
    if (ientry % 1000 == 0) {
      printf("Processing event: %15lu\n", ientry);
    }

    SortOutOverlaps(Ev);


    // Loop over all leptons
    for (int ilepton = 0; ilepton < Ev.NLeptons; ++ilepton) {
      if (Ev.Lepton_Flavor[ilepton] == Dtuple::kLeptonFlavor_Electron) {
        Hist1D["Electron_ECalIso"]->Fill(Ev.Lepton_ECalIso[ilepton]);

        if ( Ev.Lepton_Detector[ilepton] & ( 0x1 << Dtuple::kElectronDet_EE ) ) {
          Hist1D["Electron_ECalIsoEE"]->Fill(Ev.Lepton_ECalIso[ilepton]);
        } else if ( Ev.Lepton_Detector[ilepton] & ( 0x1 << Dtuple::kElectronDet_EB ) ) {
          Hist1D["Electron_ECalIsoEB"]->Fill(Ev.Lepton_ECalIso[ilepton]);
        }
        Hist1D["Electron_HCalIso"]->Fill(Ev.Lepton_HCalIso[ilepton]);
        Hist1D["Electron_TrkIso"]->Fill(Ev.Lepton_TrkIso[ilepton]);
        Hist1D["Electron_IsConvertedPhoton"]->Fill(Ev.Lepton_IsConvertedPhoton[ilepton]);
      } else if (Ev.Lepton_Flavor[ilepton] == Dtuple::kLeptonFlavor_Muon) {
        Hist1D["Muon_ECalIso"]->Fill(Ev.Lepton_ECalIso[ilepton]);
        Hist1D["Muon_HCalIso"]->Fill(Ev.Lepton_HCalIso[ilepton]);
        Hist1D["Muon_TrkIso"]->Fill(Ev.Lepton_TrkIso[ilepton]);
        Hist1D["Muon_IsConvertedPhoton"]->Fill(Ev.Lepton_IsConvertedPhoton[ilepton]);
      }
    }




  }

  return;
}



void SimpleSUSYPlots::SortOutOverlaps (DtupleReader::Event_Struct& Ev)
{
  DtupleReader::Event_Struct NewEv;

  // Quick check of lepton and jet numbers
  if (Ev.NLeptons > Dtuple::NMaxLeptons) {
    std::cerr << "WARNING: NLeptons > Dtuple::NMaxLeptons.  You are definitely missing some leptons." << std::endl;
  }
  if (Ev.NJets > Dtuple::NMaxJets) {
    std::cerr << "WARNING: NJets > Dtuple::NMaxJets.  You are definitely missing some jets." << std::endl;
  }

  std::vector<int> Muons;
  for (int i = 0; i < Ev.NLeptons; ++i) {
    if (Ev.Lepton_Flavor[i] == Dtuple::kLeptonFlavor_Muon) {
      Muons.push_back(i);
    }
  }

  std::vector<int> Electrons;
  for (int i = 0; i < Ev.NLeptons; ++i) {
    if (Ev.Lepton_Flavor[i] == Dtuple::kLeptonFlavor_Electron) {
      bool AcceptElectron = true;
      for (size_t iMu = 0; iMu != Muons.size(); ++iMu) {
        if ( TMath::Sqrt( TMath::Power( Ev.Lepton_Eta[iMu] - Ev.Lepton_Eta[i], 2) + TMath::Power( Ev.Lepton_Phi[iMu] - Ev.Lepton_Phi[i], 2) ) < 0.4 ) {
          AcceptElectron = false;
        }
      }
      if (AcceptElectron) {
        Electrons.push_back(i);
      }

    }
  }

  std::vector<int> Jets;
  for (int i = 0; i < Ev.NJets; ++i) {
    bool AcceptJet = true;
    for (size_t iLep = 0; iLep != Ev.NLeptons; ++iLep) {
      if ( TMath::Sqrt( TMath::Power( Ev.Lepton_Eta[iLep] - Ev.Jet_Eta[i], 2) + TMath::Power( Ev.Lepton_Phi[iLep] - Ev.Jet_Phi[i], 2) ) < 0.4 ) {
        AcceptJet = false;
      }
    }
    if (AcceptJet) {
      Jets.push_back(i);
    }
  }

  //printf("%4i %4i %4i %4i\n", (int) Ev.NLeptons, (int) Muons.size() + (int) Electrons.size(), (int) Ev.NJets, (int) Jets.size());


  return;
}


