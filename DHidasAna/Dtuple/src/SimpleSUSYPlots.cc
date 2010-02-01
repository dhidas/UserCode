#include "DHidasAna/Dtuple/interface/SimpleSUSYPlots.h"

#include <iostream>

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
  Hist1D["Electron_ECalIso"] = new TH1D("Electron_ECalIso", "Electron ECal Iso", 100, 0, 5);
  Hist1D["Electron_HCalIso"] = new TH1D("Electron_HCalIso", "Electron HCal Iso", 100, 0, 5);
  Hist1D["Electron_TrkIso"] = new TH1D("Electron_TrkIso", "Electron Trk Iso", 100, 0, 5);
  Hist1D["Muon_ECalIso"] = new TH1D("Muon_ECalIso", "Muon ECal Iso", 100, 0, 5);
  Hist1D["Muon_HCalIso"] = new TH1D("Muon_HCalIso", "Muon HCal Iso", 100, 0, 5);
  Hist1D["Muon_TrkIso"] = new TH1D("Muon_TrkIso", "Muon Trk Iso", 100, 0, 5);
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


    // Loop over all leptons
    for (int ilepton = 0; ilepton < Ev.NLeptons; ++ilepton) {
      if (Ev.Lepton_Flavor[ilepton] == Dtuple::kElectron) {
        Hist1D["Electron_ECalIso"]->Fill(Ev.Lepton_ECalIso[ilepton]);
        Hist1D["Electron_HCalIso"]->Fill(Ev.Lepton_HCalIso[ilepton]);
        Hist1D["Electron_TrkIso"]->Fill(Ev.Lepton_TrkIso[ilepton]);
      } else if (Ev.Lepton_Flavor[ilepton] == Dtuple::kMuon) {
        Hist1D["Muon_ECalIso"]->Fill(Ev.Lepton_ECalIso[ilepton]);
        Hist1D["Muon_HCalIso"]->Fill(Ev.Lepton_HCalIso[ilepton]);
        Hist1D["Muon_TrkIso"]->Fill(Ev.Lepton_TrkIso[ilepton]);
      }
    }




  }

  return;
}
