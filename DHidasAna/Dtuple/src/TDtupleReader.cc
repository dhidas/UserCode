#include "DHidasAna/Dtuple/interface/TDtupleReader.h"


TDtupleReader::TDtupleReader (TChain* Chain) : TDtuple(Chain)
{
}


TDtupleReader::~TDtupleReader ()
{
}


void TDtupleReader::Loop (long unsigned int Max)
{

  unsigned int const ReportEvery = 10000;
  if (Max == 0) {
    for (long unsigned int ientry = 0; GetEntry(ientry); ++ientry) {
      if (ientry % ReportEvery == 0) {
        std::cout << "Events Processed: " << ientry << std::endl;
      }
      Analyze(ientry);
    }
  } else {
    for (long unsigned int ientry = 0; GetEntry(ientry) && (ientry < Max); ++ientry) {
      if (ientry % ReportEvery == 0) {
        std::cout << "Events Processed: " << ientry << std::endl;
      }
      Analyze(ientry);
    }
  }

  return;
}




void TDtupleReader::ObjectCleaning ()
{
  std::vector<TLorentzVector> GoodObjects;

  std::vector<TLepton> NewLeptons;
  for (size_t i = 0; i != Leptons.size(); ++i) {
    if (Leptons[i].IsFlavor(TLepton::kLeptonFlavor_Muon)) {
      NewLeptons.push_back(Leptons[i]);
      GoodObjects.push_back(Leptons[i]);
    }
  }

  for (size_t i = 0; i != Leptons.size(); ++i) {
    bool Keep = true;
    if (Leptons[i].IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      for (size_t j = 0; j != NewLeptons.size(); ++j) {
        if (Leptons[i].DeltaR(NewLeptons[j]) < 0.3) {
          Keep = false;
        }
      }
    } else {
      Keep = false;
    }
    if (Keep) {
      NewLeptons.push_back(Leptons[i]);
      GoodObjects.push_back(Leptons[i]);
    }
  }
  Leptons = NewLeptons;

  std::vector<TPhoton> NewPhotons;
  for (size_t i = 0; i != Photons.size(); ++i) {
    bool Keep = true;
    for (size_t j = 0; j != Leptons.size(); ++j) {
      if (Photons[i].DeltaR( Leptons[j] ) < 0.3) {
        Keep = false;
      }
    }
    if (Keep) {
      NewPhotons.push_back(Photons[i]);
      GoodObjects.push_back(Photons[i]);
    }
  }
  Photons = NewPhotons;

  std::vector<TJet> NewJets;
  for (size_t i = 0; i != Jets.size(); ++i) {
    bool Keep = true;
    for (size_t j = 0; j != GoodObjects.size(); ++j) {
      if (Jets[i].DeltaR( GoodObjects[j] ) < 0.3) {
        Keep = false;
      }
    }
    if (Keep) {
      NewJets.push_back(Jets[i]);
    }
  }
  Jets = NewJets;

  return;
}
