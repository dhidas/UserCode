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

  std::vector<int> Leptons;
  Leptons.insert(Leptons.end(), Electrons.begin(), Electrons.end());
  Leptons.insert(Leptons.end(), Muons.begin(), Muons.end());
  std::sort(Leptons.begin(), Leptons.end());
  std::reverse(Leptons.begin(), Leptons.end());

  DtupleReader::Event_Struct NewEv;
  DefaultValues(NewEv);
  CopyILeptonFromTo(1, Ev, NewEv);
  //Ev = NewEv;
  //printf("%4i %4i %4i %4i\n", (int) Ev.NLeptons, (int) Muons.size() + (int) Electrons.size(), (int) Ev.NJets, (int) Jets.size());


  return;
}



void SimpleSUSYPlots::CopyILeptonFromTo (int const i, Dtuple::Event_Struct& From, Dtuple::Event_Struct& To)
{
  if (i >= Dtuple::NMaxLeptons || To.NLeptons >= Dtuple::NMaxLeptons || i >= From.NLeptons) {
    std::cerr << "ERROR: you're over the lepton limit.  There is something wrong with your copying" << std::endl;
    std::cout << "  i=" << i << " From.NLeptons=" << From.NLeptons << " To.NLeptons=" << To.NLeptons << " Max=" << Dtuple::NMaxLeptons << std::endl;
    return;
  }

  int const iTo = To.NLeptons;
  ++To.NLeptons;

  To.Lepton_Px[iTo]                 = From.Lepton_Px[i];
  To.Lepton_Py[iTo]                 = From.Lepton_Py[i];
  To.Lepton_Pz[iTo]                 = From.Lepton_Pz[i];
  To.Lepton_Pt[iTo]                 = From.Lepton_Pt[i];
  To.Lepton_TrkPt[iTo]              = From.Lepton_TrkPt[i];
  To.Lepton_Eta[iTo]                = From.Lepton_Eta[i];
  To.Lepton_Phi[iTo]                = From.Lepton_Phi[i];
  To.Lepton_dxy[iTo]                = From.Lepton_dxy[i];
  To.Lepton_dz[iTo]                 = From.Lepton_dz[i];
  To.Lepton_Z0[iTo]                 = From.Lepton_Z0[i];
  To.Lepton_Charge[iTo]             = From.Lepton_Charge[i];
  To.Lepton_Flavor[iTo]             = From.Lepton_Flavor[i];
  To.Lepton_TrkIso[iTo]             = From.Lepton_TrkIso[i];
  To.Lepton_CalIso[iTo]             = From.Lepton_CalIso[i];
  To.Lepton_ECalIso[iTo]            = From.Lepton_ECalIso[i];
  To.Lepton_HCalIso[iTo]            = From.Lepton_HCalIso[i];
  To.Lepton_CalE[iTo]               = From.Lepton_CalE[i];
  To.Lepton_HCalOverECal[iTo]       = From.Lepton_HCalOverECal[i];
  To.Lepton_EoverPin[iTo]           = From.Lepton_EoverPin[i];
  To.Lepton_fBrem[iTo]              = From.Lepton_fBrem[i];
  To.Lepton_IsConvertedPhoton[iTo]  = From.Lepton_IsConvertedPhoton[i];
  To.Lepton_PassSelection[iTo]      = From.Lepton_PassSelection[i];
  To.Lepton_Detector[iTo]           = From.Lepton_Detector[i];
  To.Lepton_Classification[iTo]     = From.Lepton_Classification[i];
  To.Lepton_SigmaIEtaIEta[iTo]      = From.Lepton_SigmaIEtaIEta[i];
  To.Lepton_DeltaEtaIn[iTo]         = From.Lepton_DeltaEtaIn[i];
  To.Lepton_DeltaPhiIn[iTo]         = From.Lepton_DeltaPhiIn[i];
  To.Lepton_E2x5overE5x5[iTo]       = From.Lepton_E2x5overE5x5[i];

  return;
}
