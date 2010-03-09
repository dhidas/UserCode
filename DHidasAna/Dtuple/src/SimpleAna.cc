#include "DHidasAna/Dtuple/interface/SimpleAna.h"
#include "DHidasAna/Dtuple/interface/TAnaHist.h"

SimpleAna::SimpleAna (TString const ProcName, TChain* Chain) : TDtupleReader(Chain)
{
  fProcName = ProcName;
  fOutFile = new TFile("SimpleAna_"+ProcName+".root", "recreate");
  if (!fOutFile->IsOpen()) {
    std::cerr << "ERROR: cannot open output file" << std::endl;
    exit(1);
  }
}



SimpleAna::~SimpleAna ()
{
  fOutFile->Write();
  fOutFile->Close();
  delete fOutFile;
}


void SimpleAna::Analyze (long unsigned int const ientry)
{
  Selection();
  ObjectCleaning();

  PlotLeptons();
  PlotPhotons();
  PlotJets();

  PlotDileptonMass ();
  PlotZllE();

  return;
}



void SimpleAna::PlotLeptons ()
{
  static TAnaHist Hist(fOutFile, "PlotLeptons");

  Hist.FillTH1D("NLeptons_"+fProcName, 10, 0, 10, Leptons.size());

  for (size_t i = 0; i != Leptons.size(); ++i) {
    char Name[100];

    Hist.FillTH1D("LeptonPt_"+fProcName, 100, 0, 200, Leptons[i].Perp());
    sprintf(Name, "LeptonPt%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 100, 0, 200, Leptons[i].Perp());

    Hist.FillTH1D("LeptonEta_"+fProcName, 50, -3, 3, Leptons[i].Eta());
    sprintf(Name, "LeptonEta%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -3, 3, Leptons[i].Eta());

    Hist.FillTH1D("LeptonPhi_"+fProcName, 50, -1.0*TMath::Pi(), TMath::Pi(), Leptons[i].Phi());
    sprintf(Name, "LeptonPhi%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Leptons[i].Phi());
  }
  return;
}




void SimpleAna::PlotPhotons ()
{
  static TAnaHist Hist(fOutFile, "PlotPhotons");

  Hist.FillTH1D("NPhotons_"+fProcName, 10, 0, 10, Photons.size());

  for (size_t i = 0; i != Photons.size(); ++i) {
    char Name[100];

    Hist.FillTH1D("PhotonPt_"+fProcName, 100, 0, 200, Photons[i].Perp());
    sprintf(Name, "PhotonPt%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 100, 0, 200, Photons[i].Perp());

    Hist.FillTH1D("PhotonEta_"+fProcName, 50, -3, 3, Photons[i].Eta());
    sprintf(Name, "PhotonEta%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -3, 3, Photons[i].Eta());

    Hist.FillTH1D("PhotonPhi_"+fProcName, 50, -1.0*TMath::Pi(), TMath::Pi(), Photons[i].Phi());
    sprintf(Name, "PhotonPhi%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Photons[i].Phi());
  }
  return;
}





void SimpleAna::PlotJets ()
{
  static TAnaHist Hist(fOutFile, "PlotJets");

  Hist.FillTH1D("NJets_"+fProcName, 20, 0, 20, Jets.size());

  for (size_t i = 0; i != Jets.size(); ++i) {
    Hist.FillTH1D("JetPt_"+fProcName, 100, 0, 200, Jets[i].Perp());

    char Name[100];
    sprintf(Name, "JetPt%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 100, 0, 200, Jets[i].Perp());

    Hist.FillTH1D("JetEta_"+fProcName, 50, -3, 3, Jets[i].Eta());
    sprintf(Name, "JetEta%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -3, 3, Jets[i].Eta());

    Hist.FillTH1D("JetPhi_"+fProcName, 50, -1.0*TMath::Pi(), TMath::Pi(), Jets[i].Phi());
    sprintf(Name, "JetPhi%i_%s", (int) i, fProcName.Data());
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Jets[i].Phi());
  }
  return;
}



void SimpleAna::PlotZllE ()
{
  static TAnaHist Hist(fOutFile, "PlotZllE");

  if (Leptons.size() != 3) {
    return;
  }

  std::vector<TLepton> Zll = ClosestZMatch(Leptons, false, false);
  if (Zll.size() == 2) {
    std::cout << (Zll[0] + Zll[1]).M() << std::endl;
  }

  return;
}



void SimpleAna::PlotDileptonMass ()
{
  static TAnaHist Hist(fOutFile, "PlotDileptonMass");

  if (Leptons.size() != 2) {
    return;
  }

  TString Charges;
  if (Leptons[0].GetCharge() != Leptons[1].GetCharge()) {
    Charges = "OS";
  } else {
    Charges = "SS";
  }

  TString const Flavors = GetLeptonFlavorsString(Leptons);

  Hist.FillTH1D("DileptonMass"+Charges+"_"+Flavors+"_"+fProcName, 100, 0, 200, (Leptons[0]+Leptons[1]).M());

  return;
}



std::vector<TLepton> SimpleAna::ClosestZMatch(std::vector<TLepton>& InLeptons, bool const RequireOS, bool const RequireSF)
{
  std::vector<TLepton> Zll;
  if (Leptons.size() < 2) {
    return Zll;
  }

  int BestZll[2];
  float BestDiffMass = 999999;

  for (size_t i = 0; i < InLeptons.size()-1; ++i) {
    for (size_t j = i+1; j < InLeptons.size(); ++j) {
      if (RequireOS && (InLeptons[i].GetCharge() + InLeptons[j].GetCharge() != 0)) {
        continue;
      }
      if (RequireSF && (InLeptons[i].GetFlavor() != InLeptons[j].GetFlavor())) {
        continue;
      }
        std::cout << "MassDiff = " << TMath::Abs((InLeptons[i] + InLeptons[j]).M() - 91.0) << std::endl;
      if ( TMath::Abs((InLeptons[i] + InLeptons[j]).M() - 91.0) < BestDiffMass ) {
        BestDiffMass = TMath::Abs((InLeptons[i] + InLeptons[j]).M() - 91.0);
        BestZll[0] = i;
        BestZll[1] = j;
      }
    }
  }

  return Zll;
}


void SimpleAna::Selection ()
{
  SelectionLepton();
  SelectionPhoton();
  SelectionJet();

  return;
}



void SimpleAna::SelectionLepton ()
{
  std::vector<TLepton> NewLeptons;

  for (std::vector<TLepton>::iterator lep = Leptons.begin(); lep != Leptons.end(); ++lep) {
    if (lep->IsFlavor(TLepton::kLeptonFlavor_Muon)) {
      NewLeptons.push_back(*lep);
    } else if (lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      NewLeptons.push_back(*lep);
    }
  }

  Leptons = NewLeptons;

  return;
}



void SimpleAna::SelectionPhoton ()
{
  return;
}



void SimpleAna::SelectionJet ()
{
  return;
}
