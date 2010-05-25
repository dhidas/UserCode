#include "DHidasAna/SimpleAna/interface/SimpleAna.h"
#include "DHidasAna/Dtuple/interface/TAnaHist.h"

SimpleAna::SimpleAna (TString const ProcName, TChain* Chain) : TDtupleReader(Chain)
{
  // Constructor
  fProcName = ProcName;
  fOutFile = new TFile("SimpleAna_"+ProcName+".root", "recreate");
  if (!fOutFile->IsOpen()) {
    std::cerr << "ERROR: cannot open output file" << std::endl;
    throw;
  }

}



SimpleAna::~SimpleAna ()
{
  // Destructor!
  fOutFile->Write();
  fOutFile->Close();
  delete fOutFile;

  if (fRunFakes) {
    fFakeDtuple->Write();
    delete fFakeDtuple;
  }

}


void SimpleAna::BeginJob ()
{
  // Do things you want to do here at the beginning..
  std::cout << "SimpleAna::BeginJob() called" << std::endl;

  // If you're running fakes make a new dtuple for the fakes
  if (fRunFakes) {
    std::cout << "Running fakes!" << std::endl;
    fFakeDtuple = new TDtuple( new TFile("Fakes_"+fProcName+".root", "recreate") );
  }


  return;
}



void SimpleAna::SetFakeRateFile (TString const FileName)
{
  // Open the root file with the fakerates
  fFakeRate = new TFakeRate(FileName);
  return;
}


void SimpleAna::Analyze (long unsigned int const ientry)
{
  // This is the main function which is called once per event
  // Add analysis functions here.

  // Plot some electron ID variables
  PlotElectronId();

  // Fully select objects and remove overlaps
  Selection();
  ObjectCleaning();

  // Plot a few quantities of interest
  PlotEventQuantities();
  PlotLeptons();
  PlotPhotons();
  PlotJets();

  return;
}



void SimpleAna::PlotEventQuantities ()
{
  // Just plot a few basic event quantities
  static TAnaHist Hist(fOutFile, "PlotEventQuantities");

  Hist.FillTH1D("SumEt", 50, 0, 200, GetSumEt());
  Hist.FillTH1D("Met", 50, 0, 200, GetMet());
  Hist.FillTH1D("MetPhi", 50, 0, 2*TMath::Pi(), GetMetPhi());



  return;
}




void SimpleAna::PlotLeptons ()
{
  static TAnaHist Hist(fOutFile, "PlotLeptons");

  Hist.FillTH1D("NLeptons", 10, 0, 10, Leptons.size());

  int i = 0;
  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {

    TString const Flavors = GetLeptonFlavorsString(Leptons);
    Hist.FillTH1D("LeptonTrkIso_"+Lep->GetFlavorString(), 100, -0.1, 0.5, Lep->GetTrkIso() / Lep->Perp());
    Hist.FillTH1D("LeptonCalIso_"+Lep->GetFlavorString(), 100, -0.1, 0.5, Lep->GetCalIso() / Lep->Perp());
    Hist.FillTH1D("LeptonECalIsoDep_"+Lep->GetFlavorString(), 100, -0.1, 8.5, Lep->GetECalIsoDep() / Lep->Perp());
    Hist.FillTH1D("LeptonHCalIsoDep_"+Lep->GetFlavorString(), 100, -0.1, 8.5, Lep->GetHCalIsoDep() / Lep->Perp());

    char Name[100];
    Hist.FillTH1D("LeptonPt", 100, 0, 200, Lep->Perp());
    sprintf(Name, "LeptonPt%i", (int) i);
    Hist.FillTH1D(Name, 100, 0, 200, Lep->Perp());

    Hist.FillTH1D("LeptonEta", 50, -3, 3, Lep->Eta());
    sprintf(Name, "LeptonEta%i", (int) i);
    Hist.FillTH1D(Name, 50, -3, 3, Lep->Eta());

    Hist.FillTH1D("LeptonPhi", 50, -1.0*TMath::Pi(), TMath::Pi(), Lep->Phi());
    sprintf(Name, "LeptonPhi%i", (int) i);
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Lep->Phi());

  }

  return;
}




void SimpleAna::PlotPhotons ()
{
  static TAnaHist Hist(fOutFile, "PlotPhotons");

  Hist.FillTH1D("NPhotons", 10, 0, 10, Photons.size());

  for (size_t i = 0; i != Photons.size(); ++i) {
    char Name[100];

    Hist.FillTH1D("PhotonPt", 100, 0, 200, Photons[i].Perp());
    sprintf(Name, "PhotonPt%i", (int) i);
    Hist.FillTH1D(Name, 100, 0, 200, Photons[i].Perp());

    Hist.FillTH1D("PhotonEta", 50, -3, 3, Photons[i].Eta());
    sprintf(Name, "PhotonEta%i", (int) i);
    Hist.FillTH1D(Name, 50, -3, 3, Photons[i].Eta());

    Hist.FillTH1D("PhotonPhi", 50, -1.0*TMath::Pi(), TMath::Pi(), Photons[i].Phi());
    sprintf(Name, "PhotonPhi%i", (int) i);
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Photons[i].Phi());
  }
  return;
}





void SimpleAna::PlotJets ()
{
  static TAnaHist Hist(fOutFile, "PlotJets");

  Hist.FillTH1D("NJets", 20, 0, 20, Jets.size());

  for (size_t i = 0; i != Jets.size(); ++i) {
    Hist.FillTH1D("JetPt", 100, 0, 200, Jets[i].Perp());

    char Name[100];
    sprintf(Name, "JetPt%i", (int) i);
    Hist.FillTH1D(Name, 100, 0, 200, Jets[i].Perp());

    Hist.FillTH1D("JetEta", 50, -3, 3, Jets[i].Eta());
    sprintf(Name, "JetEta%i", (int) i);
    Hist.FillTH1D(Name, 50, -3, 3, Jets[i].Eta());

    Hist.FillTH1D("JetPhi", 50, -1.0*TMath::Pi(), TMath::Pi(), Jets[i].Phi());
    sprintf(Name, "JetPhi%i", (int) i);
    Hist.FillTH1D(Name, 50, -1.0*TMath::Pi(), TMath::Pi(), Jets[i].Phi());
  }

  if (Jets.size() >= 2) {
    Hist.FillTH1D("JetDeltaR01", 100, 0, 8, Jets[0].DeltaR(Jets[1]));
    Hist.FillTH1D("JetDeltaEta01", 100, 0, 8, TMath::Abs(Jets[0].Eta() - Jets[1].Eta()));
    Hist.FillTH1D("JetDeltaPhi01", 100, -TMath::Pi(), TMath::Pi(), Jets[0].DeltaPhi(Jets[1]));
  }

  if (Jets.size() == 2) {
    Hist.FillTH1D("DiJetMass", 100, 0, 300, (Jets[0] + Jets[1]).M());
    Hist.FillTH1D("DiJetSumPt", 100, 0, 400, Jets[0].Perp() + Jets[1].Perp());
    Hist.FillTH1D("DiJetDeltaR01", 100, 0, 8, Jets[0].DeltaR(Jets[1]));
    Hist.FillTH1D("DiJetDeltaEta01", 100, 0, 8, TMath::Abs(Jets[0].Eta() - Jets[1].Eta()));
    Hist.FillTH1D("DiJetDeltaPhi01", 100, -TMath::Pi(), TMath::Pi(), Jets[0].DeltaPhi(Jets[1]));
  }

  if (Jets.size() == 3) {
    Hist.FillTH1D("TriJetMass", 100, 0, 300, (Jets[0] + Jets[1] + Jets[2]).M());
    Hist.FillTH1D("DiJetSumPt", 100, 0, 400, Jets[0].Perp() + Jets[1].Perp() + Jets[1].Perp());
    Hist.FillTH1D("TriJetDeltaR01", 100, 0, 8, Jets[0].DeltaR(Jets[1]));
    Hist.FillTH1D("TriJetDeltaEta01", 100, 0, 8, TMath::Abs(Jets[0].Eta() - Jets[1].Eta()));
    Hist.FillTH1D("TriJetDeltaPhi01", 100, -TMath::Pi(), TMath::Pi(), Jets[0].DeltaPhi(Jets[1]));
  }

  return;
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
      bool Keep = PassSelectionMuon(*lep);
      if (Keep) {
        NewLeptons.push_back(*lep);
      }
    } else if (lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      bool Keep = PassSelectionElectron(*lep);
      if (Keep) {
        NewLeptons.push_back(*lep);
      }
    }
  }

  Leptons = NewLeptons;

  return;
}



void SimpleAna::SelectionPhoton ()
{
  std::vector<TPhoton> NewPhotons;

  for (std::vector<TPhoton>::iterator photon = Photons.begin(); photon != Photons.end(); ++photon) {
    bool Keep = true;
    if (TMath::Abs(photon->Eta()) > 2.5) Keep = false;
    if (photon->Perp() < 15.0) Keep = false;
    if (photon->GetCalIso() / photon->Perp() > 0.1) Keep = false;
    if (Keep) {
      NewPhotons.push_back(*photon);
    }
  }

  Photons = NewPhotons;

  return;
}


bool SimpleAna::PassSelectionElectron (TLepton& Lep)
{
  bool Pass = true;
  if (!Lep.PassesSelection(TLepton::kElectronSel_RobustTight)) Pass = false;
  if ( (TMath::Abs(Lep.Eta()) > 2.1) ) Pass = false;
  if ( Lep.GetCalIso() / Lep.Perp() > 0.15) Pass = false;
  if ( Lep.Perp() < 8.0) Pass = false;
  if ( TMath::Abs(Lep.Getdxy()) > 0.02) Pass = false;
  return Pass || Lep.GetIsFake();
}


bool SimpleAna::PassSelectionMuon (TLepton& Lep)
{
  bool Pass = true;
  if ( !Lep.PassesSelection(TLepton::kMuonSel_GlobalMuonPromptTight)) Pass = false;
  if ( (TMath::Abs(Lep.Eta()) > 2.1) ) Pass = false;
  if ( Lep.GetCalIso() / Lep.Perp() > 0.1) Pass = false;
  return Pass || Lep.GetIsFake();
}


void SimpleAna::PlotElectronId ()
{
  static TAnaHist Hist(fOutFile, "PlotElectronId");

  TString Det;

  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {
    if (Lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      if (PassSelectionElectron(*Lep)) {
        Hist.FillTH1D("Pt", 50, 0, 200, Lep->Perp());
        Hist.FillTH1D("Eta", 50, -3, 3, Lep->Eta());
        Hist.FillTH1D("Phi", 50, -TMath::Pi(), TMath::Pi(), Lep->Phi());
        Hist.FillTH1D("HOverE", 50, 0, 0.02, Lep->GetHCalOverECal());
        Hist.FillTH1D("SigmaIEtaIEta", 50, 0, 0.04, Lep->GetSigmaIEtaIEta());
        Hist.FillTH1D("DeltaPhiIn", 50, -0.1, 0.1, Lep->GetDeltaPhiIn());
        Hist.FillTH1D("DeltaEtaIn", 50, -0.02, 0.02, Lep->GetDeltaEtaIn());
        Hist.FillTH1D("e25Maxoe55", 50, 0, 1, Lep->GetE2x5overE5x5());
        //Hist.FillTH1D("e15oe55", 50, 0, 200, Lep->Gete15oe55());


        if (Lep->IsElectronDet(TLepton::kElectronDet_EB)) {
          Det = "_EB";
        } else if (Lep->IsElectronDet(TLepton::kElectronDet_EE)) {
          Det = "_EE";
        } else {
          Det = "_Gap";
        }
        Hist.FillTH1D("Pt"+Det, 50, 0, 200, Lep->Perp());
        Hist.FillTH1D("Eta"+Det, 50, -3, 3, Lep->Eta());
        Hist.FillTH1D("Phi"+Det, 50, -TMath::Pi(), TMath::Pi(), Lep->Phi());
        Hist.FillTH1D("HOverE"+Det, 50, 0, 0.02, Lep->GetHCalOverECal());
        Hist.FillTH1D("SigmaIEtaIEta"+Det, 50, 0, 0.04, Lep->GetSigmaIEtaIEta());
        Hist.FillTH1D("DeltaPhiIn"+Det, 50, -0.1, 0.1, Lep->GetDeltaPhiIn());
        Hist.FillTH1D("DeltaEtaIn"+Det, 50, -0.02, 0.02, Lep->GetDeltaEtaIn());
        Hist.FillTH1D("e25Maxoe55"+Det, 50, 0, 1, Lep->GetE2x5overE5x5());
        //Hist.FillTH1D("e15oe55"+Det, 50, 0, 200, Lep->Gete15oe55());
      }
    }
  }

  return;
}



void SimpleAna::SelectionJet ()
{
  std::vector<TJet> NewJets;

  for (std::vector<TJet>::iterator jet = Jets.begin(); jet != Jets.end(); ++jet) {
    bool Keep = true;
    if (jet->Perp() < 15) {
      Keep = false;
    }
    if ( TMath::Abs(jet->Eta()) > 2.5) {
      Keep = false;
    }

    if (Keep) {
      NewJets.push_back(*jet);
    }
  }

  Jets = NewJets;

  return;
}




void SimpleAna::RunFakes (bool const in)
{
  fRunFakes = in;
  return;
}


bool SimpleAna::RunFakes ()
{
  return fRunFakes;
}






void SimpleAna::EndJob ()
{
  // Do things you want to do here at the end..
  std::cout << "SimpleAna::EndJob() called" << std::endl;

  return;
}
