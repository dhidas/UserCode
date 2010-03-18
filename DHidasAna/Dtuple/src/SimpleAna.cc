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

  PlotDileptonMass();
  PlotTriLeptons();
  PlotZllE();
  Plot6JetEvents();
  PlotNLepSumEtJetsMet();

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



void SimpleAna::PlotTriLeptons ()
{
  static TAnaHist Hist(fOutFile, "PlotTriLeptons");

  if (Leptons.size() != 3) {
    return;
  }

  if (GetLeptonFlavorsString(Leptons) != "emm") {
    return;
  }

  std::vector<TLepton> Zll = ClosestZMatch(Leptons, false, false, true);
  Hist.FillTH1D("ZMass_"+GetLeptonFlavorsString(Leptons)+"_"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());

  TString const FlavorOrder = Zll[0].GetFlavorString() + Zll[1].GetFlavorString() + Zll[2].GetFlavorString();
  if (Zll[0].GetCharge() == Zll[1].GetCharge()) {
    Hist.FillTH1D("ZMass_FlavorsSS_"+FlavorOrder+"_"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
  } else {
    Hist.FillTH1D("ZMass_FlavorsOS_"+FlavorOrder+"_"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
  }



  return;
}



void SimpleAna::PlotZllE ()
{
  static TAnaHist Hist(fOutFile, "PlotZllE");

  if (Leptons.size() != 3) {
    return;
  }

  std::vector<TLepton> Zll = ClosestZMatch(Leptons, true, true, true);
  if (Zll.size() != 3) {
    return;
  }
  if (Zll[0].IsFlavor(TLepton::kLeptonFlavor_Muon) && Zll[2].IsFlavor(TLepton::kLeptonFlavor_Electron)) {
    Hist.FillTH1D("ZmmEConversion_"+fProcName, 2, 0, 2, Zll[2].GetIsConvertedPhoton());
    Hist.FillTH1D("ZmmMass_"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());

    float const DeltaRMin = std::min( Zll[2].DeltaR(Zll[0]), Zll[2].DeltaR(Zll[1]) );
    Hist.FillTH1D("DeltaRMinem_"+fProcName, 25, 0, 4, DeltaRMin);
  }

  // Plot for tagged conversion and non-tagged
  if (Zll[2].GetIsConvertedPhoton()) {
    Hist.FillTH1D("ZmmMass_TaggedConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
  } else {
    Hist.FillTH1D("ZmmMass_NoTagConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
  }

  // Match to the closest photon in DeltaR
  TGenP* MatchedPhoton = FindClosestGenP(Zll[2], 22);
  if (MatchedPhoton != 0) {
    printf("MatchedPhoton DeltaR Id MotherId LepPt PhoPt: %8.5f %5i %5i %12.2f %12.2f\n",
        MatchedPhoton->DeltaR(Zll[2]),
        MatchedPhoton->GetId(),
        MatchedPhoton->GetMotherId(),
        Zll[2].Perp(),
        MatchedPhoton->Perp());

    Hist.FillTH1D("Pt_PhoMatch_"+fProcName, 100, 0, 200, Zll[2].Perp());

    // Plot for tagged conversion and non-tagged where we've matched a photon genp
    if (Zll[2].GetIsConvertedPhoton()) {
      Hist.FillTH1D("ZmmMass_PhoMatch_TaggedConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
    } else {
      Hist.FillTH1D("ZmmMass_PhoMatch_NoTagConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
    }
  } else {
    if (Zll[2].GetIsConvertedPhoton()) {
      Hist.FillTH1D("ZmmMass_NotPhoMatch_TaggedConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
    } else {
      Hist.FillTH1D("ZmmMass_NotPhoMatch_NoTagConv"+fProcName, 100, 0, 200, (Zll[0] + Zll[1]).M());
    }
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



void SimpleAna::Plot6JetEvents ()
{
  static TAnaHist Hist(fOutFile, "Plot6JetEvents");

  if (Jets.size() < 6) {
    return;
  }

  if (Jets[0].Perp() < 50) return;
  if (Jets[1].Perp() < 40) return;
  if (Jets[2].Perp() < 30) return;

  for (size_t i = 0; i < Jets.size()-2; ++i) {
    for (size_t j = i+1; j < Jets.size()-1; ++j) {
      for (size_t k = j+1; k < Jets.size(); ++k) {
        float const SumEt = Jets[i].Perp() + Jets[j].Perp() + Jets[k].Perp();
        float const Mass = (Jets[i] + Jets[j] + Jets[k]).M();
        Hist.FillTH2D("SumEtVsTripletMass", 1000, 0, 900, 1000, 0, 1000, SumEt, Mass);
      }
    }
  }



  return;
}



void SimpleAna::PlotNLepSumEtJetsMet ()
{
  static TAnaHist Hist(fOutFile, "PlotNLetSumEtJetsMet");

  float SumEtJets = 0.0;
  for (size_t i = 0; i != Jets.size(); ++i) {
    SumEtJets += Jets[i].Perp();
  }

  Hist.FillTH3D("NLepSumEtJetsMet", 5, 0, 5, 20, 0, 500, 20, 0, 200, Leptons.size(), SumEtJets, GetMet());
  Hist.FillTH3D("NLepNJetsMet", 5, 0, 5, 8, 0, 8, 20, 0, 200, Leptons.size(), Jets.size(), GetMet());

  return;
}




TGenP* SimpleAna::FindClosestGenP (TLepton& Lep, int const Type)
{
  float MinDeltaR = 999999;
  int MinIndex = -1;
  for (size_t i = 0; i != Lep.GetGenPVector()->size(); ++i) {
    TGenP* ThisGenP = Lep.GetGenP(i);
    if (Type != 0 && !ThisGenP->IsId(Type)) {
      continue;
    }

    if (TMath::Abs(Lep.DeltaR(*ThisGenP)) < MinDeltaR) {
      MinDeltaR = Lep.DeltaR(*ThisGenP);
      MinIndex = (int) i;
    }
    
  }

  if (MinIndex < 0) {
    return (TGenP*) 0x0;
  }

  return Lep.GetGenP(MinIndex);

}





std::vector<TLepton> SimpleAna::ClosestZMatch(std::vector<TLepton>& InLeptons, bool const RequireOS, bool const RequireSF, bool const AddOtherLeptonsAtEnd)
{
  // This function will return the closest Z match from the leptons given to it.
  // This function will return a vector of leptons
  // There are a few options here.  You can require opposite sign leptons, same flavor,
  // and optionally you can tack the other leptons on the end.  The latch will always be
  // indices 0 and 1.  If OS and SF are specified but no match can be found satisfying 
  // all criterion, then an empty vector is returned

  // The vector of leptons we will return
  std::vector<TLepton> Zll;

  // This isn't even worth doing if there isn't 2 leptons
  if (InLeptons.size() < 2) {
    return Zll;
  }

  // 2 indices of best z match and best mass diff (wii find in alg below)
  size_t BestZll[2];
  float BestDiffMass = 999999;

  // Loop over all lepton pairs i-j
  for (size_t i = 0; i < InLeptons.size()-1; ++i) {
    for (size_t j = i+1; j < InLeptons.size(); ++j) {

      // Check for Opposite sign requirement
      if (RequireOS && (InLeptons[i].GetCharge() + InLeptons[j].GetCharge() != 0)) {
        continue;
      }

      // Check for same flavor requirement
      if (RequireSF && (InLeptons[i].GetFlavor() != InLeptons[j].GetFlavor())) {
        continue;
      }

      // Get the mass difference and if it's the best save indices and diff
      if ( TMath::Abs((InLeptons[i] + InLeptons[j]).M() - 91.0) < BestDiffMass ) {
        BestDiffMass = TMath::Abs((InLeptons[i] + InLeptons[j]).M() - 91.0);
        BestZll[0] = i;
        BestZll[1] = j;
      }
    }
  }

  // This means no match was found so return empty vector
  if (BestDiffMass == 999999) {
    return Zll;
  }

  // A match was found so add the two matching leptons to the vector
  Zll.push_back( InLeptons[ BestZll[0] ] );
  Zll.push_back( InLeptons[ BestZll[1] ] );

  // If you want the other leptons tacked on the end add them here
  if (AddOtherLeptonsAtEnd) {
    for (size_t i = 0; i < InLeptons.size(); ++i) { 
      if (i != BestZll[0] && i != BestZll[1]) {
        Zll.push_back(InLeptons[i]);
      }
    }
  }

  // Return vector of leptons
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
      bool Keep = true;
      if (!lep->PassesSelection(TLepton::kElectronSel_RobustTight)) Keep = false;
      if ( !(TMath::Abs(lep->Eta()) < 2.4) ) Keep = false;
      if ( lep->GetCalIso() / lep->Perp() > 0.1) Keep = false;
      if (Keep) {
        NewLeptons.push_back(*lep);
      }
    } else if (lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      bool Keep = true;
      if (!lep->PassesSelection(TLepton::kMuonSel_GlobalMuonPromptTight)) Keep = false;
      if ( !(TMath::Abs(lep->Eta()) < 2.1) ) Keep = false;
      if ( lep->GetCalIso() / lep->Perp() > 0.1) Keep = false;
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
    if ( TMath::Abs(photon->Eta()) > 2.5) Keep = false;
    if (photon->Perp() < 15.0) Keep = false;
    if (photon->GetCalIso() / photon->Perp() > 0.1) Keep = false;
    if (Keep) {
      NewPhotons.push_back(*photon);
    }
  }

  Photons = NewPhotons;

  return;
}



void SimpleAna::SelectionJet ()
{
  std::vector<TJet> NewJets;

  for (std::vector<TJet>::iterator jet = Jets.begin(); jet != Jets.end(); ++jet) {
    bool Keep = true;
    if (jet->Perp() < 20) {
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
