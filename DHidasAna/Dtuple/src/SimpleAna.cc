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

  fFakeDtuple = new TDtuple( new TFile("Fakes.root", "recreate") );
  fFakeDtuple->GetTree()->Delete();
  fFakeDtuple->SetTree( Chain->CloneTree(0) );

}



SimpleAna::~SimpleAna ()
{
  fOutFile->Write();
  fOutFile->Close();
  delete fOutFile;

  fFakeDtuple->Write();
  delete fFakeDtuple;
}


void SimpleAna::BeginJob ()
{
  // Do things you want to do here at the end..
  std::cout << "SimpleAna::BeginJob() called" << std::endl;
  return;
}



void SimpleAna::SetFakeRateFile (TString const FileName)
{
  fFakeRate = new TFakeRate(FileName);
  return;
}


void SimpleAna::Analyze (long unsigned int const ientry)
{
  //ElectronJetTest();
  PlotFakes();
  PlotElectronId();

  Selection();
  ObjectCleaning();
  AddFakesToEvent();

  PlotEventQuantities();
  PlotLeptons();
  PlotPhotons();
  PlotJets();

  PlotDileptonMass();
  PlotTriLeptons();
  PlotlGamma();
  PlotllGamma();
  PlotZllE();
  Plot6JetEvents();
  PlotNLepSumEtJetsMet();

  return;
}



void SimpleAna::PlotEventQuantities ()
{
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

    // Plot conversions..
    if (Lep->GetIsConvertedPhoton() == 1) {
      Hist.FillTH1D("LeptonConvPhi", 50, 0, 2*TMath::Pi(), Lep->GetConvPhi());
      Hist.FillTH1D("LeptonConvR", 50, 0, 20, Lep->GetConvR());
      Hist.FillTH2D("LeptonConvPhiR", 1000, 0, 2*TMath::Pi(), 1000, 0, 20, Lep->GetConvPhi(), Lep->GetConvR());
    }
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



void SimpleAna::PlotTriLeptons ()
{
  static TAnaHist Hist(fOutFile, "PlotTriLeptons");

  // Require three leptons
  if (Leptons.size() != 3) {
    return;
  }

  // Get flavor string
  TString const Flavors = GetLeptonFlavorsString(Leptons);

  // Get total charge string
  char ChargeStr[2];
  sprintf(ChargeStr, "%1i", abs(Leptons[0].GetCharge() + Leptons[1].GetCharge() + Leptons[2].GetCharge()));

  // Get the closest Z match
  std::vector<TLepton> Zll = ClosestZMatch(Leptons, true, true, true);
  if (Zll.size() != 3) {
    return;
  }

  // Loop over the leptons, look for electron, and see what it's from.
  for (size_t i = 0; i != Zll.size(); ++i) {
    if (Zll[i].IsFlavor(TLepton::kLeptonFlavor_Electron)) {

      // if it roughly matches to a generator electron then skip it
      if (Zll[i].GenPMatchesTo(11, 0.2, 0.15, false)) {
        continue;
      }

      // Get the closest GenP, or 0 for no match
      TGenP* MyGenP = Zll[i].GetClosestGenP();
      if (MyGenP) {
        printf("Electron Pt: %7.2f\n", Zll[i].Perp());
        Zll[i].PrintGenP();
        Hist.FillTH1D("EMother_"+Flavors, 10000, 0, 10000, TMath::Abs(MyGenP->GetId()) );
        fPlotTriLeptons_ElectronGenPMap[TMath::Abs(MyGenP->GetId())]++;
        fPlotTriLeptons_ElectronGenPMotherMap[ std::make_pair<int, int>(TMath::Abs(MyGenP->GetId()), TMath::Abs(MyGenP->GetMotherId()))]++;
      } else {
        Hist.FillTH1D("EMother_"+Flavors, 10000, 0, 10000, 0 );
        fPlotTriLeptons_ElectronGenPMap[0]++;
        fPlotTriLeptons_ElectronGenPMotherMap[ std::make_pair<int, int>(0, 0)]++;
      }
    }
  }

  TString const FlavorOrder = Zll[0].GetFlavorString() + Zll[1].GetFlavorString() + Zll[2].GetFlavorString();
  TString const FlavorZll = Zll[0].GetFlavorString() + Zll[1].GetFlavorString();

  Hist.FillTH1D("ZMass_"+FlavorZll+"_"+Flavors, 100, 0, 200, (Zll[0] + Zll[1]).M());
  Hist.FillTH1D("ZMass_"+FlavorZll+"_"+Flavors+"_Q"+ChargeStr, 100, 0, 200, (Zll[0] + Zll[1]).M());

  return;
}




void SimpleAna::PlotlGamma ()
{
  // For histograms
  static TAnaHist Hist(fOutFile, "PlotlGamma");

  if (Leptons.size() != 2) {
    return;
  }

  TGenP* MatchedPhotonA = FindClosestGenP(Leptons[0], 22);
  TGenP* MatchedPhotonB = FindClosestGenP(Leptons[1], 22);

  int i = 1;
  if (!MatchedPhotonA) {
    MatchedPhotonA = MatchedPhotonB;
    i = 0;
  }
  if (!MatchedPhotonA) {
    return;
  }

  if (MatchedPhotonA && MatchedPhotonB) {
    if (   TMath::Abs(MatchedPhotonA->Perp() - Leptons[0].Perp())/Leptons[0].Perp()
         > TMath::Abs(MatchedPhotonB->Perp() - Leptons[1].Perp())/Leptons[1].Perp() ) {
      MatchedPhotonA = MatchedPhotonB;
      i = 0;
    }
    //if (MatchedPhotonA->DeltaR(Leptons[0]) > MatchedPhotonB->DeltaR(Leptons[1])) {
    //  MatchedPhotonA = MatchedPhotonB;
    //  i = 0;
    //}
  }

  if (false) {
    printf("Converted photon RecoPt: %9.3f  GenPt: %9.3f  DeltaR: %7.5f\n",
      Leptons[1-i].Perp(), MatchedPhotonA->Perp(), MatchedPhotonA->DeltaR(Leptons[1-i]));
  }

  if (!Leptons[1-i].IsFlavor(TLepton::kLeptonFlavor_Electron)) {
    return;
  }

  if (!Leptons[i].IsFlavor(TLepton::kLeptonFlavor_Electron)) {
    return;
  }

  Hist.FillTH1D("lg_gIsConversion", "Is Conversion", "IsConversion bit", "", 2, 0, 2, Leptons[1-i].GetIsConvertedPhoton());
  Hist.FillTH1D("lg_lIsConversion", "Is Conversion", "IsConversion bit", "", 2, 0, 2, Leptons[i].GetIsConvertedPhoton());
  Hist.FillTH1D("lg_gDZ", "Conversion DZ", "DZ", "", 50, -10, 10,  Leptons[1-i].Getdz());
  Hist.FillTH1D("lg_gDXY", "Conversion DXY", "DXY", "", 50, -0.5, 0.5,  Leptons[1-i].Getdxy());
  Hist.FillTH1D("lg_lDZ", "Electron DZ", "DZ", "", 50, -10, 10,  Leptons[i].Getdz());
  Hist.FillTH1D("lg_lDXY", "Electron DXY", "DXY", "", 50, -0.5, 0.5,  Leptons[i].Getdxy());
  Hist.FillTH1D("lg_lgDiffDZ", "Electron-Conversion DZ Diff", "DZ(|e-Conv|)", "", 50, 0, 10,  TMath::Abs(Leptons[i].Getdz() - Leptons[1-i].Getdz()) );

  return;
}



void SimpleAna::PlotllGamma ()
{
  // For histograms
  static TAnaHist Hist(fOutFile, "PlotllGamma");

  if (Leptons.size() != 3) {
    return;
  }


  // Get the closest Z-match, require they are opposiute sign and same flavor
  std::vector<TLepton> Zll = ClosestZMatch(Leptons, true, true, true);
  if (Zll.size() != 3) {
    return;
  }

  // Make sure the third one is an electron
  if (!Zll[2].IsFlavor(TLepton::kLeptonFlavor_Electron)) {
    return;
  }

  // Match the electron to a genp photon?
  TGenP* MatchedPhoton = FindClosestGenP(Zll[2], 22);
  if (!MatchedPhoton) {
    return;
  }


  TString const Flavors = GetLeptonFlavorsString(Zll);
  Hist.FillTH1D(Flavors+"_gIsConversion", "Is Conversion", "IsConversion bit", "", 2, 0, 2, Zll[2].GetIsConvertedPhoton());
  Hist.FillTH1D(Flavors+"_lIsConversion", "Is Conversion", "IsConversion bit", "", 2, 0, 2, Zll[0].GetIsConvertedPhoton(), 0.5);
  Hist.FillTH1D(Flavors+"_lIsConversion", "Is Conversion", "IsConversion bit", "", 2, 0, 2, Zll[1].GetIsConvertedPhoton(), 0.5);

  Hist.FillTH1D(Flavors+"_gDXY", "DXY", "DXY", "", 200, -0.3, 0.3, Zll[2].Getdxy());
  Hist.FillTH1D(Flavors+"_lDXY", "DXY", "DXY", "", 200, -0.3, 0.3, Zll[0].Getdxy(), 0.5);
  Hist.FillTH1D(Flavors+"_lDXY", "DXY", "DXY", "", 200, -0.3, 0.3, Zll[1].Getdxy(), 0.5);

  if (Zll[2].GetIsConvertedPhoton()) {
    Hist.FillTH1D(Flavors+"_gDXYTagged", "DXYTagged", "DXY", "", 200, -0.3, 0.3, Zll[2].Getdxy());
  } else {
    Hist.FillTH1D(Flavors+"_gDXYNotTagged", "DXYTagged", "DXY", "", 200, -0.3, 0.3, Zll[2].Getdxy());
  }

  Hist.FillTH1D(Flavors+"_gDZ", "DZ", "DZ", "", 50, -10, 10, Zll[2].Getdz());
  Hist.FillTH1D(Flavors+"_lDZ", "DZ", "DZ", "", 50, -10, 10, Zll[0].Getdz(), 0.5);
  Hist.FillTH1D(Flavors+"_lDZ", "DZ", "DZ", "", 50, -10, 10, Zll[1].Getdz(), 0.5);


  return;
}




void SimpleAna::PlotZllE ()
{
  // For histograms
  static TAnaHist Hist(fOutFile, "PlotZllE");

  // THis is only for trileptons
  if (Leptons.size() != 3) {
    return;
  }

  // Get the closest Zmatch and add the 3dr lepton on at the end as [2]
  std::vector<TLepton> Zll = ClosestZMatch(Leptons, true, true, true);
  if (Zll.size() != 3) {
    return;
  }

  // Check that the muons were the match and the electron is in [2]
  if ( !( Zll[0].IsFlavor(TLepton::kLeptonFlavor_Muon) && Zll[2].IsFlavor(TLepton::kLeptonFlavor_Electron) ) ) {
    return;
  }

  // Plot some basic quantities for this trilepton candidate
  Hist.FillTH1D("mme_eIsConversion", 2, 0, 2, Zll[2].GetIsConvertedPhoton());
  Hist.FillTH1D("mme_mmMass", 100, 0, 200, (Zll[0] + Zll[1]).M());
  Hist.FillTH1D("mme_DeltaRMin_em", 25, 0, 4, std::min( Zll[2].DeltaR(Zll[0]), Zll[2].DeltaR(Zll[1]) ) );
  Hist.FillTH1D("mme_mmAvgDZ", 50, -10, 10,  (Zll[0].Getdz() + Zll[1].Getdz())/2.0);
  Hist.FillTH1D("mme_eDZ", "", "mm Mass", "Events", 50, -10, 10,  Zll[2].Getdz());
  Hist.FillTH1D("mme_eDXY", "", "Electron DXY", "Events", 50, -0.5, 0.5,  Zll[2].Getdxy());
  Hist.FillTH1D("mme_mDXY", "", "Muon DXY", "Events", 50, -0.5, 0.5,  Zll[0].Getdxy(), 0.5);
  Hist.FillTH1D("mme_mDXY", "", "Muon DXY", "Events", 50, -0.5, 0.5,  Zll[1].Getdxy(), 0.5);

  // Plot for tagged conversion and non-tagged
  if (Zll[2].GetIsConvertedPhoton()) {
    Hist.FillTH1D("mme_mmMass_eTaggedAsConv", "Electron tagged as conversion", "mm Mass", "Events", 100, 0, 200, (Zll[0] + Zll[1]).M());
    Hist.FillTH1D("mme_eTaggedAsConv_DZ", 50, -10, 10,  Zll[2].Getdz());
  } else {
    Hist.FillTH1D("mme_mmMass_eNotTaggedAsConv", "Electron not tagged as conversion", "mm Mass", "Events", 100, 0, 200, (Zll[0] + Zll[1]).M());
    Hist.FillTH1D("mme_eNotTaggedAsConv_DZ", 50, 10, 10,  Zll[2].Getdz());
  }

  // Match to the closest photon in DeltaR
  TGenP* MatchedPhoton = FindClosestGenP(Zll[2], 22);
  if (MatchedPhoton != 0) {
    
    if (false) printf("MatchedPhoton DeltaR Id MotherId LepPt PhoPt: %8.5f %5i %5i %12.2f %12.2f\n",
        MatchedPhoton->DeltaR(Zll[2]),
        MatchedPhoton->GetId(),
        MatchedPhoton->GetMotherId(),
        Zll[2].Perp(),
        MatchedPhoton->Perp());

    Hist.FillTH1D("mme_ePt_PhoMatched", 100, 0, 200, Zll[2].Perp());
    Hist.FillTH1D("mme_eGammaMatch_DeltaR", 100, 0, 0.4, Zll[2].DeltaR(*MatchedPhoton));
    Hist.FillTH1D("mme_eGammaMatch_MotherId", 40, 0, 40, abs(MatchedPhoton->GetMotherId()) );

    // Plot for tagged conversion and non-tagged where we've matched a photon genp
    if (Zll[2].GetIsConvertedPhoton()) {
      Hist.FillTH1D("mme_mmMass_ePhoMatch_TaggedAsConv", 100, 0, 200, (Zll[0] + Zll[1]).M());
    } else {
      Hist.FillTH1D("mme_mmMass_ePhoMatch_NotTaggedAsConv", 100, 0, 200, (Zll[0] + Zll[1]).M());
    }
  } else {
    if (Zll[2].GetIsConvertedPhoton()) {
      Hist.FillTH1D("mme_mmMass_eNotPhoMatch_TaggedAsConv", 100, 0, 200, (Zll[0] + Zll[1]).M());
    } else {
      Hist.FillTH1D("mme_mmMass_eNotPhoMatch_NotTaggedAsConv", 100, 0, 200, (Zll[0] + Zll[1]).M());
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

  Hist.FillTH1D("DileptonMass"+Charges+"_"+Flavors+"", 100, 0, 200, (Leptons[0]+Leptons[1]).M());
  Hist.FillTH1D("le_leDiffDZ_"+Charges+"_"+Flavors+"", "Conversion DZ Diff for dilepton", "DZ(|l1-l2|)", "", 50, 0, 10,  TMath::Abs(Leptons[0].Getdz() - Leptons[1].Getdz()) );

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
  if (!Lep.PassesSelection(TLepton::kMuonSel_GlobalMuonPromptTight)) Pass = false;
  if ( (TMath::Abs(Lep.Eta()) > 2.4) ) Pass = false;
  if ( Lep.GetCalIso() / Lep.Perp() > 0.1) Pass = false;
  if ( Lep.Perp() < 15.0) Pass = false;
  //if ( TMath::Abs(Lep.Getdxy()) > 0.02) Pass = false;
  return Pass;
}


bool SimpleAna::PassSelectionMuon (TLepton& Lep)
{
  bool Pass = true;
  if ( !Lep.PassesSelection(TLepton::kMuonSel_GlobalMuonPromptTight)) Pass = false;
  if ( (TMath::Abs(Lep.Eta()) > 2.1) ) Pass = false;
  if ( Lep.GetCalIso() / Lep.Perp() > 0.1) Pass = false;
  return Pass;
}


void SimpleAna::PlotElectronId ()
{
  static TAnaHist Hist(fOutFile, "PlotElectronId");

  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {
    if (Lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      Hist.FillTH1D("ElectronPt", 50, 0, 200, Lep->Perp());
    }
  }

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




void SimpleAna::GetElectronFromJet (TJet const& Jet, TLepton& Lepton)
{
  return;
}




void SimpleAna::AddFakesToEvent (int const NFakesToAdd)
{
  fFakeDtuple->Fill();
  
  return;
}




void SimpleAna::ElectronJetTest ()
{
  // A test looking at jet and electron overlap

  static TAnaHist Hist(fOutFile, "ElectronJetTest");

  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {
    if (Lep->IsFlavor(TLepton::kLeptonFlavor_Electron) && Lep->Perp() > 30) {
      bool FoundJetMatch = false;
      for (std::vector<TJet>::iterator Jet = Jets.begin(); Jet != Jets.end(); ++Jet) {
        if (Lep->DeltaR(*Jet) < 0.2) {
          FoundJetMatch = true;
          Hist.FillTH2D("ElectronPtVsJetPt", 1000, 0, 200, 1000, 0, 200, Lep->Perp(), Jet->Perp());
        }
      }
      if (!FoundJetMatch) {
        std::cout << "Did not find a jet match for this electron Pt: " << Lep->Perp() << "  " << Lep->Eta() << std::endl;
      }
    }
  }

  return;
}




void SimpleAna::PlotFakes ()
{
  // A test looking at jet and electron overlap

  static TAnaHist Hist(fOutFile, "PlotFakes");
  static int const NBinsPt = 2;
  static float const PtMin = 0;
  static float const PtMax = 100;
  //static float const PtBins[NBinsPt] = {1,2,3,4,5,6};
  static int const NBinsEta = 4;
  static float const EtaMin = -3;
  static float const EtaMax = 3;
  //static float const EtaBins[NBinsPt] = {1,2,3,4,5,6};

  for (std::vector<TJet>::iterator Jet = Jets.begin(); Jet != Jets.end(); ++Jet) {
    if (TMath::Abs(Jet->Eta()) < 2.1 && Jet->GetHadF() / Jet->GetEmF() > 0.05) {
      Hist.FillTH1D("EleFakeDenomJetPt", "EleFakeDenomJetPt", "P_{T}", "", NBinsPt, PtMin, PtMax, Jet->Perp());
      Hist.FillTH2D("EleFakeDenomJetEtaPt", "EleFakeDenomJetEtaPt", "#eta", "P_{T}", NBinsEta, EtaMin, EtaMax, NBinsPt, PtMin, PtMax, Jet->Eta(), Jet->Perp());
    }
  }

  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {
    if (Lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      if (PassSelectionElectron(*Lep)) {
        Hist.FillTH1D("EleFakeNumeratorPt", "EleFakeNumeratorPt", "P_{T}", "", NBinsPt, PtMin, PtMax, Lep->Perp());
        Hist.FillTH2D("EleFakeNumeratorEtaPt", "EleFakeNumeratorEtaPt", "#eta", "P_{T}", NBinsEta, EtaMin, EtaMax, NBinsPt, PtMin, PtMax, Lep->Eta(), Lep->Perp());
      } else if (Lep->GetHCalOverECal() > 0.05) {
        Hist.FillTH1D("EleFakeDenomElePt", "EleFakeDenomElePt", "P_{T}", "", NBinsPt, PtMin, PtMax, Lep->Perp());
        Hist.FillTH2D("EleFakeDenomEleEtaPt", "EleFakeDenomEleEtaPt", "#eta", "P_{T}", NBinsEta, EtaMin, EtaMax, NBinsPt, PtMin, PtMax, Lep->Eta(), Lep->Perp());
      }
    }
    else if (Lep->IsFlavor(TLepton::kLeptonFlavor_Muon)) {
      if (PassSelectionMuon(*Lep)) {
        Hist.FillTH1D("MuonFakeNumeratorPt", "MuonFakeNumeratorPt", "P_{T}", "", NBinsPt, PtMin, PtMax, Lep->Perp());
        Hist.FillTH2D("MuonFakeNumeratorEtaPt", "MuonFakeNumeratorEtaPt", "#eta", "P_{T}", NBinsEta, EtaMin, EtaMax, NBinsPt, PtMin, PtMax, Lep->Eta(), Lep->Perp());
      } else if (false) {
        Hist.FillTH1D("MuonFakeDenomMuonPt", "MuonFakeDenomMuonPt", "P_{T}", "", NBinsEta, PtMin, PtMax, Lep->Perp());
        Hist.FillTH2D("MuonFakeDenomMuonEtaPt", "MuonFakeDenomMuonEtaPt", "#eta", "P_{T}", NBinsEta, EtaMin, EtaMax, NBinsPt, PtMin, PtMax, Lep->Eta(), Lep->Perp());
      }
    }
  }

  return;
}




void SimpleAna::EndJob ()
{
  // Do things you want to do here at the end..
  std::cout << "SimpleAna::EndJob() called" << std::endl;

  std::cout << "Summary:" <<std::endl;
  TDUtility::PrintMapIntInt(fPlotTriLeptons_ElectronGenPMap, "ElectronGenP");
  std::cout << "Detailed:" <<std::endl;
  TDUtility::PrintMapIntPairInt(fPlotTriLeptons_ElectronGenPMotherMap, "ElectronGenP");
  return;
}
