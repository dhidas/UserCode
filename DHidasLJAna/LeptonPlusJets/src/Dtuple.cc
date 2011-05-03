////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Apr 26 03:47:29 EDT 2011
//
////////////////////////////////////////////////////////////////////

#include "DHidasLJAna/LeptonPlusJets/interface/Dtuple.h"

#include <iostream>

Dtuple::Dtuple ()
{
}


Dtuple::Dtuple (std::vector<TString> const& InFileNames)
{
  TChain* c = new TChain("d", "d");
  for (size_t i = 0; i != InFileNames.size(); ++i) {
    c->Add(InFileNames[i]);
  }

  fTree = (TTree*) c;
  SetBranchAddresses(fTree);
}


Dtuple::~Dtuple ()
{
}


int Dtuple::GetEntry (long long& ientry)
{
  int ret = fTree->GetEntry(ientry);
  if (ret <= 0) {
    return ret;
  }

  fEvt.Lep.resize(fEvt.NLeptons);
  for (int i = 0; i < fEvt.NLeptons; ++i) {
    fEvt.Lep[i].SetPxPyPzE(
        fEvt.LeptonPx->at(i),
        fEvt.LeptonPy->at(i),
        fEvt.LeptonPz->at(i),
        TMath::Sqrt( TMath::Power(fEvt.LeptonPt->at(i), 2) + TMath::Power(fEvt.LeptonPz->at(i), 2) ));
  }

  fEvt.Jet.resize(fEvt.NJets);
  for (int i = 0; i < fEvt.NJets; ++i) {
    fEvt.Jet[i].SetPxPyPzE(
        fEvt.JetPx->at(i),
        fEvt.JetPy->at(i),
        fEvt.JetPz->at(i),
        TMath::Sqrt( TMath::Power(fEvt.JetPt->at(i), 2) + TMath::Power(fEvt.JetPz->at(i), 2) ));
  }

  fEvt.MET.SetMagPhi(fEvt.METMag, fEvt.METPhi);

  return ret;
}


void Dtuple::SetBranchAddresses (TTree* T)
{
  fEvt.LeptonPx = (std::vector<float>*) 0x0;
  fEvt.LeptonPy = (std::vector<float>*) 0x0;
  fEvt.LeptonPz = (std::vector<float>*) 0x0;
  fEvt.LeptonPt = (std::vector<float>*) 0x0;
  fEvt.LeptonType = (std::vector<int>*) 0x0;

  fEvt.JetPx = (std::vector<float>*) 0x0;
  fEvt.JetPy = (std::vector<float>*) 0x0;
  fEvt.JetPz = (std::vector<float>*) 0x0;
  fEvt.JetPt = (std::vector<float>*) 0x0;

  fEvt.TriJetSumPt = (std::vector<float>*) 0x0;
  fEvt.TriJetMasses = (std::vector<float>*) 0x0;


  T->SetBranchAddress("Run", &fEvt.Run);
  T->SetBranchAddress("LumiSection",  &fEvt.LumiSection);
  T->SetBranchAddress("Event",  &fEvt.Event);

  T->SetBranchAddress("LeptonPx", &fEvt.LeptonPx);
  T->SetBranchAddress("LeptonPy", &fEvt.LeptonPy);
  T->SetBranchAddress("LeptonPz", &fEvt.LeptonPz);
  T->SetBranchAddress("LeptonPt", &fEvt.LeptonPt);
  T->SetBranchAddress("LeptonType", &fEvt.LeptonType);
  T->SetBranchAddress("NLeptons", &fEvt.NLeptons);

  T->SetBranchAddress("JetPx", &fEvt.JetPx);
  T->SetBranchAddress("JetPy", &fEvt.JetPy);
  T->SetBranchAddress("JetPz", &fEvt.JetPz);
  T->SetBranchAddress("JetPt", &fEvt.JetPt);
  T->SetBranchAddress("NJets", &fEvt.NJets);

  T->SetBranchAddress("SumPtJets", &fEvt.SumPtJets);
  T->SetBranchAddress("TriJetSumPt",  &fEvt.TriJetSumPt);
  T->SetBranchAddress("TriJetMasses", &fEvt.TriJetMasses);

  T->SetBranchAddress("METMag", &fEvt.METMag);
  T->SetBranchAddress("METPhi", &fEvt.METPhi);

  return;
}
void Dtuple::SetBranches (TTree* T)
{
  fEvt.LeptonPx       = new std::vector<float>();
  fEvt.LeptonPy       = new std::vector<float>();
  fEvt.LeptonPz       = new std::vector<float>();
  fEvt.LeptonPt       = new std::vector<float>();
  fEvt.LeptonType     = new std::vector<int>();
                                         
  fEvt.JetPx          = new std::vector<float>();
  fEvt.JetPy          = new std::vector<float>();
  fEvt.JetPz          = new std::vector<float>();
  fEvt.JetPt          = new std::vector<float>();
                                         
  fEvt.TriJetSumPt    = new std::vector<float>();
  fEvt.TriJetMasses   = new std::vector<float>();


  T->Branch("Run", &fEvt.Run);
  T->Branch("LumiSection",  &fEvt.LumiSection);
  T->Branch("Event",  &fEvt.Event);

  T->Branch("LeptonPx", fEvt.LeptonPx);
  T->Branch("LeptonPy", fEvt.LeptonPy);
  T->Branch("LeptonPz", fEvt.LeptonPz);
  T->Branch("LeptonPt", fEvt.LeptonPt);
  T->Branch("LeptonType", fEvt.LeptonType);
  T->Branch("NLeptons", &fEvt.NLeptons);

  T->Branch("JetPx", fEvt.JetPx);
  T->Branch("JetPy", fEvt.JetPy);
  T->Branch("JetPz", fEvt.JetPz);
  T->Branch("JetPt", fEvt.JetPt);
  T->Branch("NJets", &fEvt.NJets);

  T->Branch("SumPtJets", &fEvt.SumPtJets);

  T->Branch("TriJetSumPt",  fEvt.TriJetSumPt);
  T->Branch("TriJetMasses", fEvt.TriJetMasses);


  T->Branch("METMag", &fEvt.METMag);
  T->Branch("METPhi", &fEvt.METPhi);

  return;
}


void Dtuple::ClearDtuple ()
{
  fEvt.Run = -1;
  fEvt.Event = -1;
  fEvt.LumiSection = -1;

  fEvt.NLeptons = 0;
  fEvt.LeptonPx->clear();
  fEvt.LeptonPy->clear();
  fEvt.LeptonPz->clear();
  fEvt.LeptonPt->clear();
  fEvt.LeptonType->clear();

  fEvt.NJets = 0;
  fEvt.JetPx->clear();
  fEvt.JetPy->clear();
  fEvt.JetPz->clear();
  fEvt.JetPt->clear();

  fEvt.SumPtJets = 0;
  fEvt.TriJetSumPt->clear();
  fEvt.TriJetMasses->clear();

  fEvt.METMag = 0;
  fEvt.METPhi = 0;

  return;
}



