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


Dtuple::~Dtuple ()
{
}



void Dtuple::SetBranches (TTree* T)
{
  TString MaxLep = ""; MaxLep += kMaxLeptons;
  TString MaxJet = ""; MaxJet += kMaxJets;
  TString MaxCom = ""; MaxCom += kMaxCombinations;

  T->Branch("Run", &fEvt.Run);
  T->Branch("LumiSection",  &fEvt.LumiSection);
  T->Branch("Event",  &fEvt.Event);

  T->Branch("LeptonPx", &fEvt.LeptonPx, "LeptonPx["+MaxLep+"]/F");
  T->Branch("LeptonPy", &fEvt.LeptonPy, "LeptonPy["+MaxLep+"]/F");
  T->Branch("LeptonPz", &fEvt.LeptonPz, "LeptonPz["+MaxLep+"]/F");
  T->Branch("LeptonPt", &fEvt.LeptonPt, "LeptonPt["+MaxLep+"]/F");
  T->Branch("LeptonType", &fEvt.LeptonType, "LeptonType["+MaxLep+"]/I");
  T->Branch("NLeptons", &fEvt.NLeptons);

  T->Branch("JetPx", &fEvt.JetPx, "JetPx["+MaxJet+"]/F");
  T->Branch("JetPy", &fEvt.JetPy, "JetPy["+MaxJet+"]/F");
  T->Branch("JetPz", &fEvt.JetPz, "JetPz["+MaxJet+"]/F");
  T->Branch("JetPt", &fEvt.JetPt, "JetPt["+MaxJet+"]/F");
  T->Branch("NJets", &fEvt.NJets);

  T->Branch("SumPtJets", &fEvt.SumPtJets);
  T->Branch("TriJetSumPt",  &fEvt.TriJetSumPt,  "TriJetSumPt["+MaxCom+"]/F");
  T->Branch("TriJetMasses", &fEvt.TriJetMasses, "TriJetMasses["+MaxCom+"]/F");

  T->Branch("MET", &fEvt.MET);

  return;
}


void Dtuple::ClearDtuple ()
{
  fEvt.Run = -1;
  fEvt.Event = -1;
  fEvt.LumiSection = -1;

  fEvt.NLeptons = 0;
  for (int i = 0; i != kMaxLeptons; ++i) {
    fEvt.LeptonPx[i] = 0;
    fEvt.LeptonPy[i] = 0;
    fEvt.LeptonPz[i] = 0;
    fEvt.LeptonPt[i] = 0;
    fEvt.LeptonType[i] = -1;
  }

  fEvt.NJets = 0;
  for (int i = 0; i != kMaxJets; ++i) {
    fEvt.JetPx[i] = 0;
    fEvt.JetPy[i] = 0;
    fEvt.JetPz[i] = 0;
    fEvt.JetPt[i] = 0;
  }

  fEvt.SumPtJets = 0;
  for (int i = 0; i != Dtuple::kMaxCombinations; ++i) {
    fEvt.TriJetSumPt[i] = 0;
    fEvt.TriJetMasses[i] = 0;
  }

  fEvt.MET = 0;

  return;
}
