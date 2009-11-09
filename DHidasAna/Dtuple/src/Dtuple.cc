////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include "DHidasAna/Dtuple/interface/Dtuple.h"

#include <iostream>




const TString Dtuple::fEvent_Format = 
  "Run/I:"
  "Event/I:"
  "MetMag/F:"
  "MetPhi/F:"
  "SumEt/F:"
  "MetSig/F:"

  "NLeptons/I:"
  "Lepton_Px[4]/F:"
  "Lepton_Py[4]/F:"
  "Lepton_Pz[4]/F:"
  "Lepton_Pt[4]/F:"
  "Lepton_TrkPt[4]/F:"
  "Lepton_Eta[4]/F:"
  "Lepton_Phi[4]/F:"
  "Lepton_D0[4]/F:"
  "Lepton_Z0[4]/F:"
  "Lepton_Charge[4]/I:"
  "Lepton_Flavor[4]/I:"
  "Lepton_TrkIso[4]/F:"
  "Lepton_CalIso[4]/F:"
  "Lepton_ECalIso[4]/F:"
  "Lepton_HCalIso[4]/F:"
  "Lepton_CalE[4]/F:"
  "Lepton_HCalOverECal[4]/F:"
  "Lepton_EoverPin[4]/F:"

  "NPhotons/I:"
  "Photon_Px[4]/F:"
  "Photon_Py[4]/F:"
  "Photon_Pz[4]/F:"
  "Photon_Pt[4]/F:"
  "Photon_Eta[4]/F:"
  "Photon_Phi[4]/F:"
  "Photon_TrkIso[4]/F:"
  "Photon_CalIso[4]/F:"

  "NJets/I:"
  "Jet_Px[10]/F:"
  "Jet_Py[10]/F:"
  "Jet_Pz[10]/F:"
  "Jet_Pt[10]/F:"
  "Jet_Eta[10]/F:"
  "Jet_Phi[10]/F:"
  "Jet_EmF[10]/F:"
  "Jet_HadF[10]/F";




int Dtuple::GetEntry (unsigned long const ientry)
{
  return fDtupleTree->GetEntry(ientry);
}


int Dtuple::Fill ()
{
  //std::cout << "Begin fill" << std::endl;
  return fDtupleTree->Fill();
  //std::cout << "End fill" << std::endl;
}


void Dtuple::DefaultValues ()
{
  //std::cout << "Begin default values" << std::endl;
  fEvent.Run    = -999999;
  fEvent.Event  = -999999;
  fEvent.MetMag = -999999;
  fEvent.MetPhi = -999999;
  fEvent.SumEt  = -999999;
  fEvent.MetSig = -999999;

  fEvent.NLeptons = -1;
  for (int i = 0; i != NMaxLeptons; ++i) {
    fEvent.Lepton_Px[i] = -999999;
    fEvent.Lepton_Py[i] = -999999;
    fEvent.Lepton_Pz[i] = -999999;
    fEvent.Lepton_Pt[i] = -999999;
    fEvent.Lepton_TrkPt[i] = -999999;
    fEvent.Lepton_Eta[i] = -999999;
    fEvent.Lepton_Phi[i] = -999999;
    fEvent.Lepton_D0[i] = -999999;
    fEvent.Lepton_Z0[i] = -999999;
    fEvent.Lepton_Charge[i] = -999999;
    fEvent.Lepton_Flavor[i] = -999999;
    fEvent.Lepton_TrkIso[i] = -999999;
    fEvent.Lepton_CalIso[i] = -999999;
    fEvent.Lepton_ECalIso[i] = -999999;
    fEvent.Lepton_HCalIso[i] = -999999;
    fEvent.Lepton_CalE[i] = -999999;
    fEvent.Lepton_HCalOverECal[i] = -999999;
    fEvent.Lepton_EoverPin[i] = -999999;
  }

  fEvent.NPhotons = -1;
  for (int i = 0; i != NMaxPhotons; ++i) {
    fEvent.Photon_Px[i] = -999999;
    fEvent.Photon_Py[i] = -999999;
    fEvent.Photon_Pz[i] = -999999;
    fEvent.Photon_Pt[i] = -999999;
    fEvent.Photon_Eta[i] = -999999;
    fEvent.Photon_Phi[i] = -999999;
    fEvent.Photon_TrkIso[i] = -999999;
    fEvent.Photon_CalIso[i] = -999999;
  }

  fEvent.NJets = -1;
  for (int i = 0; i != NMaxJets; ++i) {
    fEvent.Jet_Px[i] = -999999;
    fEvent.Jet_Py[i] = -999999;
    fEvent.Jet_Pz[i] = -999999;
    fEvent.Jet_Pt[i] = -999999;
    fEvent.Jet_Eta[i] = -999999;
    fEvent.Jet_Phi[i] = -999999;
    fEvent.Jet_EmF[i] = -999999;
    fEvent.Jet_HadF[i] = -999999;
  }
  //std::cout << "End default values" << std::endl;

  return;
}
