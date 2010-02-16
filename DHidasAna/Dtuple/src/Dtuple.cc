////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include "DHidasAna/Dtuple/interface/Dtuple.h"

#include <iostream>




// This is the format of the ntuple.  The order must be taken exactly from the struct.
// If you update this you need to update the struct.  Also, the size of the arrays
// is important.  If you update that here you must change the struct
const TString Dtuple::fEvent_Format = 
  "Run/I:"
  "Event/I:"
  "EventWeight/F:"
  "TriggerEff/F:"
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
  "Lepton_dxy[4]/F:"
  "Lepton_dz[4]/F:"
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
  "Lepton_fBrem[4]/F:"
  "Lepton_IsConvertedPhoton[4]/I:"
  "Lepton_PassSelection[4]/I:"
  "Lepton_Detector[4]/I:"
  "Lepton_Classification[4]/I:"
  "Lepton_SigmaIEtaIEta[4]/F:"
  "Lepton_DeltaEtaIn[4]/F:"
  "Lepton_DeltaPhiIn[4]/F:"
  "Lepton_E2x5overE5x5[4]/F:"

  "NPhotons/I:"
  "Photon_Px[4]/F:"
  "Photon_Py[4]/F:"
  "Photon_Pz[4]/F:"
  "Photon_Pt[4]/F:"
  "Photon_Eta[4]/F:"
  "Photon_Phi[4]/F:"
  "Photon_TrkIso[4]/F:"
  "Photon_CalIso[4]/F:"
  "Photon_HCalOverECal[4]/F:"

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
  // Get the i-th entry from the Dtuple tree.
  fDtupleTree->GetEntry(ientry);
  return fDtupleTree->GetEntry(ientry);
}


int Dtuple::Fill ()
{
  // Fill an entry!!

  //std::cout << "Begin fill" << std::endl;
  return fDtupleTree->Fill();
  //std::cout << "End fill" << std::endl;
}



void Dtuple::DefaultValues ()
{
  DefaultValues(fEvent);
  return;
}

void Dtuple::DefaultValues (Dtuple::Event_Struct& Ev)
{
  // Set default values for all entries in the dtuple

  //std::cout << "Begin default values" << std::endl;
  Ev.Run    = -999999;
  Ev.Event  = -999999;
  Ev.EventWeight  = -999999;
  Ev.TriggerEff   = -999999;
  Ev.MetMag = -999999;
  Ev.MetPhi = -999999;
  Ev.SumEt  = -999999;
  Ev.MetSig = -999999;

  Ev.NLeptons = 0;
  for (int i = 0; i != NMaxLeptons; ++i) {
    Ev.Lepton_Px[i] = -999999;
    Ev.Lepton_Py[i] = -999999;
    Ev.Lepton_Pz[i] = -999999;
    Ev.Lepton_Pt[i] = -999999;
    Ev.Lepton_TrkPt[i] = -999999;
    Ev.Lepton_Eta[i] = -999999;
    Ev.Lepton_Phi[i] = -999999;
    Ev.Lepton_dxy[i] = -999999;
    Ev.Lepton_dz[i] = -999999;
    Ev.Lepton_Z0[i] = -999999;
    Ev.Lepton_Charge[i] = -999999;
    Ev.Lepton_Flavor[i] = -999999;
    Ev.Lepton_TrkIso[i] = -999999;
    Ev.Lepton_CalIso[i] = -999999;
    Ev.Lepton_ECalIso[i] = -999999;
    Ev.Lepton_HCalIso[i] = -999999;
    Ev.Lepton_CalE[i] = -999999;
    Ev.Lepton_HCalOverECal[i] = -999999;
    Ev.Lepton_EoverPin[i] = -999999;
    Ev.Lepton_fBrem[i] = -999999;
    Ev.Lepton_IsConvertedPhoton[i] = -999999;
    Ev.Lepton_PassSelection[i] = 0x0;
    Ev.Lepton_Detector[i] = 0x0;
    Ev.Lepton_Classification[i] = -999999;
    Ev.Lepton_SigmaIEtaIEta[i] = -999999;
    Ev.Lepton_DeltaEtaIn[i] = -999999;
    Ev.Lepton_DeltaPhiIn[i] = -999999;
    Ev.Lepton_E2x5overE5x5[i] = -999999;

  }

  Ev.NPhotons = 0;
  for (int i = 0; i != NMaxPhotons; ++i) {
    Ev.Photon_Px[i] = -999999;
    Ev.Photon_Py[i] = -999999;
    Ev.Photon_Pz[i] = -999999;
    Ev.Photon_Pt[i] = -999999;
    Ev.Photon_Eta[i] = -999999;
    Ev.Photon_Phi[i] = -999999;
    Ev.Photon_TrkIso[i] = -999999;
    Ev.Photon_CalIso[i] = -999999;
    Ev.Photon_HCalOverECal[i] = -999999;
  }

  Ev.NJets = 0;
  for (int i = 0; i != NMaxJets; ++i) {
    Ev.Jet_Px[i] = -999999;
    Ev.Jet_Py[i] = -999999;
    Ev.Jet_Pz[i] = -999999;
    Ev.Jet_Pt[i] = -999999;
    Ev.Jet_Eta[i] = -999999;
    Ev.Jet_Phi[i] = -999999;
    Ev.Jet_EmF[i] = -999999;
    Ev.Jet_HadF[i] = -999999;
  }
  //std::cout << "End default values" << std::endl;

  return;
}
