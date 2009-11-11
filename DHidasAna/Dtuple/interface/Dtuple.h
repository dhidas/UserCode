////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////

#ifndef GUARD_Dtuple_h
#define GUARD_Dtuple_h


#include "TString.h"
#include "TFile.h"
#include "TTree.h"

class Dtuple
{
  public:
    Dtuple () {};
    Dtuple (TString const&) {};
    virtual ~Dtuple () {};

    int GetEntry (unsigned long const);
    int Fill ();
    void DefaultValues ();


  protected:
    TFile* fDtupleFile;
    TTree* fDtupleTree;


  public:
    static int const NMaxLeptons = 4;
    static int const NMaxPhotons = 4;
    static int const NMaxJets = 10;

    enum LeptonFlavor {
      kElectron,
      kMuon,
      kTau
    };

  public:
    struct Event_Struct {
      int   Run;
      int   Event;
      float MetMag;
      float MetPhi;
      float SumEt;
      float MetSig;

      int   NLeptons;
      float Lepton_Px[NMaxLeptons];
      float Lepton_Py[NMaxLeptons];
      float Lepton_Pz[NMaxLeptons];
      float Lepton_Pt[NMaxLeptons];
      float Lepton_TrkPt[NMaxLeptons];
      float Lepton_Eta[NMaxLeptons];
      float Lepton_Phi[NMaxLeptons];
      float Lepton_D0[NMaxLeptons];
      float Lepton_Z0[NMaxLeptons];
      int   Lepton_Charge[NMaxLeptons];
      int   Lepton_Flavor[NMaxLeptons];
      float Lepton_TrkIso[NMaxLeptons];
      float Lepton_CalIso[NMaxLeptons];
      float Lepton_ECalIso[NMaxLeptons];
      float Lepton_HCalIso[NMaxLeptons];
      float Lepton_CalE[NMaxLeptons];
      float Lepton_HCalOverECal[NMaxLeptons];
      float Lepton_EoverPin[NMaxLeptons];

      int   NPhotons;
      float Photon_Px[NMaxPhotons];
      float Photon_Py[NMaxPhotons];
      float Photon_Pz[NMaxPhotons];
      float Photon_Pt[NMaxPhotons];
      float Photon_Eta[NMaxPhotons];
      float Photon_Phi[NMaxPhotons];
      float Photon_TrkIso[NMaxPhotons];
      float Photon_CalIso[NMaxPhotons];
      float Photon_HCalOverECal[NMaxLeptons];

      int   NJets;
      float Jet_Px[NMaxJets];
      float Jet_Py[NMaxJets];
      float Jet_Pz[NMaxJets];
      float Jet_Pt[NMaxJets];
      float Jet_Eta[NMaxJets];
      float Jet_Phi[NMaxJets];
      float Jet_EmF[NMaxJets];
      float Jet_HadF[NMaxJets];
    };


  public:
    static const TString fEvent_Format;
    Event_Struct fEvent;

};





#endif
