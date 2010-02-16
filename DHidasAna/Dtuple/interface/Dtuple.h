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

    int  GetEntry (unsigned long const);
    int  Fill ();


  protected:
    TFile* fDtupleFile;
    TTree* fDtupleTree;


  public:
    static int const NMaxLeptons = 4;
    static int const NMaxPhotons = 4;
    static int const NMaxJets = 10;

    enum LeptonFlavor {
      kLeptonFlavor_Electron,
      kLeptonFlavor_Muon,
      kLeptonFlavor_Tau
    };

    enum ElectronSel {
      kElectronSel_RobustHighEnergy = 0,
      kElectronSel_RobustLoose = 1,
      kElectronSel_RobustTight = 2,
      kElectronSel_Loose = 3,
      kElectronSel_Tight = 4
    };

    enum MuonSel {
      kMuonSel_TrackerMuonArbitrated = 0,
      kMuonSel_AllArbitrated = 1,
      kMuonSel_GlobalMuonPromptTight = 2,
      kMuonSel_TMLastStationLoose = 3,
      kMuonSel_TMLastStationTight = 4,
      kMuonSel_TM2DCompatibilityLoose = 5,
      kMuonSel_TM2DCompatibilityTight = 6,
      kMuonSel_TMOneStationLoose = 7,
      kMuonSel_TMOneStationTight = 8,
      kMuonSel_TMLastStationOptimizedLowPtLoose = 9,
      kMuonSel_TMLastStationOptimizedLowPtTight = 10
    };

    enum ElectronDet {
      kElectronDet_EE = 0,
      kElectronDet_EB = 1,
      kElectronDet_EBEEGap = 2,
      kElectronDet_EBEtaGap = 3,
      kElectronDet_EBGap = 4,
      kElectronDet_EBPhiGap = 5,
      kElectronDet_EEDeeGap = 6,
      kElectronDet_EEGap = 7,
      kElectronDet_EERingGap = 8
    };

    enum MuonDet {
      kMuonDet_Global = 0,
      kMuonDet_Tracker = 1,
      kMuonDet_StandAlone = 2,
      kMuonDet_Calo = 3
    };

    enum ElectronClass {
      kElectronClass_Unknown = -1,
      kElectronClass_Golden,
      kElectronClass_BigBrem,
      kElectronClass_Narrow,
      kElectronClass_Showering,
      kElectronClass_Gap
    };

  public:
    struct Event_Struct {
      int   Run;
      int   Event;
      float EventWeight;
      float TriggerEff;
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
      float Lepton_dxy[NMaxLeptons];
      float Lepton_dz[NMaxLeptons];
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
      float Lepton_fBrem[NMaxLeptons];
      int   Lepton_IsConvertedPhoton[NMaxLeptons];
      int   Lepton_PassSelection[NMaxLeptons];
      int   Lepton_Detector[NMaxLeptons];
      int   Lepton_Classification[NMaxLeptons];
      float Lepton_SigmaIEtaIEta[NMaxLeptons];
      float Lepton_DeltaEtaIn[NMaxLeptons];
      float Lepton_DeltaPhiIn[NMaxLeptons];
      float Lepton_E2x5overE5x5[NMaxLeptons];
      float Lepton_ConvDist[NMaxLeptons];
      float Lepton_ConvdCotTheta[NMaxLeptons];


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
    void DefaultValues ();
    void DefaultValues (Dtuple::Event_Struct&);


  public:
    static const TString fEvent_Format;
    Event_Struct fEvent;

};





#endif
