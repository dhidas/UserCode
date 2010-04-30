#ifndef GUARD_TLepton_h
#define GUARD_TLepton_h

#include <iostream>
#include <string>
#include <vector>
#include <ostream>
#include <algorithm>

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "DHidasAna/Dtuple/interface/TGenP.h"





//
// Class: TLepton
//
// Purpose: This is a class to hold information about leptons.  It
//          does so by storing some information to private variables
//          and then has methods to access those variables as well as
//          methods which combine these basic variables into other quantaties
//          which one might be interested in.  The basic variables which should
//          be filled are listed below.  They are private so as to be nondestructive.
//          In other words, it's harder to mess them up...
//
class TLepton : public TLorentzVector
{
  // Constructors and destructor.
  // More of these will be added in time.
  public:
    TLepton ();
    ~TLepton ();


  // The private variables which need to be
  // set by the various methods
  private:
    float TrkPt;
    float dxy;
    float dz;
    float Z0;
    int   Charge;
    int   Flavor;
    float TrkIso;
    float CalIso;
    float ECalIso;
    float HCalIso;
    float CalE;
    float HCalOverECal;
    float EoverPin;
    float fBrem;
    int   IsConvertedPhoton;
    int   PassSelection;
    int   Detector;
    int   Classification;
    float SigmaIEtaIEta;
    float DeltaEtaIn;
    float DeltaPhiIn;
    float E2x5overE5x5;
    float ConvDist;
    float ConvdCotTheta;
    float ConvR;
    float ConvPhi;
    int   NValidTrackerHits;
    float TrackChi2;
    float TrackNDoF;
    float ECalIsoDep;
    float HCalIsoDep;
    bool  IsFake;

  public:
    //TGenP GenP;
    std::vector<TGenP> GenP;



  // Public methods for accessing the variables in
  // this class
  public:
    void SetTrkPt (float const);
    void Setdxy (float const);
    void Setdz (float const);
    void SetZ0 (float const);
    void SetCharge (int const);
    void SetFlavor (int const);
    void SetTrkIso (float const);
    void SetCalIso (float const);
    void SetECalIso (float const);
    void SetHCalIso (float const);
    void SetCalE (float const);
    void SetHCalOverECal (float const);
    void SetEoverPin (float const);
    void SetfBrem (float const);
    void SetIsConvertedPhoton (int const);
    void SetPassSelection (int const);
    void SetDetector (int const);
    void SetClassification (int const);
    void SetSigmaIEtaIEta (float const);
    void SetDeltaEtaIn (float const);
    void SetDeltaPhiIn (float const);
    void SetE2x5overE5x5 (float const);
    void SetConvDist (float const);
    void SetConvdCotTheta (float const);
    void SetConvR (float const);
    void SetConvPhi (float const);
    void SetNValidTrackerHits (int const);
    void SetTrackChi2 (float const);
    void SetTrackNDoF (float const);
    void SetECalIsoDep (float const);
    void SetHCalIsoDep (float const);
    void SetIsFake (bool const);


    float GetTrkPt ();
    float Getdxy ();
    float Getdz ();
    float GetZ0 ();
    int GetCharge ();
    int GetFlavor ();
    float GetTrkIso ();
    float GetCalIso ();
    float GetECalIso ();
    float GetHCalIso ();
    float GetCalE ();
    float GetHCalOverECal ();
    float GetEoverPin ();
    float GetfBrem ();
    int GetIsConvertedPhoton ();
    int GetPassSelection ();
    int GetDetector ();
    int GetClassification ();
    float GetSigmaIEtaIEta ();
    float GetDeltaEtaIn ();
    float GetDeltaPhiIn ();
    float GetE2x5overE5x5 ();
    float GetConvDist ();
    float GetConvdCotTheta ();
    float GetConvR ();
    float GetConvPhi ();
    int GetNValidTrackerHits ();
    float GetTrackChi2 ();
    float GetTrackNDoF ();
    float GetECalIsoDep ();
    float GetHCalIsoDep ();
    bool  GetIsFake ();

    TGenP* GetGenP (size_t const);
    std::vector<TGenP>* GetGenPVector ();
    TGenP* GetClosestGenP ();
    bool GenPMatchesTo (int const, float const, float const, bool const);
    void PrintGenP (float const MaxDeltaR = 999999, float const MinPt = -1);

    TString GetFlavorString ();
    bool PassesSelection(int const);

    void DefaultValues ();


    // operator def for TLepton
    int  operator<(TLepton&);
    bool cmp (TLepton a, TLepton b) {
      return a.Perp() < b.Perp();
    };

  public:
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


    bool IsFlavor (LeptonFlavor const);


  public:
    // Root Likes ClassDef and ClassImp.
    // Comment them out if you don't need them.
    // There should NOT be a ; since this is a root macro and not a function
    // ClassDef must be the last line of the class before the };
    ClassDef(TLepton,4) // TLepton class
};


#endif
