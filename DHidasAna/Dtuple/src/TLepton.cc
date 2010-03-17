// Author: Dean Andrew Hidas <http://www-cdf.fnal.gov/~dhidas/>

////////////////////////////////////////////////////////////////
//
// TLepton
//
// The TLepton Class!
//
////////////////////////////////////////////////////////////////

#include "DHidasAna/Dtuple/interface/TLepton.h"

#include <iostream>
#include <algorithm>
#include <iomanip>


// Root Likes ClassDef and ClassImp.
// Comment them out if you don't need them.
// There should NOT be a ; since this is a root macro and not a function
ClassImp(TLepton)


//
// Default Constructor
// 
TLepton::TLepton ()
{
  // Default constructor:  Will call SetDummyValues to 
  // set the initial TLepton values(dummy values)
  DefaultValues();
}






//
// Default Destructor
// 
TLepton::~TLepton ()
{
  // Destructor!
}




float TLepton::GetTrkPt ()
{
  return TrkPt;
}


void TLepton::SetTrkPt (float const in)
{
  TrkPt = in;
  return;
}


float TLepton::Getdxy ()
{
  return dxy;
}


void TLepton::Setdxy (float const in)
{
  dxy = in;
  return;
}


float TLepton::Getdz ()
{
  return dz;
}


void TLepton::Setdz (float const in)
{
  dz = in;
  return;
}


float TLepton::GetZ0 ()
{
  return Z0;
}


void TLepton::SetZ0 (float const in)
{
  Z0 = in;
  return;
}


int TLepton::GetCharge ()
{
  return Charge;
}


void TLepton::SetCharge (int const in)
{
  Charge = in;
  return;
}


int TLepton::GetFlavor ()
{
  return Flavor;
}


void TLepton::SetFlavor (int const in)
{
  Flavor = in;
  return;
}


float TLepton::GetTrkIso ()
{
  return TrkIso;
}


void TLepton::SetTrkIso (float const in)
{
  TrkIso = in;
  return;
}


float TLepton::GetCalIso ()
{
  return CalIso;
}


void TLepton::SetCalIso (float const in)
{
  CalIso = in;
  return;
}


float TLepton::GetECalIso ()
{
  return ECalIso;
}


void TLepton::SetECalIso (float const in)
{
  ECalIso = in;
  return;
}


float TLepton::GetHCalIso ()
{
  return HCalIso;
}


void TLepton::SetHCalIso (float const in)
{
  HCalIso = in;
  return;
}


float TLepton::GetCalE ()
{
  return CalE;
}


void TLepton::SetCalE (float const in)
{
  CalE = in;
  return;
}


float TLepton::GetHCalOverECal ()
{
  return HCalOverECal;
}


void TLepton::SetHCalOverECal (float const in)
{
  HCalOverECal = in;
  return;
}


float TLepton::GetEoverPin ()
{
  return EoverPin;
}


void TLepton::SetEoverPin (float const in)
{
  EoverPin = in;
  return;
}


float TLepton::GetfBrem ()
{
  return fBrem;
}


void TLepton::SetfBrem (float const in)
{
  fBrem = in;
  return;
}


int TLepton::GetIsConvertedPhoton ()
{
  return IsConvertedPhoton;
}


void TLepton::SetIsConvertedPhoton (int const in)
{
  IsConvertedPhoton = in;
  return;
}


int TLepton::GetPassSelection ()
{
  return PassSelection;
}


void TLepton::SetPassSelection (int const in)
{
  PassSelection = in;
  return;
}


int TLepton::GetDetector ()
{
  return Detector;
}


void TLepton::SetDetector (int const in)
{
  Detector = in;
  return;
}


int TLepton::GetClassification ()
{
  return Classification;
}


void TLepton::SetClassification (int const in)
{
  Classification = in;
  return;
}


float TLepton::GetSigmaIEtaIEta ()
{
  return SigmaIEtaIEta;
}


void TLepton::SetSigmaIEtaIEta (float const in)
{
  SigmaIEtaIEta = in;
  return;
}


float TLepton::GetDeltaEtaIn ()
{
  return DeltaEtaIn;
}


void TLepton::SetDeltaEtaIn (float const in)
{
  DeltaEtaIn = in;
  return;
}


float TLepton::GetDeltaPhiIn ()
{
  return DeltaPhiIn;
}


void TLepton::SetDeltaPhiIn (float const in)
{
  DeltaPhiIn = in;
  return;
}


float TLepton::GetE2x5overE5x5 ()
{
  return E2x5overE5x5;
}


void TLepton::SetE2x5overE5x5 (float const in)
{
  E2x5overE5x5 = in;
  return;
}


float TLepton::GetConvDist ()
{
  return ConvDist;
}


void TLepton::SetConvDist (float const in)
{
  ConvDist = in;
  return;
}


float TLepton::GetConvdCotTheta ()
{
  return ConvdCotTheta;
}


void TLepton::SetConvdCotTheta (float const in)
{
  ConvdCotTheta = in;
  return;
}


int TLepton::GetNValidTrackerHits ()
{
  return NValidTrackerHits;
}


void TLepton::SetNValidTrackerHits (int const in)
{
  NValidTrackerHits = in;
  return;
}


float TLepton::GetTrackChi2 ()
{
  return TrackChi2;
}


void TLepton::SetTrackChi2 (float const in)
{
  TrackChi2 = in;
  return;
}


float TLepton::GetTrackNDoF ()
{
  return TrackNDoF;
}


void TLepton::SetTrackNDoF (float const in)
{
  TrackNDoF = in;
  return;
}


float TLepton::GetECalIsoDep ()
{
  return ECalIsoDep;
}


void TLepton::SetECalIsoDep (float const in)
{
  ECalIsoDep = in;
  return;
}


float TLepton::GetHCalIsoDep ()
{
  return HCalIsoDep;
}


void TLepton::SetHCalIsoDep (float const in)
{
  HCalIsoDep = in;
  return;
}



TGenP* TLepton::GetGenP (size_t const i)
{
  return &GenP[i];
}




std::vector<TGenP>* TLepton::GetGenPVector ()
{
  return &GenP;
}




TString TLepton::GetFlavorString ()
{
  if (Flavor == kLeptonFlavor_Electron) {
    return "e";
  } else if (Flavor == kLeptonFlavor_Muon) {
    return "m";
  } else if (Flavor == kLeptonFlavor_Tau) {
    return "t";
  }

  return "";
}


bool TLepton::PassesSelection(int const Sel)
{
  if ( PassSelection & (0x1 << Sel) ) {
    return true;
  }

  return false;
}
















void TLepton::DefaultValues ()
{
  TrkPt             = -999999;
  dxy               = -999999;
  dz                = -999999;
  Z0                = -999999;
  Charge            = -999999;
  Flavor            = -999999;
  TrkIso            = -999999;
  CalIso            = -999999;
  ECalIso           = -999999;
  HCalIso           = -999999;
  CalE              = -999999;
  HCalOverECal      = -999999;
  EoverPin          = -999999;
  fBrem             = -999999;
  IsConvertedPhoton = -999999;
  PassSelection     = 0x0;
  Detector          = 0x0;
  Classification    = 0x0;
  SigmaIEtaIEta     = -999999;
  DeltaEtaIn        = -999999;
  DeltaPhiIn        = -999999;
  E2x5overE5x5      = -999999;
  ConvDist          = -999999;
  ConvdCotTheta     = -999999;

}



bool TLepton::IsFlavor (LeptonFlavor const i)
{
  if (i == GetFlavor()) {
    return true;
  }
  return false;
}



int TLepton::operator<(TLepton &rhs)
{
    // The < operator is defined to work on Pt().
    // See source for more info

    if (this->Perp() < rhs.Perp()) return 1;
      return 0;
}

