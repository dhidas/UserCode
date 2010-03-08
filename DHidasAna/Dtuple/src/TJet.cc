// Author: Dean Andrew Hidas <http://www-cdf.fnal.gov/~dhidas/>

////////////////////////////////////////////////////////////////
//
// TJet
//
// This is a class to hold information about Jets.  It
// does so by storing some information to private variables
// and then has methods to access those variables as well as
// methods which combine these basic variables into other quantaties
// which one might be interested in.  It is intended to be used in
// a TDtuple.  One should fill TJet object like so:
//
////////////////////////////////////////////////////////////////


#include "DHidasAna/Dtuple/interface/TJet.h"



// Root Likes ClassDef and ClassImp.
// Comment them out if you don't need them.
// There should NOT be a ; since this is a root macro and not a function
ClassImp(TJet)
  
//
// Default Constructor
// 
TJet::TJet ()
{
  // Default constructor:  This will set all variables with
  // a dummy value (-9999)

  DefaultValues();
}






//
// Default Destructor
// 
TJet::~TJet ()
{
  // Destructor!
}












float TJet::GetDetEta ()
{
  // Get detector eta of this jet from DetEta variable

  return DetEta;
}






void TJet::SetDetEta (float const in)
{
  // Set the dtuple variable for detector eta: DetEta

  DetEta = in;

  return;
}






void TJet::SetEScaleFactor (float const in)
{
  // Set the dtuple variable for the jet energy scale factor: EScaleFactor

  EScaleFactor = in;

  return;
}






float TJet::GetEScaleFactor ()
{
  // Get the energy scale factor from dtuple variable EScaleFactor

  return EScaleFactor;
}






void TJet::SetEmF (float const in)
{
  // Set the dtuple variable for the Em fraction: EmFrac

  EmF = in;

  return;
}






void TJet::SetHadF (float const in)
{
  // Set the dtuple variable for the Had fraction

  HadF = in;

  return;
}






float TJet::GetEmF ()
{
  // Get the EM fraction from dtuple variable EmFrac

  return EmF;
}






float TJet::GetHadF ()
{
  // Get the Had fraction from dtuple variable HadFrac

  return HadF;
}






float TJet::GetSysPos ()
{
  // Get the positive systematic from the applied jet corrections - SysPos

  return SysPos;
}






void TJet::SetSysPos (float const in)
{
  // Set the positive systematic from the applied jet corrections - SysPos

  SysPos = in;

  return;
}






float TJet::GetSysNeg ()
{
  // Get the negative systematic from the applied jet corrections - SysNeg

  return SysNeg;
}






void TJet::SetSysNeg (float const in)
{
  // Set the negative systematic from the applied jet corrections - SysNeg

  SysNeg = in;

  return;
}






void TJet::DefaultValues()
{
  // Set dummy values to all of the TJet variables.
  // This function is called whenever a new TJet object
  // is created.  All values are set to -9999

  // Set the dummy values
  DetEta        = -999999;
  EScaleFactor  = -999999;
  EmF           = -999999;
  HadF          = -999999;
  SysPos        = -999999;
  SysNeg        = -999999;
}






void TJet::ReScale (float const in)
{
  // Recalse the energy using the scale factor given.

  SetE(  E()  * in );
  SetPx( Px() * in );
  SetPy( Py() * in );
  SetPz( Pz() * in );

  SetEScaleFactor( GetEScaleFactor() * in );

  return;
}





int TJet::operator<(TJet &rhs)
{
  // < operator: Based on Pt.

  if (this->Perp() < rhs.Perp()) return 1;
  return 0;
}





