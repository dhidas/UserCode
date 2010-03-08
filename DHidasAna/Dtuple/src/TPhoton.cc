#include "DHidasAna/Dtuple/interface/TPhoton.h"

ClassImp(TPhoton)

TPhoton::TPhoton ()
{
}


TPhoton::~TPhoton ()
{
}


void TPhoton::SetTrkIso (float const in)
{
  TrkIso = in;
  return;
}


float TPhoton::GetTrkIso ()
{
  return TrkIso;
}



void TPhoton::SetCalIso (float const in)
{
  CalIso = in;
  return;
}


float TPhoton::GetCalIso ()
{
  return CalIso;
}



void TPhoton::SetHCalOverECal (float const in)
{
  HCalOverECal = in;
  return;
}


float TPhoton::GetHCalOverECal ()
{
  return HCalOverECal;
}
