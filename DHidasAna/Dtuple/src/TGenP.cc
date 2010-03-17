#include "DHidasAna/Dtuple/interface/TGenP.h"

ClassImp(TGenP)

TGenP::TGenP ()
{
  Id = 0;
  MotherId = 0;
}



TGenP::~TGenP ()
{
}


int TGenP::GetId ()
{
  return Id;
}


int TGenP::GetMotherId ()
{
  return MotherId;
}


void TGenP::SetId (int const in)
{
  Id = in;
  return;
}


void TGenP::SetMotherId (int const in)
{
  MotherId = in;
  return;
}


bool TGenP::IsId (int const in)
{
  return (Id == in);
}


bool TGenP::IsMotherId (int const in)
{
  return (MotherId == in);
}
