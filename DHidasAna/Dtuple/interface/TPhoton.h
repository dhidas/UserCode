#ifndef GUARD_TPhoton_h
#define GUARD_TPhoton_h

#include "TLorentzVector.h"


class TPhoton : public TLorentzVector
{
  public:
    TPhoton ();
    ~TPhoton ();

    void SetTrkIso (float const);
    void SetCalIso (float const);
    void SetHCalOverECal (float const);

    float GetTrkIso ();
    float GetCalIso ();
    float GetHCalOverECal ();


  private:
    float TrkIso;
    float CalIso;
    float HCalOverECal;



  public:
    ClassDef(TPhoton, 1)
};













#endif
