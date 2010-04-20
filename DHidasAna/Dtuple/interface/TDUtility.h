#ifndef GUARD_TDUtility_h
#define GUARD_TDUtility_h

#include <set>

#include "DHidasAna/Dtuple/interface/TLepton.h"

#include "TObject.h"

class TDUtility : public TObject
{
  public:
    TDUtility ();
    virtual ~TDUtility ();

    bool IsDuplicateEvent (int const, int const);
    TString GetLeptonFlavorsString (std::vector<TLepton>&);
    static float GetConversionR(TLorentzVector&, int, float, TLorentzVector&, int, float, float);
    static std::pair<float, float> GetConversionXY(TLorentzVector&, int, float, TLorentzVector&, int, float, float);

  public:
    ClassDef(TDUtility,1) // TDUtility class
};













#endif
