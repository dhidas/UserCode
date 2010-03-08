#ifndef GUARD_TDUtility_h
#define GUARD_TDUtility_h

#include <set>

#include "TObject.h"

class TDUtility : public TObject
{
  public:
    TDUtility ();
    virtual ~TDUtility ();

    bool IsDuplicateEvent (int const, int const);

  public:
    ClassDef(TDUtility,1) // TDUtility class
};













#endif
