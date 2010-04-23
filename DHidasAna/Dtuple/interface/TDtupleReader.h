#ifndef GUARD_TDtupleReader_h
#define GUARD_TDtupleReader_h

#include "DHidasAna/Dtuple/interface/TDtuple.h"

class TDtupleReader : public TDtuple
{
  public:
    TDtupleReader (TChain*);
    ~TDtupleReader ();

    void Loop (long unsigned int Max = 0);
    void ObjectCleaning ();

    virtual void BeginJob () = 0;
    virtual void Analyze (long unsigned int const) = 0;
    virtual void EndJob () = 0;


};






































#endif
