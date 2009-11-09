////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////

#ifndef GUARD_DtupleReader_h
#define GUARD_DtupleReader_h

#include "DHidasAna/Dtuple/interface/Dtuple.h"

class DtupleReader : public Dtuple
{
  public:
    DtupleReader ();
    DtupleReader (TString const&);
    virtual ~DtupleReader ();

    void InitializeTree ();


};





#endif
