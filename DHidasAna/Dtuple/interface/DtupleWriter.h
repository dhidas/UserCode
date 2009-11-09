////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////

#ifndef GUARD_DtupleWriter_h
#define GUARD_DtupleWriter_h

#include "DHidasAna/Dtuple/interface/Dtuple.h"

class DtupleWriter : public Dtuple
{
  public:
    DtupleWriter ();
    DtupleWriter (TString const&);
    virtual ~DtupleWriter ();

};











#endif
