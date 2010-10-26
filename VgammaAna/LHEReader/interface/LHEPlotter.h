////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 16:33:37 PDT 2010
//
////////////////////////////////////////////////////////////////////
#ifndef GUARD_LHEPlotter_h
#define GUARD_LHEPlotter_h

#include "VgammaAna/LHEReader/interface/LHEEvent.h"

#include "TFile.h"

class LHEPlotter : public LHEEvent
{
  public:
    LHEPlotter (TString const, TString const);
    LHEPlotter (std::vector<TString> const&, TString const);
    ~LHEPlotter ();

    int Loop ();

  private:
    TFile* OutFile;
};








#endif
