////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Feb  1 10:32:42 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "DHidasAna/Dtuple/interface/SimpleSUSYPlots.h"

int MakeSimpleSUSYPlots (TString const InFile)
{
  SimpleSUSYPlots MySSP(InFile);

  MySSP.Loop();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  MakeSimpleSUSYPlots(argv[1]);

  return 0;
}
