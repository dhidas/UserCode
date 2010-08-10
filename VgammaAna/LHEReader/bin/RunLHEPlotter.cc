////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 16:40:54 PDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "VgammaAna/LHEReader/interface/LHEPlotter.h"

int RunLHEPlotter (TString const LHEFileName, TString const OutFileName)
{
  LHEPlotter MyPlotter(LHEFileName, OutFileName);
  MyPlotter.Loop();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [LHEFileName] [OutFileName.root]" << std::endl;
    return 1;
  }

  TString const LHEFileName = argv[1];
  TString const OutFileName = argv[2];

  RunLHEPlotter(LHEFileName, OutFileName);

  return 0;
}
