////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 16:40:54 PDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "VgammaAna/LHEReader/interface/LHEPlotter.h"

int RunLHEPlotter (TString const OutFileName, std::vector<TString> const& LHEFileNames)
{
  LHEPlotter MyPlotter(LHEFileNames, OutFileName);
  MyPlotter.Loop();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName.root] [LHEFileName]s" << std::endl;
    return 1;
  }


  TString const OutFileName = argv[1];

  std::vector<TString> LHEFileNames;
  for (int i = 2; i < argc; ++i) {
    LHEFileNames.push_back(argv[i]);
  }

  RunLHEPlotter(OutFileName, LHEFileNames);

  return 0;
}
