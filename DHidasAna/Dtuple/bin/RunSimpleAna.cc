////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Mar  8 14:59:30 CET 2010
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include "DHidasAna/Dtuple/interface/SimpleAna.h"

#include "TChain.h"
#include "TString.h"

int RunSimpleAna (TString const ProcName, std::vector<TString> const& FileNames)
{
  TChain Chain("dtuple");
  for (size_t i = 0; i != FileNames.size(); ++i) {
    Chain.Add(FileNames[i]);
  }

  SimpleAna Ana(ProcName, &Chain);
  Ana.SetFakeRateFile("FakeRates.root");
  Ana.Loop();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [ProcName] [InFile]s" << std::endl;
    return 1;
  }

  TString const ProcName = argv[1];

  std::vector<TString> FileNames;
  for (int i = 2; i < argc; ++i) {
    FileNames.push_back(argv[i]);
  }

  RunSimpleAna(ProcName, FileNames);

  return 0;
}
