////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Mar  8 14:59:30 CET 2010
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cstdlib>

#include "DHidasAna/SimpleAna/interface/SimpleAna.h"

#include "TChain.h"
#include "TString.h"

int RunSimpleAna (TString const ProcName, std::vector<TString> const& FileNames, bool const RunFakes)
{
  TChain Chain("dtuple");
  for (size_t i = 0; i != FileNames.size(); ++i) {
    Chain.Add(FileNames[i]);
  }

  SimpleAna Ana(ProcName, &Chain);
  Ana.SetFakeRateFile("FakeRates.root");
  Ana.RunFakes(RunFakes);
  Ana.Loop();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " [RunFakes 0|1] [ProcName] [InFile]s" << std::endl;
    return 1;
  }

  bool const RunFakes = atoi( argv[1] ) == 0 ? false : true;
  TString const ProcName = argv[2];

  std::vector<TString> FileNames;
  for (int i = 3; i < argc; ++i) {
    FileNames.push_back(argv[i]);
  }

  RunSimpleAna(ProcName, FileNames, RunFakes);

  return 0;
}
