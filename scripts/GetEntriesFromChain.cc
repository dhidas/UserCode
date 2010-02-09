////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Feb  8 09:12:43 EST 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TChain.h"


int GetEntriesFromChain (TString const TreeName, std::vector<TString> const& Files)
{
  TChain MyChain(TreeName, TreeName);
  for (size_t i = 0; i != Files.size(); ++i) {
    MyChain.Add(Files[i]);
  }
  std::cout << MyChain.GetEntries() << std::endl;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [TreeName] [InputFile]s" << std::endl;
    return 1;
  }

  std::vector<TString> Files;
  for (int i = 2; i < argc; ++i) {
    Files.push_back(argv[i]);
  }

  GetEntriesFromChain(argv[1], Files);

  return 0;
}

