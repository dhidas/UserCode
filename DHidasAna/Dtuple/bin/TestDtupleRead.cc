////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Mar  1 14:08:15 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TDtuple.h"

#include "TString.h"
#include "TChain.h"

int TestDtupleRead (TString const FileName)
{
  TChain Chain("dtuple");
  TDtuple Dtuple(&Chain);

  for (long unsigned int ientry = 0; Dtuple.GetEntry(ientry); ++ientry) {
    if (ientry % 1000 == 0) {
      std::cout << "Events Processed: " << ientry << std::endl;
    }

  }


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [DtupleFile]" << std::endl;
    return 1;
  }

  TestDtupleRead(argv[1]);

  return 0;
}
