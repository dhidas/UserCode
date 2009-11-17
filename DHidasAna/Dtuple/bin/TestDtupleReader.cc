////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Mon Nov  9 19:25:03 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "DHidasAna/Dtuple/interface/DtupleReader.h"

#include "TH1D.h"

int TestDtupleReader (TString const FileName)
{
  DtupleReader DR(FileName);

  for (unsigned long ientry = 0; DR.GetEntry(ientry) != 0; ++ientry) {
    if (ientry % 1000 == 0) {
      printf("Processing event: %15lu\n", ientry);
    }
  }
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [FileName.root]" << std::endl;
    return 1;
  }

  TestDtupleReader(argv[1]);

  return 0;
}
