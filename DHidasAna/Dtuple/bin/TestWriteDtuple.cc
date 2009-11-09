////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Mon Oct 26 21:23:45 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "DHidasAna/Dtuple/interface/DtupleWriter.h"

int TestWriteDtuple ()
{
  std::cout << "Begin TestWriteDtuple" << std::endl;
  DtupleWriter A("TestOutFile.root");

  for (int i = 0; i != 200000; ++i) {
    A.DefaultValues();
    A.Fill();
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestWriteDtuple();

  return 0;
}
