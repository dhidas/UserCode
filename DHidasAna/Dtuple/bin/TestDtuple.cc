////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Fri Oct 23 18:03:20 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "DHidasAna/Dtuple/interface/Dtuple.h"


int TestDtuple ()
{
  Dtuple A;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestDtuple();

  return 0;
}
