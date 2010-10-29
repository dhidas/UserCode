////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Sun Sep 19 12:20:23 PDT 2010
//
////////////////////////////////////////////////////////////////////

#include "StackAll.h"

#include <iostream>


int RunStackAll (TString const InputFileName, TString const OutputFileName)
{
  StackAll Stacker(InputFileName, OutputFileName);
  Stacker.Run();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InputFile] [OutputFileName]" << std::endl;
    return 1;
  }

  TString const InputFileName  = argv[1];
  TString const OutputFileName = argv[2];

  RunStackAll(InputFileName, OutputFileName);

  return 0;
}
