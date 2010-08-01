////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Jul 12 11:33:18 PDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TString.h"

#include "VgammaAna/TaiwanAna/interface/VgAna.h"

int RunVgAna (TString const& OutFileName, std::vector<TString> const& InFiles)
{
  VgAna MyAna;

  MyAna.AddInFiles(InFiles);
  MyAna.SetOutFile(OutFileName);

  MyAna.Loop();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName] [InFileName]s" << std::endl;
    return 1;
  }

  TString const OutFileName = argv[1];
  printf("OutFile: %s\n", OutFileName.Data());

  std::vector<TString> InFiles;
  for (int i = 2; i < argc; ++i) {
    InFiles.push_back( argv[i] );
    printf("InFile:  %s\n", argv[i]);
  }

  RunVgAna(OutFileName, InFiles);

  return 0;
}
