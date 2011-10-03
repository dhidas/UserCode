////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Fri Sep 16 13:14:17 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "RooWorkspace.h"
#include "TFile.h"
#include "TString.h"
#include "RooStats/ModelConfig.h"

int PrintWorkspace (TString const InFileName)
{
  TFile f(InFileName, "read");
  RooWorkspace* m = (RooWorkspace*) f.Get("ws");
  m->Print();


  RooStats::ModelConfig* mc = (RooStats::ModelConfig*) m->obj("ModelConfigSB");
  mc->Print();
  mc = (RooStats::ModelConfig*) m->obj("ModelConfigBG");
  mc->Print();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  PrintWorkspace(argv[1]);

  return 0;
}
