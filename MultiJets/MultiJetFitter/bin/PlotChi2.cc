////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Feb 17 19:52:35 CET 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"

int PlotWidths (TString const InName, float const Column, float const CutBelow)
{
  std::ifstream In(InName.Data());
  if (!In) {
    std::cerr << "ERROR cannot open input file" << std::endl;
    exit(1);
  }

  // Skip first line
  std::string Line;
  std::getline(In, Line);


  float Val;
  float XS, Chi2;
  std::istringstream InLine;

  TH1F h("Chi2", "Chi2", 50, -5, 5);

  while (!In.eof()) {
    std::getline(In, Line);
    InLine.str(Line);

    for (int i = 0; i < Column; ++i) {
      InLine >> Val >> Val;
    }

    InLine >> XS >> Chi2;
    printf("XS/Chi2: %12.1f %12.1f\n", XS, Chi2);

    if (XS > CutBelow) {
      // Fill Hist
      h.Fill(Chi2);
    }

    if (!In.peek()) break;
  }

  TCanvas Can;
  Can.cd();
  h.Draw("hist");
  Can.SaveAs("Chi2.eps");



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  float const Column   = 0;
  float const CutBelow = 0;

  PlotWidths("TestFile", Column, CutBelow);

  return 0;
}
