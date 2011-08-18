////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Aug 10 16:48:00 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

int PlotAndFit (TString const InFileName)
{
  // Open input file
  std::ifstream f(InFileName.Data());
  if (!f) {
    throw;
  }

  std::vector< std::pair<float, float> > Values;

  std::pair<float, float> p;
  while (!f.eof()) {
    f >> p.first >> p.second;
    Values.push_back(p);
  }

  int const N = (int) Values.size();

  TGraph g(N);
  for (size_t i = 0; i != Values.size(); ++i) {
    g.SetPoint(i, Values[i].first, Values[i].second);
  }


  TCanvas c;
  c.cd();
  g.Draw("Al*");
  g.Fit("pol3");
  c.SaveAs("Fit.eps");


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  PlotAndFit(argv[1]);

  return 0;
}
