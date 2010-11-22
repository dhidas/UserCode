////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Oct 12 12:18:33 EDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"

int FillLandauGaus ()
{
  TFile f("Data.root", "create");
  if (!f.IsOpen()) {
    std::cerr << "ERROR: cannot open output file" << std::endl;
    return -1;
  }

  TRandom r(12043);

  float const LMPV = 90;
  float const LSigma = 15;
  int   const NData = 1440;
  int   const NPE = 10;

  TH1F h("data", "data", 45, 0, 450);
  for (int i = 0; i != NData; ++i) {
    h.Fill( r.Landau(LMPV, LSigma) );
  }

  for (int i = 0; i != 125; ++i) {
    //h.Fill( r.Gaus(210, 5) );
  }

  for (int i = 0; i != 75; ++i) {
    //h.Fill( r.Gaus(330, 10) );
  }

  h.Write();


  char BUFF[100];
  for (int ip = 0; ip != NPE; ++ip) {
    if (ip % 1000 == 0) {
      std::cout << "Filling PE: " << ip << std::endl;
    }

    sprintf(BUFF, "Pseudo_LandauOnly_%i", ip);
    TH1F h2(BUFF, BUFF, 45, 0, 450);

    int const N = r.Poisson(NData);
    for (int i = 0; i != N; ++i) {
      h2.Fill( r.Landau(LMPV, LSigma) );
    }

    h2.Write();

  }




  f.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  return FillLandauGaus();
}
