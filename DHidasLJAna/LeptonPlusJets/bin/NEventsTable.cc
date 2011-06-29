////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Jun 13 11:01:29 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"


int NEventsTable ()
{
  TString Files[9] = {
    "ZprimeM500PAT.out",
    "ZprimeM750PAT.out",
    "ZprimeM1000PAT.out",
    "ZprimeM1250PAT.out",
    "ZprimeM1500PAT.out",
    "ZprimeM2000PAT.out",
    "ZprimeM3000PAT.out",
    "ZprimeM4000PAT.out",
    "TTJets.out"
  };

  int Mass[9] = {
     500,
     750,
    1000,
    1250,
    1500,
    2000,
    3000,
    4000,
     175
  };

  int const DiagCut = 100;

  float X[9];
  float Y[9];

  for (int i = 0; i != 9; ++i) {
    std::ifstream InFile(Files[i].Data());
    if (!InFile.is_open()) {
      std::cerr << "ERROR: cannot open file: " << Files[i] << std::endl;
      throw;
    }
    TString Name;
    int Cut;
    float Acceptance, NEvents;
    while (!InFile.eof()) {
      InFile >> Name >> Cut >> Acceptance >> NEvents;
      if (Cut == DiagCut) {
        X[i] = Mass[i];
        Y[i] = NEvents;
        break;
      }
    }
    InFile.close();
  }

  float SigY[8];
  for (int i = 0; i != 9; ++i) {
    if (i < 8) {
      SigY[i] = Y[i] / TMath::Sqrt(Y[i] + Y[8]);
    }
    printf("Mass = %4i  NExpected = %9.3E  sig/sqrt(sig+bg) = %9.3E\n", (int) X[i], Y[i], Y[i] / TMath::Sqrt(Y[i] + Y[8]));
  }

  TGraph g(8, X, Y);
  TCanvas c;
  c.cd();
  c.SetLogy(1);
  g.SetTitle( TString::Format("Expected Number of Events in 1/fb for DiagCut=%i", DiagCut) );
  g.GetXaxis()->SetTitle("Z' Mass (GeV)");
  g.GetYaxis()->SetTitle("Number of Events");
  g.Draw("AC*");
  c.SaveAs("Expected.gif");

  TGraph g2(8, X, SigY);
  TCanvas c2;
  c2.cd();
  c2.SetLogy(1);
  g2.SetTitle( TString::Format("Significance for 1/fb and DiagCut=%i", DiagCut) );
  g2.GetXaxis()->SetTitle("Z' Mass (GeV)");
  g2.GetYaxis()->SetTitle("S / #sqrt{S+B}");
  g2.Draw("AC*");
  c2.SaveAs("Significance.gif");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  NEventsTable();

  return 0;
}
