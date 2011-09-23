////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Sat Mar 12 13:29:37 CET 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPaveLabel.h"

int GetPValues (TString const DataFileName, std::vector<TString> PEFileNames)
{
  // Read DataFile
  std::vector<float> Masses;
  std::vector<float> DataXS;
  std::istringstream InLine;
  TString Line;
  std::ifstream InDataFile(DataFileName.Data());

  Line.ReadLine(InDataFile);
  InLine.str(Line.Data());
  for (float tmp; InLine >> tmp; ) {
    Masses.push_back(tmp);
  }

  Line.ReadLine(InDataFile);
  std::cout << "hi " << Line << std::endl;
  InLine.clear();
  InLine.str(Line.Data());

  size_t const NMasses = Masses.size();
  printf("NMasses %i\n", (int) NMasses);

  float tmp;
  for (size_t im = 0; im != NMasses; ++im) {
    InLine >> tmp;
    printf("%7E  ", tmp);
    DataXS.push_back(tmp);
  }
  printf("\n");


  std::vector<TH1F*> hD;
  for (size_t imass = 0; imass < NMasses; ++imass) {
    char BUFF[100];
    sprintf(BUFF, "D_%i", (int) Masses[imass]);
    TH1F* h = new TH1F(BUFF, BUFF, 100, 0, 20);
    hD.push_back(h);
  }
  std::map<float, std::pair<int, int> > PassFail;


  // for each input file PE
  for (size_t ifile = 0; ifile != PEFileNames.size(); ++ifile) {
    printf("Opening file: %s\n", PEFileNames[ifile].Data());
    std::ifstream InDataFile(PEFileNames[ifile].Data());
    Line.ReadLine(InDataFile);

    while (Line.ReadLine(InDataFile)) {
      InLine.str(Line.Data());
      for (size_t i = 0; i != NMasses; ++i) {
        InLine >> tmp;
        if (tmp == -9999) {
          continue;
        }

        hD[i]->Fill(tmp);

        if (tmp >= DataXS[i]) {
          ++PassFail[Masses[i]].first;
        } else {
          ++PassFail[Masses[i]].second;
        }
      }

    }
  }

  for (size_t i = 0; i != NMasses; ++i) {
    float const PVal = ((float) PassFail[ Masses[i] ].first) / ((float) PassFail[ Masses[i] ].first + PassFail[ Masses[i] ].second);
    printf("Mass: %4i DataXS: %6.2f  Pass/Fail %10i  %10i  p-value: %12E  Sigma: %8.2f\n",
        (int) Masses[i], DataXS[i], PassFail[ Masses[i] ].first, PassFail[ Masses[i] ].second,
        PVal,
        TMath::Sqrt(2)*TMath::ErfcInverse(PVal));
  }


  for (size_t im = 0; im != NMasses; ++im) {
    TCanvas Can;
    Can.cd();
    hD[im]->Draw("hist");
    TLine MyLine(DataXS[im], 0, DataXS[im], hD[im]->GetMaximum());
    MyLine.SetLineColor(2);
    MyLine.SetLineWidth(2);
    MyLine.Draw("same");
    Can.SetLogy(true);
    char BUFF[100];

    float const PVal = ((float) PassFail[ Masses[im] ].first) / ((float) PassFail[ Masses[im] ].first + PassFail[ Masses[im] ].second);
    sprintf(BUFF, "p-value: %10E = %5.2f#sigma", PVal, (float) TMath::Sqrt(2)*TMath::ErfcInverse(PVal));
    std::cout << BUFF << std::endl;
    TPaveLabel PLabel;
    PLabel.SetLabel(BUFF);
    PLabel.SetX1NDC(0.40);
    PLabel.SetX2NDC(0.90);
    PLabel.SetY1NDC(0.51);
    PLabel.SetY2NDC(0.55);
    PLabel.SetTextSize();
    PLabel.Draw("same");

    sprintf(BUFF, "TestStatisticDist_%i.eps", (int) Masses[im]);
    Can.SaveAs(BUFF);
  }

  return 0;
}



int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [DataFile.dat] [PEFile.dat]s" << std::endl;
    return 1;
  }

  TString const DataFileName = argv[1];
  std::vector<TString> PEFiles;
  for (int i = 2; i < argc; ++i) {
    PEFiles.push_back(argv[i]);
  }

  GetPValues(DataFileName, PEFiles);

  return 0;
}

