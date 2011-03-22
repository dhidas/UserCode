////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Mar 22 10:10:15 CET 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "TString.h"

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


  std::map<float, std::pair<int, int> > PassFail;


  // for each input file PE
  for (size_t ifile = 0; ifile != PEFileNames.size(); ++ifile) {
    printf("Opening file: %s\n", PEFileNames[ifile].Data());
    std::ifstream InDataFile(PEFileNames[ifile].Data());
    Line.ReadLine(InDataFile);

    while (Line.ReadLine(InDataFile)) {
      InLine.str(Line.Data());
      InLine >> tmp;
      if (tmp == 9999) {
        continue;
      }
      for (size_t i = 0; i != NMasses; ++i) {

        if (tmp >= DataXS[i]) {
          ++PassFail[Masses[i]].first;
        } else {
          ++PassFail[Masses[i]].second;
        }
      }

    }
  }

  for (size_t i = 0; i != NMasses; ++i) {
    printf("Mass: %4i DataXS: %6.2f  Pass/Fail %10i  %10i  p-value: %12E\n", (int) Masses[i], DataXS[i], PassFail[ Masses[i] ].first, PassFail[ Masses[i] ].second,
        ((float) PassFail[ Masses[i] ].first) / ((float) PassFail[ Masses[i] ].first + PassFail[ Masses[i] ].second));
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

