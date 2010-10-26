////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 15:58:52 PDT 2010
//
////////////////////////////////////////////////////////////////////

#include "VgammaAna/LHEReader/interface/LHEEvent.h"


#include <iostream>


LHEEvent::LHEEvent (TString const In)
{
  LHEFile = 0x0;
  LHEFileNames.push_back(In);
  ifile = 0;
  OpenNextFile();
}



LHEEvent::LHEEvent (std::vector<TString> const& In)
{
  LHEFile = 0x0;
  LHEFileNames = In;
  ifile = 0;
  OpenNextFile();
}



LHEEvent::~LHEEvent ()
{
  if (!LHEFile) {
    LHEFile->close();
    delete LHEFile;
  }
}



bool LHEEvent::OpenNextFile ()
{
  if (LHEFile != 0x0) {
    LHEFile->close();
    delete LHEFile;
  }


  if (ifile < LHEFileNames.size()) {
    std::cout << "Attempting to open file " << LHEFileNames[ifile] << std::endl;
    LHEFile = new std::ifstream(LHEFileNames[ifile].Data());
    LHEFileName = LHEFileNames[ifile];
    if (!LHEFile->is_open()) {
      std::cerr << "ERROR: cannot open file " << LHEFileNames[ifile] << std::endl;
      return false;
    }
  } else {
    return false;
  }

  ++ifile;
  return true;
}





int LHEEvent::NextEvent ()
{
  Particles.clear();

  int NParticles;
  float Junk;
  int Id, Mo1, Mo2;
  float Px, Py, Pz, En;
  LHEParticle MyP;
  for (TString Line; Line.ReadLine(*LHEFile); ) {
    if (Line.BeginsWith("<event>")) {
      *LHEFile >> NParticles
               >> Junk
               >> Weight
               >> Junk
               >> Junk
               >> Junk;
      for (int i = 0; i < NParticles; ++i) {
        *LHEFile >> Id
                 >> Junk
                 >> Junk
                 >> Mo1
                 >> Mo2
                 >> Junk
                 >> Px
                 >> Py
                 >> Pz
                 >> En
                 >> Junk
                 >> Junk
                 >> Junk;

        MyP.SetPxPyPzE(Px, Py, Pz, En);
        MyP.Id = Id;
        MyP.Mo1 = Mo1;
        MyP.Mo2 = Mo2;

        Particles.push_back(MyP);
      }
    } else {
      continue;
    }
    return 1;
  }

  if (OpenNextFile()) {
    return 1;
  }

  return 0;
}
