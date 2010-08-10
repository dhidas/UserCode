////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Aug  9 15:58:52 PDT 2010
//
////////////////////////////////////////////////////////////////////

#include "VgammaAna/LHEReader/interface/LHEEvent.h"


#include <iostream>


LHEEvent::LHEEvent (TString const InFileName)
{
  LHEFile = new std::ifstream(InFileName.Data());
  LHEFileName = InFileName;
  if (!LHEFile->is_open()) {
    std::cerr << "ERROR: cannot open file " << InFileName << std::endl;
  }
}



LHEEvent::~LHEEvent ()
{
  if (LHEFile != 0x0) {
    LHEFile->close();
  }
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


  return 0;
}
