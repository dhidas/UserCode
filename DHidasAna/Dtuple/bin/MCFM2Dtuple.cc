////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Nov 18 19:24:41 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TString.h"
#include "TChain.h"
#include "TLeafD.h"
#include "TMath.h"
#include "TVector2.h"

#include "DtupleWriter.h"

int MCFM2Dtuple (TString const InFileName, TString const OutFileName)
{
  // Convert as best you can the mcfm ntuple into dtuple

  // Make a chain..
  TChain Chain("SomeTree");
  Chain.Add(InFileName);

  // Make a new DtupleWriter Object
  DtupleWriter DT(OutFileName);

  // The event struct
  Dtuple::Event_Struct Ev = DT.fEvent;


  for (int ientry = 0; Chain.GetEntry() != 0; ++ientry) {
    Ev.Lepton_Px[0] = Chain.GetLeaf("P4_Px")->GetValue();
    Ev.Lepton_Py[0] = Chain.GetLeaf("P4_Py")->GetValue();
    Ev.Lepton_Pz[0] = Chain.GetLeaf("P4_Pz")->GetValue();
    Ev.Lepton_Pt[0] = TMath::Sqrt( TMath::Power(Chain.GetLeaf("P4_Px")->GetValue(), 2)
                                 + TMath::Power(Chain.GetLeaf("P4_Py")->GetValue(), 2));

    Ev.Photon_Px[0] = Chain.GetLeaf("P5_Px")->GetValue();
    Ev.Photon_Py[0] = Chain.GetLeaf("P5_Py")->GetValue();
    Ev.Photon_Pz[0] = Chain.GetLeaf("P5_Pz")->GetValue();
    Ev.Photon_Pt[0] = TMath::Sqrt( TMath::Power(Chain.GetLeaf("P5_Px")->GetValue(), 2)
                                 + TMath::Power(Chain.GetLeaf("P5_Py")->GetValue(), 2));

    Ev.Run = 0;
    Ev.Event = ientry;
    Ev.MetMag = TMath::Sqrt( TMath::Power(Chain.GetLeaf("P3_Px")->GetValue(), 2)
                           + TMath::Power(Chain.GetLeaf("P3_Py")->GetValue(), 2));
    Ev.MetPhi = TVector2(Chain.GetLeaf("P3_Px")->GetValue(), Chain.GetLeaf("P3_Px")->GetValue()).Phi();
    Ev.SumEt = Ev.Lepton_Pt[0] + Ev.Photon_Pt[0];

    DT.Fill();
  }


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [OutFile]" << std::endl;
    return 1;
  }

  MCFM2Dtuple(argv[1], argv[2]);

  return 0;
}
