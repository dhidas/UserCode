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
    Ev.Lepton_Px[0] = Chain.GetLeaf("px4")->GetValue();
    Ev.Lepton_Py[0] = Chain.GetLeaf("py4")->GetValue();
    Ev.Lepton_Pz[0] = Chain.GetLeaf("pz4")->GetValue();
    Ev.Lepton_Pt[0] = TMath::Sqrt( TMath::Power(Chain.GetLeaf("px4")->GetValue(), 2)
                                 + TMath::Power(Chain.GetLeaf("py4")->GetValue(), 2));

    Ev.Photon_Px[0] = Chain.GetLeaf("px5")->GetValue();
    Ev.Photon_Py[0] = Chain.GetLeaf("py5")->GetValue();
    Ev.Photon_Pz[0] = Chain.GetLeaf("pz5")->GetValue();
    Ev.Photon_Pt[0] = TMath::Sqrt( TMath::Power(Chain.GetLeaf("px5")->GetValue(), 2)
                                 + TMath::Power(Chain.GetLeaf("py5")->GetValue(), 2));

    Ev.Run = 0;
    Ev.Event = ientry;
    Ev.MetMag = TMath::Sqrt( TMath::Power(Chain.GetLeaf("px3")->GetValue(), 2)
                           + TMath::Power(Chain.GetLeaf("py3")->GetValue(), 2));
    Ev.MetPhi = TVector2(Chain.GetLeaf("px3")->GetValue(), Chain.GetLeaf("px3")->GetValue()).Phi();
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
