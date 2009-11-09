////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////

#include "DHidasAna/Dtuple/interface/DtupleWriter.h"

#include <iostream>


DtupleWriter::DtupleWriter ()
{
}


DtupleWriter::DtupleWriter (TString const& OutName)
{
  fDtupleFile = new TFile(OutName, "recreate");
  if (fDtupleFile) {
    fDtupleTree = new TTree("Dtuple", "Dtuple");
    fDtupleTree->SetDirectory(fDtupleFile);
    fDtupleTree->Branch("Event", &fEvent, fEvent_Format);
  }

  //std::cout << "DtupleWriter::DtupleWriter done." << std::endl;

  // Comment this in if you want to change the default tree size
  //fDtupleTree->SetMaxTreeSize(4000000);

}

DtupleWriter::~DtupleWriter ()
{
  TFile* File = fDtupleTree->GetCurrentFile();
  File->Write();
  File->Close();
}

