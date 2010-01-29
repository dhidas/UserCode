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
  // Default constructor
}


DtupleWriter::DtupleWriter (TString const& OutName)
{
  // Constructor which takes the name of the output file.

  // Create new file
  fDtupleFile = new TFile(OutName, "recreate");
  if (fDtupleFile) {
    // Set the tree
    fDtupleTree = new TTree("Dtuple", "Dtuple");
    fDtupleTree->SetDirectory(fDtupleFile);
    fDtupleTree->Branch("Event", &fEvent, fEvent_Format);
  } else {
    std::cerr << "ERROR: cannot open output file " << OutName << std::endl;
  }

  //std::cout << "DtupleWriter::DtupleWriter done." << std::endl;

  // Comment this in if you want to change the default tree size
  //fDtupleTree->SetMaxTreeSize(4000000);

}

DtupleWriter::DtupleWriter (TTree* DtupleTree)
{
  // Constructor which takes the name of the output file.

  if (DtupleTree) {
    // Set the tree
    fDtupleTree = DtupleTree;
    fDtupleTree->Branch("Event", &fEvent, fEvent_Format);
  } else {
    std::cerr << "ERROR: input TTree* is 0x0" << std::endl;
  }

  //std::cout << "DtupleWriter::DtupleWriter done." << std::endl;

  // Comment this in if you want to change the default tree size
  //fDtupleTree->SetMaxTreeSize(4000000);

}

DtupleWriter::~DtupleWriter ()
{
  // Destructor

  // At the end, let's grab the current file, write, and close it
  TFile* File = fDtupleTree->GetCurrentFile();
  File->Write();
  File->Close();
}

