////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@cern.ch>
//
// Created on: Fri Oct 23 14:35:46 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include "DHidasAna/Dtuple/interface/DtupleReader.h"

#include "TChain.h"

DtupleReader::DtupleReader ()
{
  // Default constructor.
  fDtupleFile =0x0;
}

DtupleReader::~DtupleReader ()
{
  if (fDtupleFile && fDtupleTree) {
  }
}



DtupleReader::DtupleReader (TString const& FileName)
{
  // Set the TTree
  TChain* Chain = new TChain("Dtuple", "Dtuple");
  Chain->Add(FileName);
  fDtupleTree = (TTree*) Chain;
}

DtupleReader::DtupleReader (std::vector<TString> const& FileNames)
{
  // Set the TTree
  TChain* Chain = new TChain("Dtuple", "Dtuple");
  for (size_t i = 0; i != FileNames.size(); ++i) {
    Chain->Add(FileNames[i]);
  }
  fDtupleTree = (TTree*) Chain;
}

void DtupleReader::InitializeTree ()
{
  fDtupleTree->SetBranchAddress("Event", &fEvent);
  return;
}


