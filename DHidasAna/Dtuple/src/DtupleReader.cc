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



DtupleReader::DtupleReader (TString const& FileName)
{
  // Constructor which takes the name of the input file
  // You should be able to use \* with this..

  // Set the TTree
  TChain* Chain = new TChain("Dtuple", "Dtuple");
  Chain->Add(FileName);
  fDtupleTree = (TTree*) Chain;
  InitializeTree();
}



DtupleReader::DtupleReader (std::vector<TString> const& FileNames)
{
  // Constructor which takes a vector of input filenames

  // Set the TTree
  TChain* Chain = new TChain("Dtuple", "Dtuple");
  for (size_t i = 0; i != FileNames.size(); ++i) {
    Chain->Add(FileNames[i]);
  }
  fDtupleTree = (TTree*) Chain;
  InitializeTree();
}



DtupleReader::~DtupleReader ()
{
  // Destructor

  if (fDtupleFile && fDtupleTree) {
    // Well, for now nothing needs to be done..
  }
}



void DtupleReader::InitializeTree ()
{
  // Initialize the tree with the correct event branch
  fDtupleTree->SetBranchAddress("Event", &fEvent);
  return;
}


