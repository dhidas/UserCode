////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon May  2 04:36:37 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TString.h"
#include "TFile.h"

#include "DHidasLJAna/LeptonPlusJets/interface/Dtuple.h"
#include "DHidasLJAna/LeptonPlusJets/interface/TAnaHist.h"





void PlotDileptonEvents (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  static TAnaHist Hist(&OutFile, "DileptonEvents");

  if (Ev.NLeptons == 1) {
    Hist.FillTH1F("Single", "DileptonMass", "Dilepton Mass", "Events", 50, 0, 500, Ev.Lep[0].Pt());
  }
  if (Ev.NLeptons != 2) {
    return;
  }

  //std::cout << (Ev.Lep[0] + Ev.Lep[1]).M() << std::endl;
  Hist.FillTH1F("DileptonMass", "Dilepton Mass", "Dilepton Mass", "Events", 250, 0, 500, (Ev.Lep[0]+Ev.Lep[1]).M());
  return;
}




int RunLJDtuple (TString const OutFileName, std::vector<TString> const& InFileNames)
{
  // Open OutFile
  TFile OutFile(OutFileName, "create");
  if (!OutFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
    throw;
  }


  // Grab small ntuple.  Everything is contained in "Ev"
  Dtuple DT(InFileNames);
  Dtuple::SimpleEvent& Ev = DT.GetEvt();
  for (long long ientry = 0; DT.GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 10000 == 0) {
      std::cout << "Processing: " << ientry << std::endl;
    }

    // Analysis starts here!
    PlotDileptonEvents(Ev, OutFile);




    DT.ClearDtuple();
  }

  OutFile.Write();
  OutFile.Close();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName] [InFileName]s" << std::endl;
    return 1;
  }

  TString const OutFileName = argv[1];
  std::vector<TString> InFileNames;
  for (int i = 2; i < argc; ++i) {
    InFileNames.push_back(argv[i]);
  }

  RunLJDtuple(OutFileName, InFileNames);

  return 0;
}
