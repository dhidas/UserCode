////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Nov 30 14:37:31 EST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "../interface/DHRunTracker.h"
#include "MakeSomePlots.h"







int MakeSomePlots (TString const OutName, std::vector<TString> InFiles)
{

  // Open output file
  TFile OutFile(OutName, "create");
  if (!OutFile.IsOpen()) {
    std::cerr << "Change this if you want, but be careful..." << std::endl;
    return -1;
  }
  OutFile.cd();

  // Duplicate runs thingy
  DHRunTracker RunEvt;


  // Define histograms and set max number of jets to look at..
  int const NMaxJets = 7;
  TH1F* hJetPt[NMaxJets];
  TH1F* hJetEta[NMaxJets];
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPt[i] = new TH1F( TString::Format("JetPt%i", i), TString::Format("JetPt%i", i), 50, 0, 500);
    hJetEta[i] = new TH1F( TString::Format("JetEta%i", i), TString::Format("JetEta%i", i), 50, -4, 4);
  }


  // Define our chain of files/trees
  TChain Chain("micro", "micro");
  for (size_t i = 0; i != InFiles.size(); ++i) {
    Chain.Add(InFiles[i]);
  }

  // Init all branches
  initMicroNtuple((TTree*) &Chain);

  // Loop over all events in chain
  for (int ientry = 0; Chain.GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 1000 == 0) {
      std::cout << "ientry: " << ientry << std::endl;
    }

    // Check for duplicate runs
    if (RunEvt.IsDuplicate(runNumber, eventNumber)) {
      std::cout << "Duplicate event!!  Skiping event: " << runNumber << "   " << eventNumber << std::endl;
      continue;
    }

    // Loop over the jets
    float LastJetPt;
    for (int ijet = 0; ijet != numberOfJets; ++ijet) {
      TLorentzVector Jet(JetPx[ijet], JetPy[ijet], JetPz[ijet], JetE[ijet]);

      if (ijet == 0) {
        LastJetPt = Jet.Pt() + 1;
      }

      // if this jet is less than max number of jets histogram it
      if (ijet < NMaxJets) {
        hJetPt[ijet]->Fill(Jet.Pt(), weight);
        hJetEta[ijet]->Fill(Jet.Eta(), weight);

        // This is just a sanity check...
      }
      if (Jet.Pt() > LastJetPt) {
        std::cerr << "I don't know what's going on with the jets, but they need to be sorted" << std::endl;
      }
      LastJetPt = Jet.Pt();
    }

  }

  OutFile.cd();
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPt[i]->Write();
    hJetEta[i]->Write();
  }

  OutFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutName.root] [InFile]s" << std::endl;
    return 1;
  }

  TString const OutName = argv[1];
  std::vector<TString> InFiles;
  for (int i = 2; i < argc; ++i) {
    InFiles.push_back(argv[i]);
  }

  MakeSomePlots(OutName, InFiles);

  return 0;
}





