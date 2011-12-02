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
  TH1F* hJetPtMu[NMaxJets];
  TH1F* hJetEtaMu[NMaxJets];
  TH1F* hJetPtEl[NMaxJets];
  TH1F* hJetEtaEl[NMaxJets];
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPtMu[i] = new TH1F( TString::Format("JetPt%iMu", i), TString::Format("JetPt%iMu", i), 50, 0, 500);
    hJetEtaMu[i] = new TH1F( TString::Format("JetEta%iMu", i), TString::Format("JetEta%iMu", i), 50, -4, 4);
    hJetPtEl[i] = new TH1F( TString::Format("JetPt%iEl", i), TString::Format("JetPt%iEl", i), 50, 0, 500);
    hJetEtaEl[i] = new TH1F( TString::Format("JetEta%iEl", i), TString::Format("JetEta%iEl", i), 50, -4, 4);
  }
  TH1F* LeptonPtMu = new TH1F("LeptonPtMu","LeptonPtMu", 50, 0, 500);
  TH1F* LeptonPtEl= new TH1F("LeptonPtEl","LeptonPtEl", 50, 0, 500);


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
	if(muSEL > 1 && eSEL==0){
        hJetPtMu[ijet]->Fill(Jet.Pt(), weight);
        hJetEtaMu[ijet]->Fill(Jet.Eta(), weight);
	LeptonPtMu->Fill(leptonPtRec,weight);
	}
	if(eSEL > 1 && muSEL ==0){
	  hJetPtEl[ijet]->Fill(Jet.Pt(), weight);
	  hJetEtaEl[ijet]->Fill(Jet.Eta(), weight);
	  LeptonPtEl->Fill(leptonPtRec,weight);
        }

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
    hJetPtMu[i]->Write();
    hJetEtaMu[i]->Write();
    hJetPtEl[i]->Write();
    hJetEtaEl[i]->Write();
  }
  LeptonPtMu->Write();
  LeptonPtEl->Write();
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





