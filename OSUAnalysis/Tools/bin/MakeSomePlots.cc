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
  TH1F* hJetPtMu_NoPile[NMaxJets];
  TH1F* hJetEtaMu_NoPile[NMaxJets];
  TH1F* hJetPtEl_NoPile[NMaxJets];
  TH1F* hJetEtaEl_NoPile[NMaxJets];
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPtMu[i] = new TH1F( TString::Format("JetPt%iMu", i), TString::Format("JetPt%iMu", i), 50, 0, 500);
    hJetEtaMu[i] = new TH1F( TString::Format("JetEta%iMu", i), TString::Format("JetEta%iMu", i), 50, -4, 4);
    hJetPtEl[i] = new TH1F( TString::Format("JetPt%iEl", i), TString::Format("JetPt%iEl", i), 50, 0, 500);
    hJetEtaEl[i] = new TH1F( TString::Format("JetEta%iEl", i), TString::Format("JetEta%iEl", i), 50, -4, 4);
    hJetPtMu_NoPile[i] = new TH1F( TString::Format("JetPt%iMu_NoPile", i), TString::Format("JetPt%iMu_NoPile", i), 50, 0, 500);
    hJetEtaMu_NoPile[i] = new TH1F( TString::Format("JetEta%iMu_NoPile", i), TString::Format("JetEta%iMu_NoPile", i), 50, -4, 4);
    hJetPtEl_NoPile[i] = new TH1F( TString::Format("JetPt%iEl_NoPile", i), TString::Format("JetPt%iEl_NoPile", i), 50, 0, 500);
    hJetEtaEl_NoPile[i] = new TH1F( TString::Format("JetEta%iEl_NoPile", i), TString::Format("JetEta%iEl_NoPile", i), 50, -4, 4);

  }
  TH1F* LeptonPtMu = new TH1F("LeptonPtMu","LeptonPtMu", 50, 0, 500);
  TH1F* LeptonPtEl= new TH1F("LeptonPtEl","LeptonPtEl", 50, 0, 500);
  TH1F* LeptonPtMu_NoPile = new TH1F("LeptonPtMu_NoPile","LeptonPtMu_NoPile", 50, 0, 500);
  TH1F* LeptonPtEl_NoPile= new TH1F("LeptonPtEl_NoPile","LeptonPtEl_NoPile", 50, 0, 500);

  TH1F* HTMu = new TH1F("HTMu","HTMu", 50, 0, 2000);
  TH1F* HTEl= new TH1F("HTEl","HTEl", 50, 0, 2000);
  TH1F* HTMu_NoPile = new TH1F("HTMu_NoPile","HTMu_NoPile", 50, 0, 2000);
  TH1F* HTEl_NoPile= new TH1F("HTEl_NoPile","HTEl_NoPile", 50, 0, 2000);

  TH1F* METMu = new TH1F("METMu","METMu", 50, 0, 500);
  TH1F* METEl= new TH1F("METEl","METEl", 50, 0, 500);
  TH1F* METMu_NoPile = new TH1F("METMu_NoPile","METMu_NoPile", 50, 0, 500);
  TH1F* METEl_NoPile= new TH1F("METEl_NoPile","METel_NoPile", 50, 0,500);

  TH1F* HnPileUpVtx= new TH1F("nPileUpVtx","nPileUpVtx", 50, 0, 50);
  TH1F* HnPileUpVtx_NoPile = new TH1F("nPileUpVtx_NoPile","nPileUpVtx_NoPile", 50, 0, 50);

  TH1F*  hNumberOfJetsMu = new TH1F("NumberOfJetsMu","NumberOfJetsMu", 10, 5, 15);
  TH1F*  hNumberOfBJetsMu = new TH1F("NumberOfBJetsMu","NumberOfBJetsMu", 10, 0, 10);
  TH1F*  hNumberOfJetsEl = new TH1F("NumberOfJetsEl","NumberOfJetsEl", 10, 5, 15);
  TH1F*  hNumberOfBJetsEl = new TH1F("NumberOfBJetsEl","NumberOfBJetsEl", 10, 0, 10);

  TH1F* hBjetPtMu  = new TH1F("hBjetPtMu","hBjetPtMu", 50, 0, 500);
  TH1F* hBjetPtEl = new TH1F("hBjetPtEl","hBjetPtEl", 50, 0, 500);
  TH1F* hBjetEtaMu  = new TH1F("hBjetEtaMu","hBjetEtaMu", 50, -4, 4);
  TH1F* hBjetEtaEl = new TH1F("hBjetEtaEl","hBjetEtaEl", 50, -4, 4);
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
    //define weight without PileUp reweighting

    float BtagScaleFactorSSVHEM=1;
    if(type>0)    BtagScaleFactorSSVHEM=0.95;
    std::cout<<weight<<std::endl;
    weight=weight*BtagScaleFactorSSVHEM;
    std::cout<<weight<<std::endl;
    std::cout<<"--------------------"<<std::endl;
    float  weightNoPile=weight/PileUpWeight;
    //make lepton plots
    HnPileUpVtx->Fill(nPileUpVtx,weight);
    HnPileUpVtx_NoPile->Fill(nPileUpVtx,weightNoPile);

    if(HT>0.0){
      TLorentzVector BJet(BJetPx[0], BJetPy[0], BJetPz[0], BJetE[0]);
    if(muSEL > 1 && eSEL==0){

      LeptonPtMu->Fill(leptonPtRec,weight);
      LeptonPtMu_NoPile->Fill(leptonPtRec,weightNoPile);
      HTMu->Fill(HT,weight);
      HTMu_NoPile->Fill(HT,weightNoPile);
      METMu->Fill(MET,weight);
      METMu_NoPile->Fill(MET,weightNoPile);

      hNumberOfJetsMu->Fill(numberOfJets,weight);
      hNumberOfBJetsMu->Fill(numberOfBJets,weight);
      hBjetPtMu->Fill(BJet.Pt(), weight);
      hBjetEtaMu ->Fill(BJet.Eta(), weight);


}
    if(eSEL > 1 && muSEL ==0){

      LeptonPtEl->Fill(leptonPtRec,weight);
      LeptonPtEl_NoPile->Fill(leptonPtRec,weightNoPile);
      HTEl->Fill(HT,weight);
      HTEl_NoPile->Fill(HT,weightNoPile);
      METEl->Fill(MET,weight);
      METEl_NoPile->Fill(MET,weightNoPile);
	hNumberOfJetsEl->Fill(numberOfJets,weight);
      hNumberOfBJetsEl->Fill(numberOfBJets,weight);
      hBjetPtEl->Fill(BJet.Pt(), weight);
      hBjetEtaEl ->Fill(BJet.Eta(), weight);

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
	hJetPtMu_NoPile[ijet]->Fill(Jet.Pt(), weightNoPile);
        hJetEtaMu_NoPile[ijet]->Fill(Jet.Eta(), weightNoPile);
	}
	if(eSEL > 1 && muSEL ==0){
	  hJetPtEl[ijet]->Fill(Jet.Pt(), weight);
	  hJetEtaEl[ijet]->Fill(Jet.Eta(), weight);
	  hJetPtEl_NoPile[ijet]->Fill(Jet.Pt(), weightNoPile);
          hJetEtaEl_NoPile[ijet]->Fill(Jet.Eta(), weightNoPile);


        }

        // This is just a sanity check...
      }
      if (Jet.Pt() > LastJetPt) {
        std::cerr << "I don't know what's going on with the jets, but they need to be sorted" << std::endl;
      }
      LastJetPt = Jet.Pt();
    }

  }
  }
  OutFile.cd();
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPtMu[i]->Write();
    hJetEtaMu[i]->Write();
    hJetPtEl[i]->Write();
    hJetEtaEl[i]->Write();
    //   hJetPtMu_NoPile[i]->Write();
    //hJetEtaMu_NoPile[i]->Write();
    //hJetPtEl_NoPile[i]->Write();
    //hJetEtaEl_NoPile[i]->Write();
  }
  LeptonPtMu->Write();
  LeptonPtEl->Write();
  //LeptonPtMu_NoPile->Write();
  //LeptonPtEl_NoPile->Write();

  HTMu->Write();
  HTEl->Write();
  //HTMu_NoPile->Write();
  //HTEl_NoPile->Write();
  METMu->Write();
  METEl->Write();
  //METMu_NoPile->Write();
  //METEl_NoPile->Write();
  HnPileUpVtx->Write();
  //HnPileUpVtx_NoPile->Write();

  hNumberOfJetsMu->Write();
  hNumberOfBJetsMu->Write();
  hNumberOfJetsEl->Write();
  hNumberOfBJetsEl->Write();
  hBjetPtMu->Write();
  hBjetEtaMu->Write();
  hBjetPtEl->Write();
  hBjetEtaEl->Write();

  std::cout<<"----------------------------------------------------------------"<<std::endl;
  std::cout<<"MuonEvents (noPileUpReweightung): "<<HTMu_NoPile->Integral(0,50)<<" MuonEvents (afterPileUpReweightung): "<<HTMu->Integral(0,50)<<std::endl;
  std::cout<<"ElectronEvents (noPileUpReweightung): "<<HTEl_NoPile->Integral(0,50)<<" ElectronEvents (afterPileUpReweightung): "<<HTEl->Integral(0,50)<<std::endl;
  std::cout<<"----------------------------------------------------------------"<<std::endl;
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





