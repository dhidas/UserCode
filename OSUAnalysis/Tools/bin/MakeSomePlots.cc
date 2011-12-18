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

float GetMuonTriggerWeight (float MuonPt, float MuonEta)
{
  if (MuonPt>=35 && MuonPt <40)
    {
      if (MuonEta < -1.5)
	return 0.856;
      if (MuonEta <= -1.0 && MuonEta > -1.5)
        return 0.879;
      if (MuonEta <= -0.5 && MuonEta > -1.0)
        return 0.913;
      if (MuonEta <= 0 && MuonEta > -0.5)
        return 0.932;
      if (MuonEta > 0 && MuonEta < 0.5)
        return 0.934;
      if (MuonEta >= 0.5 && MuonEta < 1.0)
        return 0.914;
      if (MuonEta >= 1.0 && MuonEta < 1.5)
        return 0.866;
      if (MuonEta > 1.5)
        return 0.865;
    }
  if (MuonPt>=40 && MuonPt <45)
    {
      if (MuonEta < -1.5)
        return 0.883;
      if (MuonEta <= -1.0 && MuonEta > -1.5)
        return 0.914;
      if (MuonEta <= -0.5 && MuonEta > -1.0)
        return 0.946;
      if (MuonEta <= 0 && MuonEta > -0.5)
        return 0.962;
      if (MuonEta > 0 && MuonEta < 0.5)
        return 0.962;
      if (MuonEta >= 0.5 && MuonEta < 1.0)
        return 0.945;
      if (MuonEta >= 1.0 && MuonEta < 1.5)
        return 0.898;
      if (MuonEta > 1.5)
        return 0.885;
    }
  if (MuonPt>=45 && MuonPt <50)
    {
      if (MuonEta < -1.5)
        return 0.888;
      if (MuonEta <= -1.0 && MuonEta > -1.5)
        return 0.913;
      if (MuonEta <= -0.5 && MuonEta > -1.0)
        return 0.946;
      if (MuonEta <= 0 && MuonEta > -0.5)
        return 0.963;
      if (MuonEta > 0 && MuonEta < 0.5)
        return 0.963;
      if (MuonEta >= 0.5 && MuonEta < 1.0)
        return 0.946;
      if (MuonEta >= 1.0 && MuonEta < 1.5)
        return 0.900;
      if (MuonEta > 1.5)
        return 0.890;
    }

  if (MuonPt>=50)
    {
      if (MuonEta < -1.5)
        return 0.889;
      if (MuonEta <= -1.0 && MuonEta > -1.5)
        return 0.814;
      if (MuonEta <= -0.5 && MuonEta > -1.0)
        return 0.945;
      if (MuonEta <= 0 && MuonEta > -0.5)
        return 0.962;
     if (MuonEta > 0 && MuonEta < 0.5)
        return 0.962;
      if (MuonEta >= 0.5 && MuonEta < 1.0)
        return 0.944;
      if (MuonEta >= 1.0 && MuonEta < 1.5)
        return 0.900;
      if (MuonEta > 1.5)
        return 0.887;
    }


  
  return 1.0;  
}






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

  int PtBinSize=10;
  int PtBinLow=0;
  int PtBinHigh=500;
  int PtNBins=(PtBinHigh-PtBinLow)/PtBinSize;

  int MassBinSize=60;
  int MassBinLow=0;
  int MassBinHigh=2000;
  int MassNBins=(MassBinHigh-MassBinLow)/MassBinSize;

  // Define histograms and set max number of jets to look at..
  int const NMaxJets = 7;
  TH1F* hJetPtMu[NMaxJets];
  TH1F* hJetEtaMu[NMaxJets];
  TH1F* hJetPhiMu[NMaxJets];
  TH1F* hJetPtEl[NMaxJets];
  TH1F* hJetEtaEl[NMaxJets];
  TH1F* hJetPhiEl[NMaxJets];
  TH1F* hJetPtMu_NoPile[NMaxJets];
  TH1F* hJetEtaMu_NoPile[NMaxJets];
  TH1F* hJetPtEl_NoPile[NMaxJets];
  TH1F* hJetEtaEl_NoPile[NMaxJets];

  for (int i = 0; i != NMaxJets; ++i) {
    hJetPtMu[i] = new TH1F( TString::Format("JetPt%iMu", i), TString::Format("JetPt%iMu", i), PtNBins, PtBinLow, PtBinHigh);
    hJetEtaMu[i] = new TH1F( TString::Format("JetEta%iMu", i), TString::Format("JetEta%iMu", i), 50, -4, 4);
    hJetPhiMu[i] = new TH1F( TString::Format("JetPhi%iMu", i), TString::Format("JetPhi%iMu", i), 50, -7, 7);
    hJetPtEl[i] = new TH1F( TString::Format("JetPt%iEl", i), TString::Format("JetPt%iEl", i), PtNBins, PtBinLow, PtBinHigh);

    hJetEtaEl[i] = new TH1F( TString::Format("JetEta%iEl", i), TString::Format("JetEta%iEl", i), 50, -4, 4);
    hJetPhiEl[i] = new TH1F( TString::Format("JetPhi%iEl", i), TString::Format("JetPhi%iEl", i), 50, -7, 7);
    hJetPtMu_NoPile[i] = new TH1F( TString::Format("JetPt%iMu_NoPile", i), TString::Format("JetPt%iMu_NoPile", i), PtNBins, PtBinLow, PtBinHigh);
    hJetEtaMu_NoPile[i] = new TH1F( TString::Format("JetEta%iMu_NoPile", i), TString::Format("JetEta%iMu_NoPile", i), 50, -4, 4);
    hJetPtEl_NoPile[i] = new TH1F( TString::Format("JetPt%iEl_NoPile", i), TString::Format("JetPt%iEl_NoPile", i), PtNBins, PtBinLow, PtBinHigh);
    hJetEtaEl_NoPile[i] = new TH1F( TString::Format("JetEta%iEl_NoPile", i), TString::Format("JetEta%iEl_NoPile", i), 50, -4, 4);

  }
  TH1F* LeptonPtMu = new TH1F("LeptonPtMu","LeptonPtMu", PtNBins, PtBinLow, PtBinHigh);
  TH1F* LeptonPtEl= new TH1F("LeptonPtEl","LeptonPtEl", PtNBins, PtBinLow, PtBinHigh);
  TH1F* LeptonPtMu_NoPile = new TH1F("LeptonPtMu_NoPile","LeptonPtMu_NoPile", PtNBins, PtBinLow, PtBinHigh);
  TH1F* LeptonPtEl_NoPile= new TH1F("LeptonPtEl_NoPile","LeptonPtEl_NoPile", PtNBins, PtBinLow, PtBinHigh);

  TH1F* LeptonPhiMu = new TH1F("LeptonPhiMu","LeptonPhiMu", 50, -7, 7);
  TH1F* LeptonPhiEl= new TH1F("LeptonPhiEl","LeptonPhiEl", 50, -7, 7);
  TH1F* LeptonEtaMu = new TH1F("LeptonEtaMu","LeptonEtaMu", 50, -4, 4);
  TH1F* LeptonEtaEl= new TH1F("LeptonEtaEl","LeptonEtaEl", 50, -4,4);

  TH1F* HTMu = new TH1F("HTMu","HTMu", 50, 0, 2000);
  TH1F* HTEl= new TH1F("HTEl","HTEl", 50, 0, 2000);
  TH1F* HTMu_NoPile = new TH1F("HTMu_NoPile","HTMu_NoPile", 50, 0, 2000);
  TH1F* HTEl_NoPile= new TH1F("HTEl_NoPile","HTEl_NoPile", 50, 0, 2000);

  TH1F* METMu = new TH1F("METMu","METMu",PtNBins, PtBinLow, PtBinHigh);
  TH1F* METEl= new TH1F("METEl","METEl", PtNBins, PtBinLow, PtBinHigh);
  TH1F* METMu_NoPile = new TH1F("METMu_NoPile","METMu_NoPile", PtNBins, PtBinLow, PtBinHigh);
  TH1F* METEl_NoPile= new TH1F("METEl_NoPile","METel_NoPile", PtNBins, PtBinLow, PtBinHigh);

  TH1F* HnPileUpVtx= new TH1F("nPileUpVtx","nPileUpVtx", PtNBins, PtBinLow, PtBinHigh);
  TH1F* HnPileUpVtx_NoPile = new TH1F("nPileUpVtx_NoPile","nPileUpVtx_NoPile", PtNBins, PtBinLow, PtBinHigh);

  TH1F*  hNumberOfJetsMu = new TH1F("NumberOfJetsMu","NumberOfJetsMu", 10, 5, 15);
  TH1F*  hNumberOfBJetsMu = new TH1F("NumberOfBJetsMu","NumberOfBJetsMu", 10, 0, 10);
  TH1F*  hNumberOfJetsEl = new TH1F("NumberOfJetsEl","NumberOfJetsEl", 10, 5, 15);
  TH1F*  hNumberOfBJetsEl = new TH1F("NumberOfBJetsEl","NumberOfBJetsEl", 10, 0, 10);

  TH1F* hBjetPtMu  = new TH1F("hBjetPtMu","hBjetPtMu", PtNBins, PtBinLow, PtBinHigh);
  TH1F* hBjetPtEl = new TH1F("hBjetPtEl","hBjetPtEl", PtNBins, PtBinLow, PtBinHigh);
  TH1F* hBjetEtaMu  = new TH1F("hBjetEtaMu","hBjetEtaMu", 50, -4, 4);
  TH1F* hBjetEtaEl = new TH1F("hBjetEtaEl","hBjetEtaEl", 50, -4, 4);
  // Define our chain of files/trees
  //Add the aysmmetr stuff
  TH1F* WprimeMassGoodMu = new TH1F("WprimeMassGoodMu","WprimeMassGoodMu", MassNBins, MassBinLow, MassBinHigh);
  TH1F* WprimeMassGoodEl = new TH1F("WprimeMassGoodEl","WprimeMassGoodEl", MassNBins, MassBinLow, MassBinHigh);
  TH1F* WprimeMassBadMu = new TH1F("WprimeMassBadMu","WprimeMassBadMu", MassNBins, MassBinLow, MassBinHigh);
  TH1F* WprimeMassBadEl = new TH1F("WprimeMassBadEl","WprimeMassBadEl", MassNBins, MassBinLow, MassBinHigh);
  TH1F* WprimeMassGoodLepton = new TH1F("WprimeMassGoodLepton","WprimeMassGoodLepton", MassNBins, MassBinLow, MassBinHigh);
  TH1F* WprimeMassBadLepton = new TH1F("WprimeMassBadLepton","WprimeMassBadLepton", MassNBins, MassBinLow, MassBinHigh);
  TChain Chain("micro", "micro");
  for (size_t i = 0; i != InFiles.size(); ++i) {
    Chain.Add(InFiles[i]);
  }

  // Init all branches
  initMicroNtuple((TTree*) &Chain);

  // Loop over all events in chain
  for (int ientry = 0; Chain.GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 1000 == 0) {
      //std::cout << "ientry: " << ientry << std::endl;
    }

    // Check for duplicate runs
    if(type==0){
    if (RunEvt.IsDuplicate(runNumber, eventNumber)) {
      std::cout << "Duplicate event!!  Skiping event: " << runNumber << "   " << eventNumber << std::endl;
      continue;
    }
    }
    //define weight without PileUp reweighting

    bool MuonEvent=false;
    bool ElectronEvent=false;
    if (muSEL > 1 and eSEL==0) MuonEvent=true;
    if (muSEL == 1 and eSEL>0) ElectronEvent=true;
    float BtagScaleFactorSSVHEM=1;
    float TriggerEff=1;
    //if(type>0)    BtagScaleFactorSSVHEM=0.95;
    //for muons get the trigger efficiency
    if (type>0){
      BtagScaleFactorSSVHEM=0.95;
      if (MuonEvent) TriggerEff=0.985*0.976;//GetMuonTriggerWeight(leptonPtRec,leptonEtaRec);
      if (ElectronEvent) TriggerEff=1;
    }

    // std::cout<<weight<<" "<<BtagScaleFactorSSVHEM<<" "<<TriggerEff<<std::endl;
      weight=weight*BtagScaleFactorSSVHEM*TriggerEff;

    float  weightNoPile=weight/PileUpWeight;
    //make lepton plots
    HnPileUpVtx->Fill(nPileUpVtx,weight);
    HnPileUpVtx_NoPile->Fill(nPileUpVtx,weightNoPile);
    TLorentzVector Jet0(JetPx[0], JetPy[0], JetPz[0], JetE[0]);
    TLorentzVector Jet1(JetPx[1], JetPy[1], JetPz[1], JetE[1]);
    if(HT>700.0 && Jet0.Pt()>180.0 && Jet1.Pt()>90.0){
      TLorentzVector BJet(BJetPx[0], BJetPy[0], BJetPz[0], BJetE[0]);
    if(muSEL > 1 && eSEL==0){

      LeptonPtMu->Fill(leptonPtRec,weight);
      LeptonPhiMu->Fill(leptonPhiRec,weight);
      LeptonEtaMu->Fill(leptonEtaRec,weight);
      LeptonPtMu_NoPile->Fill(leptonPtRec,weightNoPile);
      HTMu->Fill(HT,weight);
      HTMu_NoPile->Fill(HT,weightNoPile);
      METMu->Fill(MET,weight);
      METMu_NoPile->Fill(MET,weightNoPile);

      hNumberOfJetsMu->Fill(numberOfJets,weight);
      hNumberOfBJetsMu->Fill(numberOfBJets,weight);
      hBjetPtMu->Fill(BJet.Pt(), weight);
      hBjetEtaMu ->Fill(BJet.Eta(), weight);
      WprimeMassGoodMu->Fill(tNegJetMassRec,weight);
      WprimeMassBadMu->Fill(tPosJetMassRec,weight);


}
    if(eSEL > 1 && muSEL ==0){

      LeptonPtEl->Fill(leptonPtRec,weight);
      LeptonPhiEl->Fill(leptonPhiRec,weight);
      LeptonEtaEl->Fill(leptonEtaRec,weight);

      LeptonPtEl_NoPile->Fill(leptonPtRec,weightNoPile);
      HTEl->Fill(HT,weight);
      HTEl_NoPile->Fill(HT,weightNoPile);
      METEl->Fill(MET,weight);
      METEl_NoPile->Fill(MET,weightNoPile);
	hNumberOfJetsEl->Fill(numberOfJets,weight);
      hNumberOfBJetsEl->Fill(numberOfBJets,weight);
      hBjetPtEl->Fill(BJet.Pt(), weight);
      hBjetEtaEl ->Fill(BJet.Eta(), weight);
      WprimeMassGoodEl->Fill(tNegJetMassRec,weight);
      WprimeMassBadEl->Fill(tPosJetMassRec,weight);

    }
      WprimeMassGoodLepton->Fill(tNegJetMassRec,weight);
      WprimeMassBadLepton->Fill(tPosJetMassRec,weight);
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
        hJetPhiMu[ijet]->Fill(Jet.Phi(), weight);
	hJetPtMu_NoPile[ijet]->Fill(Jet.Pt(), weightNoPile);
        hJetEtaMu_NoPile[ijet]->Fill(Jet.Eta(), weightNoPile);
	}
	if(eSEL > 1 && muSEL ==0){
	  hJetPtEl[ijet]->Fill(Jet.Pt(), weight);
	  hJetEtaEl[ijet]->Fill(Jet.Eta(), weight);
	  hJetPhiEl[ijet]->Fill(Jet.Phi(), weight);
	  hJetPtEl_NoPile[ijet]->Fill(Jet.Pt(), weightNoPile);
          hJetEtaEl_NoPile[ijet]->Fill(Jet.Eta(), weightNoPile);


        }

        // This is just a sanity check...
      }
      if (Jet.Pt() > LastJetPt) {
        //std::cerr << "I don't know what's going on with the jets, but they need to be sorted" << std::endl;
      }
      LastJetPt = Jet.Pt();
    }

  }
  }
  OutFile.cd();
  for (int i = 0; i != NMaxJets; ++i) {
    hJetPtMu[i]->Write();
    hJetEtaMu[i]->Write();
    hJetPhiMu[i]->Write();
    hJetPtEl[i]->Write();
    hJetEtaEl[i]->Write();
    hJetPhiEl[i]->Write();
    //   hJetPtMu_NoPile[i]->Write();
    //hJetEtaMu_NoPile[i]->Write();
    //hJetPtEl_NoPile[i]->Write();
    //hJetEtaEl_NoPile[i]->Write();
  }
  LeptonPtMu->Write();
  LeptonPtEl->Write();
  LeptonPhiMu->Write();
  LeptonPhiEl->Write();
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
  WprimeMassGoodMu->Write();
  WprimeMassBadMu->Write();
  WprimeMassGoodEl->Write();
  WprimeMassBadEl->Write();
  WprimeMassGoodLepton->Write();
  WprimeMassBadLepton->Write();
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





