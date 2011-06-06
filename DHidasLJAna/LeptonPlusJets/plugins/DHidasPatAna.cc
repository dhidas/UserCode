
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h" // based on DataFormats/Candidate/interface/Particle.h
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/PatCandidates/interface/Photon.h" 


#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TTree.h"
#include "DHidasLJAna/LeptonPlusJets/plugins/DHidasPatAna.h"
#include "DHidasLJAna/LeptonPlusJets/interface/TAnaHist.h"

#include "TLorentzVector.h"





DHidasPatAna::DHidasPatAna(const edm::ParameterSet& iConfig)
{
  //outputFilename  = iConfig.getUntrackedParameter<string>("outputFilename","Pat_test.root");


  fIsData  = iConfig.getUntrackedParameter<bool>("IsData",false);
  fMakeDtuple  = iConfig.getUntrackedParameter<bool>("MakeDtuple",false);
  fOutFileName = iConfig.getUntrackedParameter<std::string>("OutFileName", "OutFile.root");
  fJSONFilename  = iConfig.getUntrackedParameter<std::string>("JSONFilename","");
  fTriggerNames = iConfig.getUntrackedParameter<std::vector<std::string> >("TriggerNames", std::vector<std::string>());

  for (std::vector<std::string>::iterator It = fTriggerNames.begin(); It != fTriggerNames.end(); ++It) {
    fTriggerMap[*It] = false;
  }

  // JSON file for data
  if (fIsData) {
    std::cout << "I see that This IS Data!!" << std::endl;
    std::cout << "Tryng JSON file: " << fJSONFilename << std::endl;
    fJSON.ReadFile(fJSONFilename);
  }


}


float DHidasPatAna::correct_met(float et, float oc, float phi, float eta, float metin, float metphiin, float met, float &cormetphi)
{

  if(fabs(eta)>1.45) return met;     // no overcleaning in EE
  float ratio = oc/et + 0.036*eta*eta;  // result of the fit of oc/et vs eta
  if(ratio > 0.85) return met;        
  if(acos(cos(phi-metphiin)) > 0.5) return met;
  float cet = et - (oc/(1.02-0.036*eta*eta));
  float metx = met*cos(metphiin) - cet*cos(phi);
  float mety = met*sin(metphiin) - cet*sin(phi);
  cormetphi = atan2(mety,metx);
  return sqrt(metx*metx+mety*mety);
}

//void DHidasPatAna::beginJob(const edm::EventSetup&) {
void DHidasPatAna::beginJob()
{
  // Open the output file
  fOutFile = new TFile(fOutFileName, "recreate");
  if (!fOutFile->IsOpen()) {
    std::cerr << "ERROR: cannot open output file: " << fOutFileName << std::endl;
    throw;
  }

  fTree = 0x0;
  if (fMakeDtuple) {
    fTree = new TTree("d", "Simple event tree");
    SetBranches(fTree);
    fTree->SetDirectory(fOutFile);
  }

  //////////////////
  ///  HISTOGRAMS
  //////////////////



}

bool DHidasPatAna::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (fIsData) {
    // Event info
    fRun         = iEvent.id().run();
    fLumiSection = iEvent.id().luminosityBlock();

    // If we're looking at data check for a good lumi section
    if (!fJSON.IsGoodLumiSection(fRun, fLumiSection)) {
      std::cout << "JSON Rejecting: ";
      std::cout << fRun << " " << fLumiSection << " " << fJSON.IsGoodLumiSection(fRun, fLumiSection) << std::endl;
      return false;
    }

    // Look for one of the triggers we care about
    bool HasTrigger = false;
    getTriggerDecision(iEvent, fTriggerMap);
    for (std::map<std::string, bool>::iterator It = fTriggerMap.begin(); It != fTriggerMap.end(); ++It) {
      if (It->second) {
        HasTrigger = true;
        break;
      }
    }
    if (!HasTrigger) {
      return false;
    }
  }

  // Look for a good lepton
  GetObjects(iEvent);
  if (fGoodElectrons.size() + fGoodMuons.size() > 0) {
    return true;
  }

  // no dice, return false.  don't take this event
  return false;
}

void DHidasPatAna::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Histograms for analyze
  static TAnaHist Hist(fOutFile, "analyze");

  // Event info
  fRun         = iEvent.id().run();
  fEvent       = iEvent.id().event();
  fLumiSection = iEvent.id().luminosityBlock();

  Hist.FillTH1F("RunNumber", 100, 160500, 161500, fRun);


  // If we're looking at data check for a good lumi section
  if (fIsData && !fJSON.IsGoodLumiSection(fRun, fLumiSection)) {
    return;
  }

  GetObjects(iEvent);

  /*
  std::vector<std::string> ElectronModule;
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronLWEt10PixelMatchFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt15CaloEleIdPixelMatchFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt17CaloEleIdPixelMatchFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt17TightEleIdDphiFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt22TighterEleIdDphiFilter");
  ElectronModule.push_back("hltL1NonIsoHLTNonIsoSingleElectronEt22TighterEleIdDphiFilter");
  ElectronModule.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
  ElectronModule.push_back("cleanElectronTriggerMatchHLTEle20SWL1R");
  matchElectrons(*TriggerEvent, ElectronModule, fCleanElectrons);
  */


  PlotObjects();
  PlotDileptonEvents();
  PlotMultiJetLeptonEvents();

  if (fTree) {
    FillTree();
  }

  // MET
  TVector2 const MET = MetColl->empty() ? TVector2(0,0) : TVector2((*MetColl)[0].px(), (*MetColl)[0].py());




  return;
}



void DHidasPatAna::endJob()
{
  if (fOutFile) {
    fOutFile->Write();
    fOutFile->Close();
  }
}



void DHidasPatAna::GetObjects (const edm::Event& iEvent)
{
  // Reset vectors
  fGoodElectrons.clear();
  fGoodPhotons.clear();
  fGoodMuons.clear();
  fGoodJets.clear();
  fCleanElectrons.clear();
  fCleanPhotons.clear();
  fCleanMuons.clear();
  fCleanJets.clear();

  // Get pat objects
  iEvent.getByLabel("selectedPatElectrons", PatElectrons); 
  //iEvent.getByLabel("generalTracks", recoTracks);
  iEvent.getByLabel("selectedPatPhotons", PatPhotons); 
  iEvent.getByLabel("selectedPatMuons", PatMuons); 
  iEvent.getByLabel("selectedPatJetsAK5PF", PatJets); 
  iEvent.getByLabel("patMETsPF",  MetColl);
  //iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
  iEvent.getByLabel("hltTriggerSummaryAOD", TriggerEvent);


  // Electron Selection
  for (size_t i = 0; i != PatElectrons->size(); ++i) {

    // Check electron ID
    int const eleid = (*PatElectrons)[i].electronID("simpleEleId80relIso");
    //        0: fails
    //        1: passes electron ID only
    //        2: passes electron Isolation only
    //        3: passes electron ID and Isolation only
    //        4: passes conversion rejection
    //        5: passes conversion rejection and ID
    //        6: passes conversion rejection and Isolation
    //        7: passes the whole selection

    if ( eleid == 7 && (*PatElectrons)[i].pt() > 30.0 && fabs((*PatElectrons)[i].eta())<2.1) {
      fGoodElectrons.push_back(&(*PatElectrons)[i]);
    }
  }


  // Photon Selection
  for (size_t i = 0; i != PatPhotons->size(); ++i) {
    if ((*PatPhotons)[i].et()>30.0 && fabs((*PatPhotons)[i].superCluster()->position().eta())<1.45) {
        // tight photons
        if ( ((*PatPhotons)[i].ecalRecHitSumEtConeDR04()< 4.2+0.006*(*PatPhotons)[i].et()) &&
            ((*PatPhotons)[i].hcalTowerSumEtConeDR04()< 2.2+0.0025*(*PatPhotons)[i].et()) && 
            ((*PatPhotons)[i].hadronicOverEm() < 0.05) &&
            ((*PatPhotons)[i].trkSumPtHollowConeDR04() < 2+0.001*(*PatPhotons)[i].et()) &&
            ((*PatPhotons)[i].sigmaIetaIeta() <0.013) &&
            (!((*PatPhotons)[i].hasPixelSeed()))) { 
          fGoodPhotons.push_back(&(*PatPhotons)[i]);
        }
      }
  }

  // Muon Selection
  for (size_t i = 0; i != PatMuons->size(); ++i) {

    double relIso = ((*PatMuons)[i].trackIso()  +
        (*PatMuons)[i].ecalIso()   +
        (*PatMuons)[i].hcalIso()) / (*PatMuons)[i].pt();

    int    nValidHits        = -1;
    int    nValidTrackerHits = -1;
    int    nValidPixelHits   = -1;

    if ((*PatMuons)[i].globalTrack().isNonnull()) {
      nValidHits = (*PatMuons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
    }

    if ((*PatMuons)[i].innerTrack().isNonnull()) {
      nValidTrackerHits = (*PatMuons)[i].innerTrack()->numberOfValidHits();
      nValidPixelHits   = (*PatMuons)[i].innerTrack()->hitPattern().pixelLayersWithMeasurement();
    }

    int stations = 0;
    unsigned stationMask((*PatMuons)[i].stationMask());
    for(unsigned j = 0; j < 8; ++j) {
      if(stationMask & 1 << j) {
        ++stations;
      }
    }

    if ((*PatMuons)[i].pt()>20.0 && fabs((*PatMuons)[i].eta())<2.1) {
      if((*PatMuons)[i].isGlobalMuon()  &&
         (*PatMuons)[i].isTrackerMuon() && 
         nValidHits                >  0 && 
         nValidTrackerHits         > 10 &&
         nValidPixelHits           >  0 &&
         (*PatMuons)[i].dB()       <  0.02 &&
         (*PatMuons)[i].globalTrack()->normalizedChi2() < 10 && 
         stations                > 1) {

        if( relIso  <  0.15) { // is good muon
          fGoodMuons.push_back(&(*PatMuons)[i]);
        }
      }
    }
  }

  // Jet Selection
  for (size_t i = 0; i != PatJets->size(); ++i) {
    if ((*PatJets)[i].pt()>30.0 && fabs((*PatJets)[i].eta())<2.6) {
      if 
        ((*PatJets)[i].correctedJet("Uncorrected").neutralHadronEnergyFraction()   < 0.99 && 
         (*PatJets)[i].correctedJet("Uncorrected").neutralEmEnergyFraction()       < 0.99 &&
         (*PatJets)[i].correctedJet("Uncorrected").numberOfDaughters()             > 1    &&
         (fabs((*PatJets)[i].eta())                    > 2.4  ||
          ((*PatJets)[i].correctedJet("Uncorrected").chargedHadronEnergyFraction() > 0.   &&
           (*PatJets)[i].correctedJet("Uncorrected").chargedEmEnergyFraction()     < 0.99 &&
           (*PatJets)[i].correctedJet("Uncorrected").chargedMultiplicity()         > 0.))) {
          fGoodJets.push_back(&(*PatJets)[i]);
        }




    }
  }


  // Sort out overlapping objects (Doesn't pat do this for you!?)
  for (size_t im = 0; im != fGoodMuons.size(); ++im) {
    fCleanMuons.push_back( fGoodMuons[im] );
  }

  // Keep non-overlapping electrons
  for (size_t ie = 0; ie != fGoodElectrons.size(); ++ie) {
    bool HasOverlap = false;
    TLorentzVector Electron(fGoodElectrons[ie]->px(), fGoodElectrons[ie]->py(), fGoodElectrons[ie]->pz(), fGoodElectrons[ie]->energy()); 
    for (size_t im = 0; im != fGoodMuons.size(); ++im) {
      TLorentzVector Muon(fGoodMuons[im]->px(), fGoodMuons[im]->py(), fGoodMuons[im]->pz(), fGoodMuons[im]->p()); 
      if (Muon.DeltaR( Electron ) < 0.4) {
        HasOverlap = true;
      }
    }
    if (!HasOverlap) {
      fCleanElectrons.push_back( fGoodElectrons[ie] );
    }
  }

  // Keep non-overlapping photons
  for (size_t ip = 0; ip != fGoodPhotons.size(); ++ip) {
    bool HasOverlap = false;
    TLorentzVector Photon(fGoodPhotons[ip]->px(), fGoodPhotons[ip]->py(), fGoodPhotons[ip]->pz(), fGoodPhotons[ip]->energy()); 
    for (size_t ie = 0; ie != fGoodElectrons.size(); ++ie) {
      TLorentzVector Electron(fGoodElectrons[ie]->px(), fGoodElectrons[ie]->py(), fGoodElectrons[ie]->pz(), fGoodElectrons[ie]->energy());
      if (Electron.DeltaR(Photon) < 0.4) {
        HasOverlap = true;
      }
    }
    if (!HasOverlap) {
      fCleanPhotons.push_back( fGoodPhotons[ip] );
    }
  }

  // Keep non-overlapping jets
  for (size_t ij = 0; ij != fGoodJets.size(); ++ij) {
    bool HasOverlap = false;
    TLorentzVector Jet(fGoodJets[ij]->px(), fGoodJets[ij]->py(), fGoodJets[ij]->pz(), fGoodJets[ij]->energy()); 
    for (size_t ie = 0; ie != fCleanElectrons.size(); ++ie) {
      TLorentzVector Electron(fCleanElectrons[ie]->px(), fCleanElectrons[ie]->py(), fCleanElectrons[ie]->pz(), fCleanElectrons[ie]->energy());
      if (Electron.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t ip = 0; ip != fCleanPhotons.size(); ++ip) {
      TLorentzVector Photon(fCleanPhotons[ip]->px(), fCleanPhotons[ip]->py(), fCleanPhotons[ip]->pz(), fCleanPhotons[ip]->energy());
      if (Photon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t im = 0; im != fCleanMuons.size(); ++im) {
      TLorentzVector Muon(fCleanMuons[im]->px(), fCleanMuons[im]->py(), fCleanMuons[im]->pz(), fCleanMuons[im]->energy());
      if (Muon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }

    if (!HasOverlap) {
      fCleanJets.push_back( fGoodJets[ij] );
    }
  }
  //printf("NElectrons=%i  NPhoton=%i  NMuons=%i  NJets=%i\n", (int) fCleanElectrons.size(), (int) fCleanPhotons.size(), (int) fCleanMuons.size(), (int) fCleanJets.size());
  return;
}





void DHidasPatAna::PlotObjects ()
{
  static TAnaHist Hist(fOutFile, "PlotObjects");

  Hist.FillTH1D("NElectrons", 20, 0, 20, fCleanElectrons.size());
  Hist.FillTH1D("NPhotons", 20, 0, 20, fCleanPhotons.size());
  Hist.FillTH1D("NMuons", 20, 0, 20, fCleanMuons.size());
  Hist.FillTH1D("NJets", 20, 0, 20, fCleanJets.size());

  for (size_t i = 0; i != fCleanElectrons.size(); ++i) {
    Hist.FillTH1D("ElectronEt", 30, 0, 300, fCleanElectrons[i]->et());
    Hist.FillTH1D("ElectronEta", 30, -3, 3, fCleanElectrons[i]->eta());
    Hist.FillTH1D("ElectronPhi", 30, -TMath::Pi(), TMath::Pi(), fCleanElectrons[i]->phi());
  }

  for (size_t i = 0; i != fCleanPhotons.size(); ++i) {
    Hist.FillTH1D("PhotonEt", 30, 0, 300, fCleanPhotons[i]->et());
    Hist.FillTH1D("PhotonEta", 30, -3, 3, fCleanPhotons[i]->eta());
    Hist.FillTH1D("PhotonPhi", 30, -TMath::Pi(), TMath::Pi(), fCleanPhotons[i]->phi());
  }

  for (size_t i = 0; i != fCleanMuons.size(); ++i) {
    Hist.FillTH1D("MuonEt", 30, 0, 300, fCleanMuons[i]->et());
    Hist.FillTH1D("MuonEta", 30, -3, 3, fCleanMuons[i]->eta());
    Hist.FillTH1D("MuonPhi", 30, -TMath::Pi(), TMath::Pi(), fCleanMuons[i]->phi());
  }

  for (size_t i = 0; i != fCleanJets.size(); ++i) {
    Hist.FillTH1D("JetEt", 30, 0, 300, fCleanJets[i]->et());
    Hist.FillTH1D("JetEta", 30, -3, 3, fCleanJets[i]->eta());
    Hist.FillTH1D("JetPhi", 30, -TMath::Pi(), TMath::Pi(), fCleanJets[i]->phi());
  }

  return;
}





void DHidasPatAna::PlotMultiJetLeptonEvents ()
{
  static TAnaHist Hist(fOutFile, "PlotMultiJetLeptonEvents");

  if (fCleanJets.size() < 3) {
    return;
  }
  if (fCleanElectrons.size() + fCleanMuons.size() != 1) {
    return;
  }

  TLorentzVector Lepton;
  TString LeptonType = "";
  if (fCleanElectrons.size() == 1) {
    Lepton.SetPxPyPzE(fCleanElectrons[0]->px(), fCleanElectrons[0]->py(), fCleanElectrons[0]->pz(), fCleanElectrons[0]->energy());
    LeptonType = "e";
  } else if (fCleanMuons.size() == 1) {
    Lepton.SetPxPyPzE(fCleanMuons[0]->px(), fCleanMuons[0]->py(), fCleanMuons[0]->pz(), fCleanMuons[0]->energy());
    LeptonType = "m";
  } else {
    return;
  }

  size_t const NJets = fCleanJets.size();
  std::vector<TLorentzVector> Jet(NJets);
  for (size_t i = 0; i != NJets; ++i) {
    Jet[i].SetPxPyPzE(fCleanJets[i]->px(), fCleanJets[i]->py(), fCleanJets[i]->pz(), fCleanJets[i]->energy());
  }
  for (size_t i = 0; i < NJets - 2; ++i) {
    for (size_t j = i+1; j < NJets - 1; ++j) {
      for (size_t k = j+1; k < NJets; ++k) {
        float const Mass      = (Jet[i]+Jet[j]+Jet[k]).M();
        float const SumPtJets = Jet[i].Pt() + Jet[j].Pt() + Jet[k].Pt();
        //printf("%i %i %i M=%7.2f  SumPtJets=%7.2f\n", (int) i, (int) j, (int) k, Mass, SumPtJets);
        Hist.FillTH2D("TrijetSumPt_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
        Hist.FillTH1D("TriJetMass", 100, 0, 800, Mass);
      }
    }
  }

  return;
}




void DHidasPatAna::PlotDileptonEvents ()
{
  static TAnaHist Hist(fOutFile, "PlotDileptonEvents");

  if ( fCleanElectrons.size() + fCleanMuons.size() < 2) {
    return;
  }

  TLorentzVector A, B;
  TString Type = "";
  if (fCleanElectrons.size() == 2 && fCleanElectrons[0]->charge() * fCleanElectrons[1]->charge() < 0) {
    A.SetPxPyPzE(fCleanElectrons[0]->px(), fCleanElectrons[0]->py(), fCleanElectrons[0]->pz(), fCleanElectrons[0]->energy());
    B.SetPxPyPzE(fCleanElectrons[1]->px(), fCleanElectrons[1]->py(), fCleanElectrons[1]->pz(), fCleanElectrons[1]->energy());
    Type = "ee";
  } else if (fCleanMuons.size() == 2 && fCleanMuons[0]->charge() * fCleanMuons[1]->charge() < 0) {
    A.SetPxPyPzE(fCleanMuons[0]->px(), fCleanMuons[0]->py(), fCleanMuons[0]->pz(), fCleanMuons[0]->energy());
    B.SetPxPyPzE(fCleanMuons[1]->px(), fCleanMuons[1]->py(), fCleanMuons[1]->pz(), fCleanMuons[1]->energy());
    Type = "mm";
  } else if (fCleanElectrons.size() == 1 && fCleanMuons.size() == 1 && fCleanElectrons[0]->charge() * fCleanMuons[0]->charge() < 0) {
    A.SetPxPyPzE(fCleanElectrons[0]->px(), fCleanElectrons[0]->py(), fCleanElectrons[0]->pz(), fCleanElectrons[0]->energy());
    B.SetPxPyPzE(fCleanMuons[0]->px(), fCleanMuons[0]->py(), fCleanMuons[0]->pz(), fCleanMuons[0]->energy());
    Type = "em";
  } else {
    return;
  }

  // Order by Pt
  if (A.Pt() < B.Pt()) {
    TLorentzVector C(A);
    A = B;
    B = C;
  }

  Hist.FillTH1D("M"+Type, 30, 0, 300, (A+B).M());
  Hist.FillTH1D("Pt0_"+Type, 30, 0, 300, A.Pt());
  Hist.FillTH1D("Pt1_"+Type, 30, 0, 300, B.Pt());

  TString Name = "M" + Type + "_";
  Name += (int) fCleanJets.size();
  Name += "j";
  Hist.FillTH1D(Name, 30, 0, 300, (A+B).M());

  return;
}




void DHidasPatAna::getTriggerDecision(const edm::Event& iEvent, std::map<std::string, bool>& TriggerMap)
{
  edm::Handle<edm::TriggerResults> triggerResults;

  std::string menu = "HLT";
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", menu), triggerResults);

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(* triggerResults);

  for (std::map<std::string, bool>::iterator It = TriggerMap.begin(); It != TriggerMap.end(); ++It) {
    It->second = false;
    unsigned int triggerIndex = triggerNames.triggerIndex(It->first);

    if (triggerIndex < triggerResults->size()) {
      if (triggerResults->accept(triggerIndex)) {
        It->second = true;
      }
    }

  }
  return;
}




void DHidasPatAna::FillTree ()
{
  ClearDtuple();

  fEvt.Run = fRun;
  fEvt.Event = fEvent;
  fEvt.LumiSection = fLumiSection;

  fEvt.NLeptons = 0;
  for (size_t i = 0; i != fCleanMuons.size(); ++i) {
    fEvt.LeptonPx->push_back(fCleanMuons[i]->px());
    fEvt.LeptonPy->push_back(fCleanMuons[i]->py());
    fEvt.LeptonPz->push_back(fCleanMuons[i]->pz());
    fEvt.LeptonPt->push_back(fCleanMuons[i]->pt());
    fEvt.LeptonType->push_back(Dtuple::kMuon);
    ++fEvt.NLeptons;
  }

  for (size_t i = 0; i != fCleanElectrons.size(); ++i) {
    fEvt.LeptonPx->push_back(fCleanElectrons[i]->px());
    fEvt.LeptonPy->push_back(fCleanElectrons[i]->py());
    fEvt.LeptonPz->push_back(fCleanElectrons[i]->pz());
    fEvt.LeptonPt->push_back(fCleanElectrons[i]->pt());
    fEvt.LeptonType->push_back(Dtuple::kElectron);
    ++fEvt.NLeptons;
  }

  fEvt.NPhotons = 0;
  for (size_t i = 0; i != fCleanPhotons.size(); ++i) {
    fEvt.PhotonPx->push_back(fCleanPhotons[i]->px());
    fEvt.PhotonPy->push_back(fCleanPhotons[i]->py());
    fEvt.PhotonPz->push_back(fCleanPhotons[i]->pz());
    fEvt.PhotonPt->push_back(fCleanPhotons[i]->pt());
    ++fEvt.NPhotons;
  }


  fEvt.NJets = 0;
  fEvt.SumPtJets = 0;
  for (size_t i = 0; i != fCleanJets.size(); ++i) {
    fEvt.JetPx->push_back(fCleanJets[i]->px());
    fEvt.JetPy->push_back(fCleanJets[i]->py());
    fEvt.JetPz->push_back(fCleanJets[i]->pz());
    fEvt.JetPt->push_back(fCleanJets[i]->pt());
    fEvt.SumPtJets += fCleanJets[i]->pt();
    ++fEvt.NJets;
  }


  size_t const NJets = fCleanJets.size();
  if (NJets >= 3) {
    std::vector<TLorentzVector> Jet(NJets);
    for (size_t i = 0; i != NJets; ++i) {
      Jet[i].SetPxPyPzE(fCleanJets[i]->px(), fCleanJets[i]->py(), fCleanJets[i]->py(), fCleanJets[i]->energy());
    }
    for (size_t i = 0; i < NJets - 2; ++i) {
      for (size_t j = i+1; j < NJets - 1; ++j) {
        for (size_t k = j+1; k < NJets; ++k) {
          fEvt.TriJetMasses->push_back( (Jet[i]+Jet[j]+Jet[k]).M() );
          fEvt.TriJetSumPt->push_back( Jet[i].Pt() + Jet[j].Pt() + Jet[k].Pt() );
        }
      }
    }
  }

  TVector2 const MET = MetColl->empty() ? TVector2(0,0) : TVector2((*MetColl)[0].px(), (*MetColl)[0].py());
  fEvt.METMag = MET.Mod();
  fEvt.METPhi = MET.Phi();

  fTree->Fill();

  return;
}











DEFINE_FWK_MODULE(DHidasPatAna); // define this as a plug-in
