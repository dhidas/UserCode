// -*- C++ -*-
//
// Package:    FillDtuple
// Class:      FillDtuple
// 
/**\class FillDtuple FillDtuple.cc DHidasAna/FillDtuple/src/FillDtuple.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dean Andrew HIDAS
//         Created:  Mon Oct 26 11:59:20 CET 2009
// $Id: FillDtuple.cc,v 1.24 2010/03/17 16:50:40 dhidas Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DHidasAna/Dtuple/interface/TDtuple.h"


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"


#include "TString.h"
#include "TH1D.h"
#include <string>
#include <vector>




#include <algorithm>


//
// class decleration
//

class FillDtuple : public edm::EDAnalyzer {
  public:
    explicit FillDtuple(const edm::ParameterSet&);
    ~FillDtuple();


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void GetHandles (const edm::Event&);
    void FillBasicEventQuantities (const edm::Event&);
    void FillLeptons (const edm::Event&);
    void FillPhotons (const edm::Event&);
    void FillJets (const edm::Event&);
    void DoMCParticleMatching (const edm::Event&);

    void PrintGenP (const edm::Event&, size_t const NMax = 999999999);

    void EventSummary ();
    bool KeepEvent ();

    TDtuple* fDtuple;

    edm::Handle< edm::View<pat::Electron> > fElectrons;
    edm::Handle< edm::View<pat::Muon> > fMuons;
    edm::Handle< edm::View<pat::Jet> > fJets;
    edm::Handle< edm::View<pat::Photon> > fPhotons;
    edm::Handle< edm::View<pat::MET> > fMETs;
    edm::Handle<reco::BeamSpot> fBeamSpot;
    edm::Handle<reco::TrackCollection> fTrackCollection;
    edm::Handle<reco::GenParticleCollection> fGenParticleCollection;

    edm::Handle< pat::TriggerEvent > fTriggerEvent;

    // ----------member data ---------------------------

    // ---------- Histograms ---------------------------
    std::map<TString, TH1D*> Hist1D;
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FillDtuple::FillDtuple(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


FillDtuple::~FillDtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  // Bye bye writer
  delete fDtuple;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
FillDtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Better get all the labels you want right here!
  GetHandles(iEvent);


  // Set Dtuple to default values
  fDtuple->DefaultValues();


  // Load this event in the dtuple with stuff of interest
  FillBasicEventQuantities(iEvent);
  FillLeptons(iEvent);
  FillPhotons(iEvent);
  FillJets(iEvent);

  DoMCParticleMatching(iEvent);
  //PrintGenP(iEvent);

  // Print the even summary
  //EventSummary();

  // Actually save this event!!
  if (KeepEvent()) {
    fDtuple->Fill();
  }


  return;
}




// ------------ Get Handles  ------------
void 
FillDtuple::GetHandles(const edm::Event& iEvent)
{
  // This function gets the physics objects of interest
  iEvent.getByLabel("selectedLayer1Electrons", fElectrons);
  iEvent.getByLabel("selectedLayer1Muons", fMuons);
  iEvent.getByLabel("selectedLayer1Jets", fJets);
  iEvent.getByLabel("selectedLayer1Photons", fPhotons);
  iEvent.getByLabel("layer1METs", fMETs);
  iEvent.getByLabel("offlineBeamSpot", fBeamSpot);
  iEvent.getByLabel("generalTracks", fTrackCollection);
  iEvent.getByLabel("genParticles", fGenParticleCollection);
  //iEvent.getByLabel("patTrigger", fTriggerEvent );
  //iEvent.getByLabel( "patTrigger", fTriggerPaths );
  //iEvent.getByLabel( "patTrigger", fTriggerFilters );
  //iEvent.getByLabel( "patTrigger", fTriggerObjects );
  return;
}


// ------------ Fill basic event quantities  ------------
void 
FillDtuple::FillBasicEventQuantities (const edm::Event& iEvent)
{
  // Fill some event quantities

  fDtuple->SetRun    ( iEvent.id().run() );
  fDtuple->SetEvent  ( iEvent.id().event() );
  fDtuple->SetMetX   ( fMETs->front().px() );
  fDtuple->SetMetY   ( fMETs->front().py() );
  fDtuple->SetSumEt  ( fMETs->front().sumEt() );
  //fDtuple->SetMetSig ( fMETs->front().mEtSig() );

  return;
}


// ------------ Fill the electrons  ------------
void 
FillDtuple::FillLeptons(const edm::Event& iEvent)
{

  // The leptons we find here!
  std::vector<TLepton> Leptons;

  for (size_t i = 0; i != fElectrons->size(); ++i) {

    // This lepton
    TLepton Lep;

    // This electron
    pat::Electron electron = fElectrons->at(i);
    //if (electron.electronID("eidRobustHighEnergy") != 1) {
    //  continue;
    //}
    //std::cout << "Electron Pt: " << electron.et() << std::endl;

    Lep.SetPx( electron.px() );
    Lep.SetPy( electron.py() );
    Lep.SetPz( electron.pz() );
    Lep.SetE( electron.p() );
    Lep.SetTrkPt( sqrt(electron.trackMomentumAtVtx().perp2()) );
    Lep.Setdxy( electron.gsfTrack()->dxy(fBeamSpot->position()) );
    Lep.Setdz( electron.gsfTrack()->dz(fBeamSpot->position()) );

    Lep.SetTrackChi2( electron.gsfTrack()->chi2() );
    Lep.SetTrackNDoF( electron.gsfTrack()->ndof() );
    Lep.SetNValidTrackerHits( electron.gsfTrack()->hitPattern().numberOfValidTrackerHits() );
    Lep.SetZ0( electron.vz() );
    Lep.SetCharge( electron.charge() );
    Lep.SetFlavor( TLepton::kLeptonFlavor_Electron );
    //Lep.SetTrkIso( electron.trackIso() );
    Lep.SetTrkIso( electron.dr03TkSumPt() );
    Lep.SetCalIso( electron.caloIso() );
    //Lep.SetECalIso( electron.ecalIso() );
    Lep.SetECalIso( electron.dr03EcalRecHitSumEt() );
    //Lep.SetHCalIso( electron.hcalIso() );
    Lep.SetHCalIso( electron.dr03HcalTowerSumEt() );
    Lep.SetECalIsoDep( electron.ecalIsoDeposit()->candEnergy() );
    Lep.SetHCalIsoDep( electron.hcalIsoDeposit()->candEnergy() );
    Lep.SetCalE( electron.caloEnergy() );
    Lep.SetHCalOverECal( electron.hadronicOverEm() );
    Lep.SetEoverPin( electron.eSuperClusterOverP() );
    Lep.SetfBrem( electron.fbrem() );
    //Lep.SetIsConvertedPhoton( electron.isConvertedPhoton() );
    Lep.SetIsConvertedPhoton( ConversionFinder::isElFromConversion( (reco::GsfElectron) electron, fTrackCollection, 3.8) ? 1 : 0 );
    int PassSelection = 0x0;
    if (electron.electronID("eidRobustHighEnergy") == 1) {
      PassSelection |= (0x1 << TLepton::kElectronSel_RobustHighEnergy);
    }
    if (electron.electronID("eidRobustLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kElectronSel_RobustLoose);
    }
    if (electron.electronID("eidRobustTight") == 1) {
      PassSelection |= (0x1 << TLepton::kElectronSel_RobustTight);
    }
    if (electron.electronID("eidLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kElectronSel_Loose);
    }
    if (electron.electronID("eidTight") == 1) {
      PassSelection |= (0x1 << TLepton::kElectronSel_Tight);
    }
    Lep.SetPassSelection(PassSelection);

    // Electron regions
    int Detector = 0x0;
    if (electron.isEE())        { Detector |= (0x1 << TLepton::kElectronDet_EE); }
    if (electron.isEB())        { Detector |= (0x1 << TLepton::kElectronDet_EB); }
    if (electron.isEBEEGap())   { Detector |= (0x1 << TLepton::kElectronDet_EBEEGap); }
    if (electron.isEBEtaGap())  { Detector |= (0x1 << TLepton::kElectronDet_EBEtaGap); }
    if (electron.isEBGap())     { Detector |= (0x1 << TLepton::kElectronDet_EBGap); }
    if (electron.isEBPhiGap())  { Detector |= (0x1 << TLepton::kElectronDet_EBPhiGap); }
    if (electron.isEEDeeGap())  { Detector |= (0x1 << TLepton::kElectronDet_EEDeeGap); }
    if (electron.isEEGap())     { Detector |= (0x1 << TLepton::kElectronDet_EEGap); }
    if (electron.isEERingGap()) { Detector |= (0x1 << TLepton::kElectronDet_EERingGap); }
    Lep.SetDetector(Detector);

    // Electron Classification
    switch (electron.classification()) {
      case reco::GsfElectron::UNKNOWN:
        Lep.SetClassification( TLepton::kElectronClass_Unknown );
        break;
      case reco::GsfElectron::GOLDEN:
        Lep.SetClassification( TLepton::kElectronClass_Golden );
        break;
      case reco::GsfElectron::BIGBREM:
        Lep.SetClassification( TLepton::kElectronClass_BigBrem );
        break;
      case reco::GsfElectron::NARROW:
        Lep.SetClassification( TLepton::kElectronClass_Narrow );
        break;
      case reco::GsfElectron::SHOWERING:
        Lep.SetClassification( TLepton::kElectronClass_Showering );
        break;
      case reco::GsfElectron::GAP:
        Lep.SetClassification( TLepton::kElectronClass_Gap );
        break;
      default:
        std::cerr << "ERROR in classifying lepton.  classification(): " << electron.classification() << std::endl;
        break;
    }

    Lep.SetSigmaIEtaIEta( electron.scSigmaIEtaIEta() );
    Lep.SetDeltaEtaIn( electron.deltaEtaSuperClusterTrackAtVtx() );
    Lep.SetDeltaPhiIn( electron.deltaPhiSuperClusterTrackAtVtx() );
    Lep.SetE2x5overE5x5( electron.scE2x5Max() / electron.scE5x5() );

    // This is to study conversions.

    // Get the reco track which is the closest match to the gsf track
    const reco::Track* electronCTFTrack = ConversionFinder::getElectronTrack(electron, 0.45);
    if (electronCTFTrack) {
      // Get the next closest track to the gsf electron
      reco::TrackRef Partner  = ConversionFinder::getConversionPartnerTrack((reco::GsfElectron) electron, fTrackCollection, 3.8, 99999, 99999);
      if (Partner.isNonnull()) {
        math::XYZTLorentzVector electron4V(electronCTFTrack->px(), electronCTFTrack->py(), electronCTFTrack->pz(), electronCTFTrack->p());
        math::XYZTLorentzVector Partnet4V(Partner->px(), Partner->py(), Partner->pz(), Partner->p());

        // Calculate the distance and delta-cot-theta
        std::pair<double, double> convInfo =  ConversionFinder::getConversionInfo(electron4V, electronCTFTrack->charge(), electronCTFTrack->d0(),
            Partnet4V, Partner->charge(), Partner->d0(), 3.8);

        // Save the Dist and dCotTheta to dtuple
        Lep.SetConvDist( convInfo.first );
        Lep.SetConvdCotTheta( convInfo.second );
      }
    }
    // end conversion test

    Leptons.push_back(Lep);
  }




  for (size_t i = 0; i != fMuons->size(); ++i) {

    // This lepton
    TLepton Lep;

    // This muon
    pat::Muon muon = fMuons->at(i);

    Lep.SetPx( muon.px() );
    Lep.SetPy( muon.py() );
    Lep.SetPz( muon.pz() );
    Lep.SetE( muon.p() );
    Lep.SetTrkPt( muon.pt() );
    reco::TrackRef MyTrackRef = muon.innerTrack();
    if (MyTrackRef.isNonnull()) {
      Lep.Setdxy( MyTrackRef->dxy(fBeamSpot->position()) );
      Lep.Setdz( MyTrackRef->dz(fBeamSpot->position()) );
      Lep.SetTrackChi2( MyTrackRef->chi2() );
      Lep.SetTrackNDoF( MyTrackRef->ndof() );
    }
    //Lep.SetNValidTrackerHits( MyTrackRef->hitPattern().numberOfValidTrackerHits() );
    Lep.SetZ0( muon.vz() );
    Lep.SetCharge( muon.charge() );
    Lep.SetFlavor( TLepton::kLeptonFlavor_Muon );
    Lep.SetTrkIso( muon.trackIso() );
    Lep.SetCalIso( muon.caloIso() );
    Lep.SetECalIsoDep( muon.ecalIsoDeposit()->candEnergy() );
    Lep.SetHCalIsoDep( muon.hcalIsoDeposit()->candEnergy() );
    Lep.SetECalIso( muon.ecalIso() );
    Lep.SetHCalIso( muon.hcalIso() );
    Lep.SetCalE( muon.calEnergy().em + muon.calEnergy().had );
    Lep.SetHCalOverECal( muon.calEnergy().had / muon.calEnergy().em );

    // Muon selections
    int PassSelection = 0x0;
    if (muon.muonID("TrackerMuonArbitrated") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TrackerMuonArbitrated);
    }
    if (muon.muonID("AllArbitrated") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_AllArbitrated);
    }
    if (muon.muonID("GlobalMuonPromptTight") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_GlobalMuonPromptTight);
    }
    if (muon.muonID("TMLastStationLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMLastStationLoose);
    }
    if (muon.muonID("TMLastStationTight") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMLastStationTight);
    }
    if (muon.muonID("TM2DCompatibilityLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TM2DCompatibilityLoose);
    }
    if (muon.muonID("TM2DCompatibilityTight") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TM2DCompatibilityTight);
    }
    if (muon.muonID("TMOneStationLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMOneStationLoose);
    }
    if (muon.muonID("TMOneStationTight") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMOneStationTight);
    }
    if (muon.muonID("TMLastStationOptimizedLowPtLoose") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMLastStationOptimizedLowPtLoose);
    }
    if (muon.muonID("TMLastStationOptimizedLowPtTight") == 1) {
      PassSelection |= (0x1 << TLepton::kMuonSel_TMLastStationOptimizedLowPtTight);
    }
    Lep.SetPassSelection( PassSelection );

    // Muon types
    int Detector = 0x0;
    if (muon.isGlobalMuon())     { Detector |= (0x1 << TLepton::kMuonDet_Global); }
    if (muon.isTrackerMuon())    { Detector |= (0x1 << TLepton::kMuonDet_Tracker); }
    if (muon.isStandAloneMuon()) { Detector |= (0x1 << TLepton::kMuonDet_StandAlone); }
    if (muon.isCaloMuon())       { Detector |= (0x1 << TLepton::kMuonDet_Calo); }
    Lep.SetDetector( Detector );

    Leptons.push_back(Lep);
  }


  fDtuple->SortLeptons(Leptons.begin(), Leptons.end());
  fDtuple->AddLeptons(Leptons);




  return;
}


// ------------ Fill the jets  ------------
void 
FillDtuple::FillJets(const edm::Event& iEvent)
{

  std::vector<TJet> Jets;
  for (size_t i = 0; i != fJets->size(); ++i) {

    // This jet
    TJet MyJet;

    pat::Jet jet = fJets->at(i);

    MyJet.SetPx( jet.px() );
    MyJet.SetPy( jet.py() );
    MyJet.SetPz( jet.pz() );
    MyJet.SetE( jet.p() );
    //MyJet.SetPt( jet.pt() );
    //MyJet.SetEta( jet.eta() );
    //MyJet.SetPhi( jet.phi() );
    MyJet.SetEmF( jet.emEnergyFraction() );
    MyJet.SetHadF( jet.energyFractionHadronic() );

    Jets.push_back(MyJet);

  }

  fDtuple->SortJets(Jets.begin(), Jets.end());
  fDtuple->AddJets(Jets);

  return;
}



// ------------ Fill the photons ------------
void 
FillDtuple::FillPhotons (const edm::Event& iEvent)
{
  std::vector<TPhoton> Photons;

  for (size_t i = 0; i != fPhotons->size(); ++i) {

    TPhoton MyPhoton;

    pat::Photon photon = fPhotons->at(i);

    MyPhoton.SetPx( photon.px() );
    MyPhoton.SetPy( photon.py() );
    MyPhoton.SetPz( photon.pz() );
    MyPhoton.SetE(  photon.p() );
    MyPhoton.SetTrkIso( photon.trackIso() );
    MyPhoton.SetCalIso( photon.caloIso() );
    MyPhoton.SetHCalOverECal( photon.hadronicOverEm() );

    Photons.push_back(MyPhoton);
  }

  fDtuple->SortPhotons(Photons.begin(), Photons.end());
  fDtuple->AddPhotons(Photons);

  return;
}




void FillDtuple::DoMCParticleMatching(const edm::Event& iEvent)
{
  //std::cout << "Number of GenParticles: " << fGenParticleCollection->size() << std::endl;
  for(size_t iGen = 0; iGen < fGenParticleCollection->size(); ++iGen) {
    reco::GenParticle const GenP = fGenParticleCollection->at(iGen);
    int const Status = GenP.status();
    if (Status != 1) {
      continue;
    }
    if (GenP.pt() == 0) {
      continue;
    }
    int const Id = GenP.pdgId();
    //int const MotherId = GenP.mother() != 0x0 ? GenP.mother()->pdgId() : 0;
    TGenP ThisGenP;
    ThisGenP.SetPxPyPzE(GenP.px(), GenP.py(), GenP.pz(), GenP.p());
    ThisGenP.SetId(Id);

    reco::GenParticle* AMother = (reco::GenParticle*) GenP.mother();
    while (AMother && AMother->pdgId() == GenP.pdgId()) {
      AMother = (reco::GenParticle*) AMother->mother();
    }
    //if (Id == AMother->pdgId()) {
    //  std::cout << "Me - Mother.Pt():  " << GenP.pt() - AMother->pt() << std::endl;
    //  std::cout << "Status MoStatus:   " << GenP.status() << "   " <<  AMother->status() << std::endl;
    //  continue;
    //}
    //if (AMother == 0x0) {
    //  continue;
    //}
    int const MotherId = AMother != 0x0 ? AMother->pdgId() : 0;
    ThisGenP.SetMotherId(MotherId);
    //std::cout << "GenPid Motherid  " << Id << "  " << AMother->pdgId() << std::endl;

    // This needs to be worked on.  consider what you might match here...
    float MinDeltaR = 99999;
    for (std::vector<TLepton>::iterator Lep = fDtuple->GetLeptons()->begin(); Lep != fDtuple->GetLeptons()->end(); ++Lep) {
      if ( TMath::Abs(Lep->DeltaR(ThisGenP)) < 0.4 ) {

        Lep->GenP.push_back(ThisGenP);
        //printf("LeptonMatch Id MotherId Status DeltaR Pt GenPt MoGenPt: %7i %7i %7i %9.3f %9.3f %9.3f %9.3f\n",
        //    Id,
        //    MotherId,
        //    Status,
        //    Lep->DeltaR(ThisGenP),
        //    Lep->Perp(),
        //    ThisGenP.Perp(),
        //    AMother->pt());

      }
    }

  }

  return;
}





void FillDtuple::PrintGenP(const edm::Event& iEvent, size_t const NMax)
{
  std::cout << "Number of GenParticles: " << fGenParticleCollection->size() << std::endl;
  for(size_t iGen = 0; iGen < fGenParticleCollection->size() && iGen < NMax; ++iGen) {
    reco::GenParticle const GenP = fGenParticleCollection->at(iGen);

    reco::GenParticle* Mother = (reco::GenParticle*) GenP.mother();
    if (Mother) {
      printf("%6i %6i %3i %8.2f %6i %3i %8.2f\n",
          iGen,
          GenP.pdgId(),
          GenP.status(),
          GenP.pt(),
          Mother->pdgId(),
          Mother->status(),
          Mother->pt());
    }
  }

  return;
}



// ------------ print brief event summary ------------
void 
FillDtuple::EventSummary ()
{
  // Print a very basic event summary
  printf("Leptons: %5i Photons: %5i Jets: %5i Met: %8.1f\n",
      (int) fDtuple->GetLeptons()->size(),
      (int) fDtuple->GetPhotons()->size(),
      (int) fDtuple->GetJets()->size(),
      fDtuple->GetMet());
  return;
}




// ------------ Keep this event or not? ------------
bool
FillDtuple::KeepEvent ()
{
  bool keep = true;

  return keep;
}




// ------------ method called once each job just before starting event loop  ------------
void 
FillDtuple::beginJob()
{
  edm::Service<TFileService> fs;
  //TTree* DtupleTree = fs->make<TTree>("Dtuple", "Dtuple");
  //DtupleTree->SetDirectory( &(fs->file()) );

  fDtuple = new TDtuple( &(fs->file()) );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
FillDtuple::endJob() {
  //delete fDtuple;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FillDtuple);
