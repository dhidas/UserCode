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
// $Id: FillDtuple.cc,v 1.13 2010/02/02 16:18:39 dhidas Exp $
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


#include "DHidasAna/Dtuple/interface/DtupleWriter.h"


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include "TString.h"
#include <string>




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
    void FillBasicEventQuantities (const edm::Event&, DtupleWriter::Event_Struct&);
    void FillLeptons (const edm::Event&, DtupleWriter::Event_Struct&);
    void FillPhotons (const edm::Event&, DtupleWriter::Event_Struct&);
    void FillJets (const edm::Event&, DtupleWriter::Event_Struct&);

    void EventSummary (DtupleWriter::Event_Struct&);
    bool KeepEvent (DtupleWriter::Event_Struct&);

    Dtuple* fDtupleWriter;

    edm::Handle< edm::View<pat::Electron> > fElectrons;
    edm::Handle< edm::View<pat::Muon> > fMuons;
    edm::Handle< edm::View<pat::Jet> > fJets;
    edm::Handle< edm::View<pat::Photon> > fPhotons;
    edm::Handle< edm::View<pat::MET> > fMETs;
    edm::Handle<reco::BeamSpot> fBeamSpot;

    edm::Handle< pat::TriggerEvent > fTriggerEvent;
    //edm::Handle< pat::TriggerPathCollection > fTriggerPaths;
    //edm::Handle< pat::TriggerFilterCollection > fTriggerFilters;
    //edm::Handle< pat::TriggerObjectCollection > fTriggerObjects;

    struct ReallyBasicLepton {
      float Pt;
      int   Flavor;
      int   Id;
    };
    static bool CompareReallyBasicLeptonPt (FillDtuple::ReallyBasicLepton A, FillDtuple::ReallyBasicLepton B)
    {
      return (A.Pt > B.Pt);
    };

    // ----------member data ---------------------------
    std::string fOutFileName;
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
FillDtuple::FillDtuple(const edm::ParameterSet& iConfig) :
fOutFileName( iConfig.getUntrackedParameter< std::string >( "OutFileName" ) )
{
   //now do what ever initialization is needed

}


FillDtuple::~FillDtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  // Bye bye writer
  delete fDtupleWriter;

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
  fDtupleWriter->DefaultValues();

  static DtupleWriter::Event_Struct& Ev = fDtupleWriter->fEvent;


  // Load this event in the dtuple with stuff of interest
  FillBasicEventQuantities(iEvent, Ev);
  FillLeptons(iEvent, Ev);
  FillPhotons(iEvent, Ev);
  FillJets(iEvent, Ev);

  // Print the even summary
  EventSummary(Ev);

  // Actually save this event!!
  if (KeepEvent(Ev)) {
    fDtupleWriter->Fill();
  }


  return;
}




// ------------ Get Handles  ------------
void 
FillDtuple::GetHandles(const edm::Event& iEvent)
{
  // This function gets the physics objects of interest
  iEvent.getByLabel("allLayer1Electrons", fElectrons);
  iEvent.getByLabel("allLayer1Muons", fMuons);
  iEvent.getByLabel("allLayer1Jets", fJets);
  iEvent.getByLabel("allLayer1Photons", fPhotons);
  iEvent.getByLabel("layer1METs", fMETs);
  iEvent.getByLabel("offlineBeamSpot", fBeamSpot);
  //iEvent.getByLabel("patTrigger", fTriggerEvent );
  //iEvent.getByLabel( "patTrigger", fTriggerPaths );
  //iEvent.getByLabel( "patTrigger", fTriggerFilters );
  //iEvent.getByLabel( "patTrigger", fTriggerObjects );
  return;
}


// ------------ Fill basic event quantities  ------------
void 
FillDtuple::FillBasicEventQuantities (const edm::Event& iEvent, DtupleWriter::Event_Struct& Ev)
{
  // Fill some event quantities

  Ev.Run    = iEvent.id().run();
  Ev.Event  = iEvent.id().event();
  Ev.MetMag = fMETs->front().et();
  Ev.MetPhi = fMETs->front().phi();
  Ev.SumEt  = fMETs->front().sumEt();
  Ev.MetSig = fMETs->front().mEtSig();

  return;
}


// ------------ Fill the electrons  ------------
void 
FillDtuple::FillLeptons(const edm::Event& iEvent, DtupleWriter::Event_Struct& Ev)
{
  std::vector<ReallyBasicLepton> Leptons;
  static ReallyBasicLepton ThisLep;

  // PAT object collection
  //edm::Handle< pat::MuonCollection > muons;
  //iEvent.getByLabel( "selectedLayer1Muons", muons );
  //const pat::helper::TriggerMatchHelper matchHelper;
  //const pat::TriggerObjectMatch* triggerMatch( fTriggerEvent->triggerObjectMatchResult( "muonTriggerMatchHLTMuons" ) );
  //for (size_t i = 0; i != muons->size(); ++i) {
  //  const reco::CandidateBaseRef candBaseRef( pat::MuonRef( muons, i ) );
  //  const pat::TriggerObjectRef trigRef( matchHelper.triggerMatchObject( candBaseRef, triggerMatch, iEvent, *fTriggerEvent ) );
  //  if ( trigRef.isAvailable() ) { // check references (necessary!)
  //     printf("pt  %10.2f %10.2f\n", candBaseRef->pt(), trigRef->pt() );
  //     printf("eta %10.2f %10.2f\n", candBaseRef->eta(), trigRef->eta() );
  //     printf("phi %10.2f %10.2f\n", candBaseRef->phi(), trigRef->phi() );
  //  }
  //}


  for (size_t i = 0; i != fElectrons->size(); ++i) {

    pat::Electron electron = fElectrons->at(i);
    //if (electron.electronID("eidRobustHighEnergy") != 1) {
    //  continue;
    //}
    //std::cout << "Electron Pt: " << electron.et() << std::endl;

    ThisLep.Pt = electron.et();
    ThisLep.Flavor = Dtuple::kLeptonFlavor_Electron;
    ThisLep.Id = i;
    Leptons.push_back(ThisLep);


  }



  for (size_t i = 0; i != fMuons->size(); ++i) {
    pat::Muon muon = fMuons->at(i);

    ThisLep.Pt = muon.pt();
    ThisLep.Flavor = Dtuple::kLeptonFlavor_Muon;
    ThisLep.Id = i;
    Leptons.push_back(ThisLep);

  }

  std::sort(Leptons.begin(), Leptons.end(), FillDtuple::CompareReallyBasicLeptonPt);
  //for (size_t i = 0; i != Leptons.size(); ++i) {
  //  printf("Leptons: Pt=%7.3f Flavor=%3i Id=%3i\n", Leptons[i].Pt, Leptons[i].Flavor, Leptons[i].Id);
  //}

  Ev.NLeptons = Leptons.size();
  for (size_t i = 0; (i < (int) Dtuple::NMaxLeptons) && (i != Leptons.size()); ++i) {
    switch (Leptons[i].Flavor) {
      case Dtuple::kLeptonFlavor_Electron: {
        pat::Electron electron = fElectrons->at( Leptons[i].Id );
        electron.caloPosition().eta();
        if (electron.isConvertedPhoton() != electron.isPhoton()) {
          std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
          exit(0);
        }
        // Conversion
        // IsEE IsEB ...
        Ev.Lepton_Px[i] = electron.px();
        Ev.Lepton_Py[i] = electron.py();
        Ev.Lepton_Pz[i] = electron.pz();
        Ev.Lepton_Pt[i] = electron.et();
        Ev.Lepton_TrkPt[i] = sqrt(electron.trackMomentumAtVtx().perp2());
        Ev.Lepton_Eta[i] = electron.eta();
        Ev.Lepton_Phi[i] = electron.phi();
        Ev.Lepton_dxy[i] = electron.gsfTrack()->dxy(fBeamSpot->position());
        Ev.Lepton_dz[i] = electron.gsfTrack()->dz(fBeamSpot->position());
        Ev.Lepton_Z0[i] = electron.vz();
        Ev.Lepton_Charge[i] = electron.charge();
        Ev.Lepton_Flavor[i] = Dtuple::kLeptonFlavor_Electron;
        //Ev.Lepton_TrkIso[i] = electron.trackIso();
        Ev.Lepton_TrkIso[i] = electron.dr03TkSumPt();
        Ev.Lepton_CalIso[i] = electron.caloIso();
        //Ev.Lepton_ECalIso[i] = electron.ecalIso();
        Ev.Lepton_ECalIso[i] = electron.dr03EcalRecHitSumEt();
        //Ev.Lepton_HCalIso[i] = electron.hcalIso();
        Ev.Lepton_HCalIso[i] = electron.dr03HcalTowerSumEt();
        Ev.Lepton_CalE[i] = electron.caloEnergy();
        Ev.Lepton_HCalOverECal[i] = electron.hadronicOverEm();
        Ev.Lepton_EoverPin[i] = electron.eSuperClusterOverP();
        Ev.Lepton_fBrem[i] = electron.fbrem();
        Ev.Lepton_IsConvertedPhoton[i] = electron.isConvertedPhoton();
        if (electron.electronID("eidRobustHighEnergy") == 1) {
          Ev.Lepton_PassSelection[i] |= (0x1 << Dtuple::kElectronSel_RobustHighEnergy);
        }
        if (electron.electronID("eidRobustLoose") == 1) {
          Ev.Lepton_PassSelection[i] |= (0x1 << Dtuple::kElectronSel_RobustLoose);
        }
        if (electron.electronID("eidRobustTight") == 1) {
          Ev.Lepton_PassSelection[i] |= (0x1 << Dtuple::kElectronSel_RobustTight);
        }
        if (electron.electronID("eidLoose") == 1) {
          Ev.Lepton_PassSelection[i] |= (0x1 << Dtuple::kElectronSel_Loose);
        }
        if (electron.electronID("eidTight") == 1) {
          Ev.Lepton_PassSelection[i] |= (0x1 << Dtuple::kElectronSel_Tight);
        }

        // Electron regions
        if (electron.isEE())        { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EE); }
        if (electron.isEB())        { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EB); }
        if (electron.isEBEEGap())   { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EBEEGap); }
        if (electron.isEBEtaGap())  { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EBEtaGap); }
        if (electron.isEBGap())     { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EBGap); }
        if (electron.isEBPhiGap())  { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EBPhiGap); }
        if (electron.isEEDeeGap())  { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EEDeeGap); }
        if (electron.isEEGap())     { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EEGap); }
        if (electron.isEERingGap()) { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kElectronDet_EERingGap); }

        // Electron Classification
        switch (electron.classification()) {
          case reco::GsfElectron::UNKNOWN:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_Unknown;
            break;
          case reco::GsfElectron::GOLDEN:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_Golden;
            break;
          case reco::GsfElectron::BIGBREM:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_BigBrem;
            break;
          case reco::GsfElectron::NARROW:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_Narrow;
            break;
          case reco::GsfElectron::SHOWERING:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_Showering;
            break;
          case reco::GsfElectron::GAP:
            Ev.Lepton_Classification[i] = Dtuple::kElectronClass_Gap;
            break;
          default:
            std::cerr << "ERROR in classifying lepton.  classification(): " << electron.classification() << std::endl;
            break;
        }

        Ev.Lepton_SigmaIEtaIEta[i] = electron.scSigmaIEtaIEta();
        Ev.Lepton_DeltaEtaIn[i] = electron.deltaEtaSuperClusterTrackAtVtx();
        Ev.Lepton_DeltaPhiIn[i] = electron.deltaPhiSuperClusterTrackAtVtx();
        Ev.Lepton_E2x5overE5x5[i] = electron.scE2x5Max() / electron.scE5x5();

        break;
      } case Dtuple::kLeptonFlavor_Muon: {
        pat::Muon muon = fMuons->at( Leptons[i].Id );
        Ev.Lepton_Px[i] = muon.px();
        Ev.Lepton_Py[i] = muon.py();
        Ev.Lepton_Pz[i] = muon.pz();
        Ev.Lepton_Pt[i] = muon.pt();
        Ev.Lepton_TrkPt[i] = muon.pt();
        Ev.Lepton_Eta[i] = muon.eta();
        Ev.Lepton_Phi[i] = muon.phi();
        Ev.Lepton_dxy[i] = muon.innerTrack()->dxy();
        Ev.Lepton_dz[i] = muon.innerTrack()->dz();
        Ev.Lepton_Z0[i] = muon.vz();
        Ev.Lepton_Charge[i] = muon.charge();
        Ev.Lepton_Flavor[i] = Dtuple::kLeptonFlavor_Muon;
        Ev.Lepton_TrkIso[i] = muon.trackIso();
        Ev.Lepton_CalIso[i] = muon.caloIso();
        Ev.Lepton_ECalIso[i] = muon.ecalIso();
        Ev.Lepton_HCalIso[i] = muon.hcalIso();
        Ev.Lepton_CalE[i] = muon.calEnergy().em + muon.calEnergy().had;
        Ev.Lepton_HCalOverECal[i] = muon.calEnergy().had / muon.calEnergy().em;

        // Muon types
        if (muon.isGlobalMuon())     { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kMuonDet_Global); }
        if (muon.isTrackerMuon())    { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kMuonDet_Tracker); }
        if (muon.isStandAloneMuon()) { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kMuonDet_StandAlone); }
        if (muon.isCaloMuon())       { Ev.Lepton_Detector[i] |= (0x1 << Dtuple::kMuonDet_Calo); }
        break;
      } default: {
        std::cerr << "ERROR: lepton type not known" << std::endl;
        break;
      }
    }
  }


  return;
}


// ------------ Fill the jets  ------------
void 
FillDtuple::FillJets(const edm::Event& iEvent, DtupleWriter::Event_Struct& Ev)
{

  Ev.NJets = fJets->size();
  for (size_t i = 0; i != fJets->size() && i < (size_t) Dtuple::NMaxJets; ++i) {

    pat::Jet jet = fJets->at(i);

    Ev.Jet_Px[i] = jet.px();
    Ev.Jet_Py[i] = jet.py();
    Ev.Jet_Pz[i] = jet.pz();
    Ev.Jet_Pt[i] = jet.pt();
    Ev.Jet_Eta[i] = jet.eta();
    Ev.Jet_Phi[i] = jet.phi();
    Ev.Jet_EmF[i] = jet.emEnergyFraction();
    Ev.Jet_HadF[i] = jet.energyFractionHadronic();


  }

  return;
}



// ------------ Fill the electrons  ------------
void 
FillDtuple::FillPhotons (const edm::Event& iEvent, DtupleWriter::Event_Struct& Ev)
{
  Ev.NPhotons = fPhotons->size();
  for (size_t i = 0; i != fPhotons->size() && i < (size_t) Dtuple::NMaxPhotons; ++i) {

    pat::Photon photon = fPhotons->at(i);

    Ev.Photon_Px[i] = photon.px();
    Ev.Photon_Py[i] = photon.py();
    Ev.Photon_Pz[i] = photon.pz();
    Ev.Photon_Pt[i] = photon.pt();
    Ev.Photon_Eta[i] = photon.eta();
    Ev.Photon_Phi[i] = photon.phi();
    Ev.Photon_TrkIso[i] = photon.trackIso();
    Ev.Photon_CalIso[i] = photon.caloIso();
    Ev.Photon_HCalOverECal[i] = photon.hadronicOverEm();
  }

  return;
}



// ------------ print brief event summary ------------
void 
FillDtuple::EventSummary (DtupleWriter::Event_Struct& Ev)
{
  // Print a very basic event summary
  printf("Leptons: %5i Photons: %5i Jets: %5i Met: %8.1f\n", Ev.NLeptons, Ev.NPhotons, Ev.NJets, Ev.MetMag);
  return;
}




// ------------ Keep this event or not? ------------
bool
FillDtuple::KeepEvent (DtupleWriter::Event_Struct& Ev)
{
  bool keep = true;
  if (Ev.NLeptons < 1) {
    keep = false;
  }

  return keep;
}




// ------------ method called once each job just before starting event loop  ------------
void 
FillDtuple::beginJob()
{
  edm::Service<TFileService> fs;
  TTree* DtupleTree = fs->make<TTree>("Dtuple", "Dtuple");
  DtupleTree->SetDirectory( &(fs->file()) );

  fDtupleWriter = new DtupleWriter(DtupleTree);
  //fDtupleWriter = new DtupleWriter(fOutFileName);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FillDtuple::endJob() {
  //delete fDtupleWriter;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FillDtuple);
