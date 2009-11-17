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
// $Id: FillDtuple.cc,v 1.6 2009/11/17 15:09:24 dhidas Exp $
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
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"





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
  iEvent.getByLabel("cleanLayer1Electrons", fElectrons);
  iEvent.getByLabel("cleanLayer1Muons", fMuons);
  iEvent.getByLabel("cleanLayer1Jets", fJets);
  iEvent.getByLabel("cleanLayer1Photons", fPhotons);
  iEvent.getByLabel("layer1METs", fMETs);
  iEvent.getByLabel("patTrigger", fTriggerEvent );
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
  edm::Handle< pat::MuonCollection > muons;
  iEvent.getByLabel( "selectedLayer1Muons", muons );
  const pat::helper::TriggerMatchHelper matchHelper;
  const pat::TriggerObjectMatch* triggerMatch( fTriggerEvent->triggerObjectMatchResult( "muonTriggerMatchHLTMuons" ) );
  for (size_t i = 0; i != muons->size(); ++i) {
    const reco::CandidateBaseRef candBaseRef( pat::MuonRef( muons, i ) );
    const pat::TriggerObjectRef trigRef( matchHelper.triggerMatchObject( candBaseRef, triggerMatch, iEvent, *fTriggerEvent ) );
    if ( trigRef.isAvailable() ) { // check references (necessary!)
       printf("pt  %10.2f %10.2f\n", candBaseRef->pt(), trigRef->pt() );
       printf("eta %10.2f %10.2f\n", candBaseRef->eta(), trigRef->eta() );
       printf("phi %10.2f %10.2f\n", candBaseRef->phi(), trigRef->phi() );
    }
  }


  for (size_t i = 0; i != fElectrons->size(); ++i) {

    pat::Electron electron = fElectrons->at(i);
    //if (electron.electronID("eidRobustHighEnergy") != 1) {
    //  continue;
    //}
    //std::cout << "Electron Pt: " << electron.et() << std::endl;

    ThisLep.Pt = electron.et();
    ThisLep.Flavor = Dtuple::kElectron;
    ThisLep.Id = i;
    Leptons.push_back(ThisLep);


  }



  for (size_t i = 0; i != fMuons->size(); ++i) {
    pat::Muon muon = fMuons->at(i);

    ThisLep.Pt = muon.pt();
    ThisLep.Flavor = Dtuple::kMuon;
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
      case Dtuple::kElectron: {
        pat::Electron electron = fElectrons->at( Leptons[i].Id );
        Ev.Lepton_Px[i] = electron.px();
        Ev.Lepton_Py[i] = electron.py();
        Ev.Lepton_Pz[i] = electron.pz();
        Ev.Lepton_Pt[i] = electron.et();
        Ev.Lepton_TrkPt[i] = sqrt(electron.trackMomentumAtVtx().perp2());
        Ev.Lepton_Eta[i] = electron.eta();
        Ev.Lepton_Phi[i] = electron.phi();
        Ev.Lepton_D0[i] = 0;
        Ev.Lepton_Z0[i] = 0;
        Ev.Lepton_Charge[i] = electron.charge();
        Ev.Lepton_Flavor[i] = Dtuple::kElectron;
        Ev.Lepton_TrkIso[i] = electron.trackIso();
        Ev.Lepton_CalIso[i] = electron.caloIso();
        Ev.Lepton_ECalIso[i] = electron.ecalIso();
        Ev.Lepton_HCalIso[i] = electron.hcalIso();
        Ev.Lepton_CalE[i] = electron.caloEnergy();
        Ev.Lepton_HCalOverECal[i] = electron.hadronicOverEm();
        Ev.Lepton_EoverPin[i] = electron.eSuperClusterOverP();
        break;
      } case Dtuple::kMuon: {
        pat::Muon muon = fMuons->at( Leptons[i].Id );
        Ev.Lepton_Px[i] = muon.px();
        Ev.Lepton_Py[i] = muon.py();
        Ev.Lepton_Pz[i] = muon.pz();
        Ev.Lepton_Pt[i] = muon.pt();
        Ev.Lepton_TrkPt[i] = muon.pt();
        Ev.Lepton_Eta[i] = muon.eta();
        Ev.Lepton_Phi[i] = muon.phi();
        Ev.Lepton_D0[i] = 0;
        Ev.Lepton_Z0[i] = 0;
        Ev.Lepton_Charge[i] = muon.charge();
        Ev.Lepton_Flavor[i] = Dtuple::kMuon;
        Ev.Lepton_TrkIso[i] = muon.trackIso();
        Ev.Lepton_CalIso[i] = muon.caloIso();
        Ev.Lepton_ECalIso[i] = muon.ecalIso();
        Ev.Lepton_HCalIso[i] = muon.hcalIso();
        Ev.Lepton_CalE[i] = muon.calEnergy().em + muon.calEnergy().had;
        Ev.Lepton_HCalOverECal[i] = muon.calEnergy().had / muon.calEnergy().em;
        Ev.Lepton_EoverPin[i] = -1;
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
  if (Ev.NLeptons != 1) {
    keep = false;
  }
  if (Ev.NPhotons != 1) {
    keep = false;
  }
  return keep;
}




// ------------ method called once each job just before starting event loop  ------------
void 
FillDtuple::beginJob()
{
  fDtupleWriter = new DtupleWriter("TestDtupleFile.root");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FillDtuple::endJob() {
  //delete fDtupleWriter;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FillDtuple);
