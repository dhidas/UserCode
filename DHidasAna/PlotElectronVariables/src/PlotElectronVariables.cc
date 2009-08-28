// -*- C++ -*-
//
// Package:    PlotElectronVariables
// Class:      PlotElectronVariables
// 
/**\class PlotElectronVariables PlotElectronVariables.cc DHidasAna/PlotElectronVariables/src/PlotElectronVariables.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Dean Andrew Hidas
//         Created:  Thu Aug 20 15:45:17 CEST 2009
// $Id$
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

// dhidas includes
#include <string>

#include "TH1D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


//
// class decleration
//

class PlotElectronVariables : public edm::EDAnalyzer {
  public:
    explicit PlotElectronVariables(const edm::ParameterSet&);
    ~PlotElectronVariables();


  private:
    virtual void beginJob(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // Input Tags
    edm::InputTag electronSrc_;

    // Histograms
    std::map<std::string, TH1D*> Hist1D;

    // ----------member data ---------------------------
};


#include "DataFormats/PatCandidates/interface/Electron.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PlotElectronVariables::PlotElectronVariables(const edm::ParameterSet& iConfig) :
electronSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))
{
  //now do what ever initialization is needed

}


PlotElectronVariables::~PlotElectronVariables()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
  void
PlotElectronVariables::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


  // Get the damn electrons
  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);

  for (edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {
    Hist1D["et"]->Fill(electron->et());
    Hist1D["eta"]->Fill(electron->eta());
    Hist1D["eSuperClusterOverP"]->Fill(electron->eSuperClusterOverP());
    //Hist1D["eSeed"]->Fill(electron->superCluster()->seed()->energy());
    Hist1D["pin"]->Fill(electron->trackMomentumAtVtx().R());
    //Hist1D["eSeedOverPin"]->Fill(electron->superCluster()->seed()->energy()/electron->trackMomentumAtVtx().R());
    Hist1D["pout"]->Fill(electron->trackMomentumOut().R());
    Hist1D["fBrem"]->Fill( (electron->trackMomentumAtVtx().R() - electron->trackMomentumOut().R())/electron->trackMomentumAtVtx().R());
    Hist1D["hOverE"]->Fill(electron->hadronicOverEm());
    Hist1D["deltaPhiIn"]->Fill(electron->deltaPhiSuperClusterTrackAtVtx());
    Hist1D["deltaEtaIn"]->Fill(electron->deltaEtaSuperClusterTrackAtVtx());
    Hist1D["trackIso"]->Fill(electron->trackIso());
    Hist1D["caloIso"]->Fill(electron->caloIso());
    Hist1D["ecalIso"]->Fill(electron->ecalIso());
    Hist1D["hcalIso"]->Fill(electron->hcalIso());
    Hist1D["userIso"]->Fill(electron->userIso());
  }


}


// ------------ method called once each job just before starting event loop  ------------
  void 
PlotElectronVariables::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  Hist1D["et"] = fs->make<TH1D>("et", "et", 100, 0, 200);
  Hist1D["eta"] = fs->make<TH1D>("eta", "eta", 100, -3, 3);
  Hist1D["eSuperClusterOverP"] = fs->make<TH1D>("eSuperClusterOverP", "eSuperClusterOverP", 100, 0, 3);
  Hist1D["eSeed"] = fs->make<TH1D>("eSeed", "eSeed", 100, 0, 200);
  Hist1D["pin"] = fs->make<TH1D>("pin", "pin", 100, 0, 200);
  Hist1D["peSeedOverPin"] = fs->make<TH1D>("peSeedOverPin", "peSeedOverPin", 100, 0, 3);
  Hist1D["pout"] = fs->make<TH1D>("pout", "pout", 100, 0, 200);
  Hist1D["fBrem"] = fs->make<TH1D>("fBrem", "fBrem", 100, 0, 3);
  Hist1D["hOverE"] = fs->make<TH1D>("hOverE", "hOverE", 100, 0, 2);
  Hist1D["deltaPhiIn"] = fs->make<TH1D>("deltaPhiIn", "deltaPhiIn", 100, -1, 1);
  Hist1D["deltaEtaIn"] = fs->make<TH1D>("deltaEtaIn", "deltaEtaIn", 100, -1, 1);
  Hist1D["trackIso"] = fs->make<TH1D>("trackIso", "trackIso", 100, 0, 100);
  Hist1D["caloIso"] = fs->make<TH1D>("caloIso", "caloIso", 100, 0, 100);
  Hist1D["ecalIso"] = fs->make<TH1D>("ecalIso", "ecalIso", 100, 0, 100);
  Hist1D["hcalIso"] = fs->make<TH1D>("hcalIso", "hcalIso", 100, 0, 100);
  Hist1D["userIso"] = fs->make<TH1D>("userIso", "hcalIso", 100, 0, 100);

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PlotElectronVariables::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PlotElectronVariables);
