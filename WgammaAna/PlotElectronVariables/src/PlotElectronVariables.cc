// -*- C++ -*-
//
// Package:    PlotElectronVariables
// Class:      PlotElectronVariables
// 
/**\class PlotElectronVariables PlotElectronVariables.cc DHidasAna/PlotElectronVariables/src/PlotElectronVariables.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dean Andrew Hidas,8 R-020,+41227673494,
//         Created:  Mon Jun 28 17:00:35 CEST 2010
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TH1D.h"
#include <string>


//
// class declaration
//

class PlotElectronVariables : public edm::EDAnalyzer {
  public:
    explicit PlotElectronVariables(const edm::ParameterSet&);
    ~PlotElectronVariables();


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    // Input Tags
    edm::InputTag electronSrc_;

    // Histograms
    std::map<std::string, TH1D*> Hist1D;
};


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/Ref.h"


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
  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);


  for (edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {

    if (electron->electronID("eidRobustHighEnergy") != 1) {
      continue;
    }


    if (electron->isEB()) {
      Hist1D["Barrel_hOverE"]->Fill(electron->hadronicOverEm());
      Hist1D["Barrel_deltaPhiIn"]->Fill( electron->deltaPhiSuperClusterTrackAtVtx() );
      Hist1D["Barrel_deltaEtaIn"]->Fill( electron->deltaEtaSuperClusterTrackAtVtx() );

      Hist1D["Barrel_et"]->Fill(electron->et());
      Hist1D["Barrel_eta"]->Fill(electron->eta());
      Hist1D["Barrel_eSuperClusterOverP"]->Fill(electron->eSuperClusterOverP());

      Hist1D["Barrel_pin"]->Fill(electron->trackMomentumAtVtx().R());
      Hist1D["Barrel_pout"]->Fill(electron->trackMomentumOut().R());
      Hist1D["Barrel_fBrem"]->Fill( (electron->trackMomentumAtVtx().R() - electron->trackMomentumOut().R())/electron->trackMomentumAtVtx().R());
      Hist1D["Barrel_trackIso"]->Fill(electron->trackIso());
      Hist1D["Barrel_caloIso"]->Fill(electron->caloIso());
      Hist1D["Barrel_ecalIso"]->Fill(electron->ecalIso());
      Hist1D["Barrel_hcalIso"]->Fill(electron->hcalIso());
      Hist1D["Barrel_userIso"]->Fill(electron->userIso());
    } else if (false) {
      Hist1D["Endcap_hOverE"]->Fill(electron->hadronicOverEm());
      Hist1D["Endcap_deltaPhiIn"]->Fill( electron->deltaPhiSuperClusterTrackAtVtx() );
      Hist1D["Endcap_deltaEtaIn"]->Fill( electron->deltaEtaSuperClusterTrackAtVtx() );

      Hist1D["Endcap_et"]->Fill(electron->et());
      Hist1D["Endcap_eta"]->Fill(electron->eta());
      Hist1D["Endcap_eSuperClusterOverP"]->Fill(electron->eSuperClusterOverP());
      Hist1D["Endcap_eSeed"]->Fill(electron->superCluster()->seed()->energy());
      Hist1D["Endcap_pin"]->Fill(electron->trackMomentumAtVtx().R());
      Hist1D["Endcap_eSeedOverPin"]->Fill(electron->superCluster()->seed()->energy()/electron->trackMomentumAtVtx().R());
      Hist1D["Endcap_pout"]->Fill(electron->trackMomentumOut().R());
      Hist1D["Endcap_fBrem"]->Fill( (electron->trackMomentumAtVtx().R() - electron->trackMomentumOut().R())/electron->trackMomentumAtVtx().R());
      Hist1D["Endcap_trackIso"]->Fill(electron->trackIso());
      Hist1D["Endcap_caloIso"]->Fill(electron->caloIso());
      Hist1D["Endcap_ecalIso"]->Fill(electron->ecalIso());
      Hist1D["Endcap_hcalIso"]->Fill(electron->hcalIso());
      Hist1D["Endcap_userIso"]->Fill(electron->userIso());
    }
  }



}


// ------------ method called once each job just before starting event loop  ------------
void 
PlotElectronVariables::beginJob()
{
  std::cout << "PlotElectronVariables::beginJob" << std::endl;
  edm::Service<TFileService> fs;

  Hist1D["Barrel_deltaPhiIn"] = fs->make<TH1D>("Barrel_deltaPhiIn", "deltaPhiIn", 100, -0.2, 0.2);
  Hist1D["Barrel_deltaEtaIn"] = fs->make<TH1D>("Barrel_deltaEtaIn", "deltaEtaIn", 100, -0.05, 0.05);

  Hist1D["Barrel_et"] = fs->make<TH1D>("Barrel_et", "et", 100, 0, 150);
  Hist1D["Barrel_eta"] = fs->make<TH1D>("Barrel_eta", "eta", 100, -3, 3);
  Hist1D["Barrel_eSuperClusterOverP"] = fs->make<TH1D>("Barrel_eSuperClusterOverP", "eSuperClusterOverP", 100, 0, 3);
  Hist1D["Barrel_eSeed"] = fs->make<TH1D>("Barrel_eSeed", "eSeed", 100, 0, 200);
  Hist1D["Barrel_pin"] = fs->make<TH1D>("Barrel_pin", "pin", 100, 0, 200);
  Hist1D["Barrel_eSeedOverPin"] = fs->make<TH1D>("Barrel_eSeedOverPin", "eSeedOverPin", 100, 0, 3);
  Hist1D["Barrel_pout"] = fs->make<TH1D>("Barrel_pout", "pout", 100, 0, 200);
  Hist1D["Barrel_fBrem"] = fs->make<TH1D>("Barrel_fBrem", "fBrem", 100, 0, 1.1);
  Hist1D["Barrel_hOverE"] = fs->make<TH1D>("Barrel_hOverE", "hOverE", 100, 0, 0.2);
  Hist1D["Barrel_trackIso"] = fs->make<TH1D>("Barrel_trackIso", "trackIso", 100, 0, 20);
  Hist1D["Barrel_caloIso"] = fs->make<TH1D>("Barrel_caloIso", "caloIso", 100, 0, 20);
  Hist1D["Barrel_ecalIso"] = fs->make<TH1D>("Barrel_ecalIso", "ecalIso", 100, 0, 20);
  Hist1D["Barrel_hcalIso"] = fs->make<TH1D>("Barrel_hcalIso", "hcalIso", 100, 0, 20);
  Hist1D["Barrel_userIso"] = fs->make<TH1D>("Barrel_userIso", "userIso", 100, 0, 20);

  Hist1D["Endcap_deltaPhiIn"] = fs->make<TH1D>("Endcap_deltaPhiIn", "deltaPhiIn", 100, -0.2, 0.2);
  Hist1D["Endcap_deltaEtaIn"] = fs->make<TH1D>("Endcap_deltaEtaIn", "deltaEtaIn", 100, -0.05, 0.05);

  Hist1D["Endcap_et"] = fs->make<TH1D>("Endcap_et", "et", 100, 0, 150);
  Hist1D["Endcap_eta"] = fs->make<TH1D>("Endcap_eta", "eta", 100, -3, 3);
  Hist1D["Endcap_eSuperClusterOverP"] = fs->make<TH1D>("Endcap_eSuperClusterOverP", "eSuperClusterOverP", 100, 0, 3);
  Hist1D["Endcap_eSeed"] = fs->make<TH1D>("Endcap_eSeed", "eSeed", 100, 0, 200);
  Hist1D["Endcap_pin"] = fs->make<TH1D>("Endcap_pin", "pin", 100, 0, 200);
  Hist1D["Endcap_eSeedOverPin"] = fs->make<TH1D>("Endcap_eSeedOverPin", "eSeedOverPin", 100, 0, 3);
  Hist1D["Endcap_pout"] = fs->make<TH1D>("Endcap_pout", "pout", 100, 0, 200);
  Hist1D["Endcap_fBrem"] = fs->make<TH1D>("Endcap_fBrem", "fBrem", 100, 0, 1.1);
  Hist1D["Endcap_hOverE"] = fs->make<TH1D>("Endcap_hOverE", "hOverE", 100, 0, 0.2);
  Hist1D["Endcap_trackIso"] = fs->make<TH1D>("Endcap_trackIso", "trackIso", 100, 0, 20);
  Hist1D["Endcap_caloIso"] = fs->make<TH1D>("Endcap_caloIso", "caloIso", 100, 0, 20);
  Hist1D["Endcap_ecalIso"] = fs->make<TH1D>("Endcap_ecalIso", "ecalIso", 100, 0, 20);
  Hist1D["Endcap_hcalIso"] = fs->make<TH1D>("Endcap_hcalIso", "hcalIso", 100, 0, 20);
  Hist1D["Endcap_userIso"] = fs->make<TH1D>("Endcap_userIso", "userIso", 100, 0, 20);

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PlotElectronVariables::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PlotElectronVariables);
