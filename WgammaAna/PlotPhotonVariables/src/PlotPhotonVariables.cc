// -*- C++ -*-
//
// Package:    PlotPhotonVariables
// Class:      PlotPhotonVariables
// 
/**\class PlotPhotonVariables PlotPhotonVariables.cc WgammaAna/PlotPhotonVariables/src/PlotPhotonVariables.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dean Andrew Hidas,8 R-020,+41227673494,
//         Created:  Tue Jun 29 17:07:12 CEST 2010
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

class PlotPhotonVariables : public edm::EDAnalyzer {
  public:
    explicit PlotPhotonVariables(const edm::ParameterSet&);
    ~PlotPhotonVariables();


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    // Input Tags
    edm::InputTag photonSrc_;

    // Histograms
    std::map<std::string, TH1D*> Hist1D;
};

#include "DataFormats/PatCandidates/interface/Photon.h"
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
PlotPhotonVariables::PlotPhotonVariables(const edm::ParameterSet& iConfig) :
photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))

{
   //now do what ever initialization is needed

}


PlotPhotonVariables::~PlotPhotonVariables()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PlotPhotonVariables::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByLabel(photonSrc_, photons);


  for (edm::View<pat::Photon>::const_iterator photon = photons->begin(); photon != photons->end(); ++photon) {
    Hist1D["Photon_eta"]->Fill(photon->eta());
  }


  return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
PlotPhotonVariables::beginJob()
{
  std::cout << "PlotPhotonVariables::beginJob" << std::endl;
  edm::Service<TFileService> fs;

  Hist1D["Photon_eta"] = fs->make<TH1D>("Photon_eta", "eta", 100, -3, 3);

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PlotPhotonVariables::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PlotPhotonVariables);
