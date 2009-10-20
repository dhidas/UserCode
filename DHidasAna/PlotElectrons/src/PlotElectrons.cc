// -*- C++ -*-
//
// Package:    PlotElectrons
// Class:      PlotElectrons
// 
/**\class PlotElectrons PlotElectrons.cc DHidasAna/PlotElectrons/src/PlotElectrons.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dean Andrew Hidas
//         Created:  Wed Sep 23 10:33:20 CEST 2009
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
#include "TH2D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


//
// class decleration
//

class PlotElectrons : public edm::EDAnalyzer {
   public:
      explicit PlotElectrons(const edm::ParameterSet&);
      ~PlotElectrons();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

    // Input Tags
    edm::InputTag electronSrc_;
    edm::InputTag metSrc_;

    // Histograms
    std::map<std::string, TH1D*> Hist1D;
    std::map<std::string, TH2D*> Hist2D;

};



#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TMath.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PlotElectrons::PlotElectrons(const edm::ParameterSet& iConfig) :
electronSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))
{
   //now do what ever initialization is needed

}


PlotElectrons::~PlotElectrons()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PlotElectrons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Get the damn electrons
  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);

  // MET object that corrects the basic calorimeter MET for muons
  edm::Handle< edm::View<pat::MET> > MET;
  iEvent.getByLabel(metSrc_, MET);


  Hist1D["NElectrons"]->Fill(electrons->size());

  if (electrons->size() > 0) {
    Hist1D["LeadingElectron_et"]->Fill( electrons->begin()->et() );
    Hist1D["Mt"]->Fill( TMath::Sqrt( 2.0 * electrons->begin()->et() * (MET->front()).et() * (1.0 - TMath::Cos(electrons->begin()->phi() - (MET->front()).phi() ))) );
  }

  // Loop over electrons
  int iElectron = 0;
  for (edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {
    Hist1D["et"]->Fill(electron->et());
    Hist1D["eta"]->Fill(electron->eta());
    Hist1D["phi"]->Fill(electron->phi());
    Hist1D["numberOfTracks"]->Fill(electron->numberOfTracks());
    Hist1D["d0"]->Fill(electron->gsfTrack()->d0());
    Hist1D["met"]->Fill((MET->front()).et());
    Hist1D["trackIso"]->Fill(electron->trackIso());
    Hist1D["caloIso"]->Fill(electron->caloIso());

    Hist2D["Electron_EtaPhi"]->Fill( electron->eta(), electron->phi() );
    Hist2D["Electron_etVsN"]->Fill( electron->et(), ++iElectron );

    
  }
  

  return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
PlotElectrons::beginJob()
{
  std::cout << "PlotElectrons::beginJob" << std::endl;
  edm::Service<TFileService> fs;

  // 1D plots
  Hist1D["NElectrons"]          = fs->make<TH1D>("NElectrons", "NElectrons", 10, 0, 10);
  Hist1D["et"]                  = fs->make<TH1D>("et", "et", 100, 0, 150);
  Hist1D["LeadingElectron_et"]  = fs->make<TH1D>("LeadingElectron_et", "LeadingElectron_et", 100, 0, 150);
  Hist1D["eta"]                 = fs->make<TH1D>("eta", "eta", 100, -3.0, 3.0);
  Hist1D["phi"]                 = fs->make<TH1D>("phi", "phi", 100, -3.2, 3.2);
  Hist1D["numberOfTracks"]      = fs->make<TH1D>("numberOfTracks", "numberOfTracks", 20, 0, 20);
  Hist1D["d0"]                  = fs->make<TH1D>("d0", "d0", 100, -0.3, 0.3);
  Hist1D["trackIso"]            = fs->make<TH1D>("trackIso", "trackIso", 100, 0, 20);
  Hist1D["caloIso"]             = fs->make<TH1D>("caloIso", "caloIso", 100, 0, 20);
  Hist1D["met"]                 = fs->make<TH1D>("met", "met", 100, 0, 100);
  Hist1D["Mt"]                  = fs->make<TH1D>("Mt", "Mt", 100, 0, 100);

  // 2D plots
  Hist2D["Electron_EtaPhi"]     = fs->make<TH2D>("Electron_EtaPhi", "Electron_EtaPhi", 200, -3.0, 3.0, 200, -3.2, 3.2);
  Hist2D["Electron_etVsN"]      = fs->make<TH2D>("Electron_etVsN", "Electron_etVsN", 100, 0, 150, 10, 0, 10);

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PlotElectrons::endJob() {
  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PlotElectrons);
