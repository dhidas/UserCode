// -*- C++ -*-
//
// Package:    Filter
// Class:      Filter
// 
/**\class Filter Filter.cc DHidasLJAna/Filter/src/Filter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dean Hidas
//         Created:  Tue Apr 19 12:07:58 EDT 2011
// $Id: Filter.cc,v 1.2 2011/05/03 09:18:15 dhidas Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DHidasLJAna/LeptonPlusJets/plugins/DHidasPatAna.h"

//
// class declaration
//

class Filter : public edm::EDFilter {
   public:
      explicit Filter(const edm::ParameterSet&);
      ~Filter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      DHidasPatAna* Ana;

      unsigned long fTotalEvents;
      unsigned long fPassedTrigger;
      unsigned long fFailedTrigger;
      
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
Filter::Filter(const edm::ParameterSet& iConfig)
{
  Ana = new DHidasPatAna(iConfig);
   //now do what ever initialization is needed

}


Filter::~Filter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete Ana;
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool const Pass = Ana->filter(iEvent, iSetup);
  ++fTotalEvents;
  if (Pass) {
    ++fPassedTrigger;
  } else {
    ++fFailedTrigger;
  }

  return Pass;
  //return Ana->filter(iEvent, iSetup);
}

// ------------ method called once each job just before starting event loop  ------------
void 
Filter::beginJob()
{
  // Test for triggers...
  fTotalEvents = 0;
  fPassedTrigger = 0;
  fFailedTrigger = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Filter::endJob() {
  printf("Trigger: Total Pass Fail  %12lu %12lu %12lu\n", fTotalEvents, fPassedTrigger, fFailedTrigger);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Filter);
