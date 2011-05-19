#ifndef GUARD_DHidasPatAna_h
#define GUARD_DHidasPatAna_h

#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "TString.h" 
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH3D.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TBranch.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"


#include "DHidasLJAna/LeptonPlusJets/interface/Dtuple.h"
#include "DHidasLJAna/LeptonPlusJets/interface/DHidasJSON.h"


class DHidasPatAna : public edm::EDAnalyzer, public Dtuple
{
  public:
    explicit DHidasPatAna (const edm::ParameterSet&);
    ~DHidasPatAna() {};

    bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup);

  private:
    //virtual void beginJob(const edm::EventSetup&) ;
    float correct_met (float et, float oc, float phi, float eta, float metin, float metphiin, float met, float &cormetphi);
    virtual void beginJob ();
    virtual void analyze (const edm::Event&, const edm::EventSetup&);
    virtual void endJob ();

    void GetObjects (const edm::Event& iEvent);
    void PlotObjects ();
    void PlotDileptonEvents ();
    void PlotMultiJetLeptonEvents ();

    void FillTree ();


    // Declaration of leaf types
    Int_t           fRun;
    Int_t           fEvent;
    Int_t           fLumiSection;
    bool            fIsData;
    bool            fMakeDtuple;

    // For output
    TFile* fOutFile;
    TTree* fTree;
    TString fOutFileName;

    // For JSON
    std::string fJSONFilename;
    DHidasJSON  fJSON;

    // Object handles, object, and clean objects
    edm::Handle< std::vector<pat::Electron> > PatElectrons; 
    //edm::Handle< reco::TrackCollection >      recoTracks;
    edm::Handle< std::vector<pat::Photon> >   PatPhotons; 
    edm::Handle< std::vector<pat::Muon> >     PatMuons; 
    edm::Handle< std::vector<pat::Jet> >      PatJets;
    edm::Handle< std::vector<pat::MET> >      MetColl;  
    //edm::Handle<reco::VertexCollection>       recVtxs;
    edm::Handle< trigger::TriggerEvent >      TriggerEvent;

    std::vector<const pat::Electron*> fGoodElectrons;
    std::vector<const pat::Photon*>   fGoodPhotons;
    std::vector<const pat::Muon*>     fGoodMuons;
    std::vector<const pat::Jet*>      fGoodJets;

    std::vector<const pat::Muon*>     fCleanMuons;
    std::vector<const pat::Electron*> fCleanElectrons;
    std::vector<const pat::Photon*>   fCleanPhotons;
    std::vector<const pat::Jet*>      fCleanJets;

    std::vector<std::string>          fTriggerNames;
    std::map<std::string, bool>       fTriggerMap;

    void getTriggerDecision(const edm::Event&, std::map<std::string, bool>&);
    int
      matchElectrons(const trigger::TriggerEvent&    triggerEvent,
          const std::vector<std::string>& module,
          std::vector<const pat::Electron*>& Electrons)
      {
        // what DR cut to use?
        int const DeltaRMatch = 0.4;

        // Value to return.  bits which are 1 means fired trigger
        int Ret = 0x0;

        // which menu. This shoudl really go elsewhere..
        std::string menu = "HLT";

        // Loop over these leptons
        for (size_t i = 0; i < Electrons.size(); i++) {
        std::cout << "HERE" << std::endl;

          // Loop over the possible paths
          for (size_t j = 0; j < module.size(); j++) {
            edm::InputTag filterTag(module.at(j), "", menu);

            std::vector<trigger::TriggerObject> triggerObjects;

            const trigger::TriggerObjectCollection triggerObjectCollection = triggerEvent.getObjects();

            const trigger::Keys& keys = triggerEvent.collectionKeys();
            for (size_t k = 0; k < keys.size(); k++) {
              //std::cout << triggerObjectCollection[keys[k]].energy() << std::endl;
            }
            if (triggerEvent.filterIndex(filterTag) < triggerEvent.sizeFilters()) {
              const trigger::Keys& keys = triggerEvent.filterKeys(triggerEvent.filterIndex(filterTag));

              for (size_t k = 0; k < keys.size(); k++) {
                triggerObjects.push_back(triggerObjectCollection[keys[k]]);
              }
            }

            for (size_t k = 0; k < triggerObjects.size(); k++) {
              double deltaR = reco::deltaR(Electrons[i]->eta(), Electrons[i]->phi(), triggerObjects.at(k).eta(), triggerObjects.at(k).phi());

              if (deltaR < DeltaRMatch) {
                Ret |= (0x1 << i);
                printf("i Ret  %i  %i\n", (int) i, Ret);
                break;
              }
            }
          }
        }
        return Ret;
      }

};



#endif
