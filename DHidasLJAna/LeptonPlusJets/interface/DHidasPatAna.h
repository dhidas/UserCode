#include <memory>
#include <iostream>

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

#include "DHidasLJAna/LeptonPlusJets/interface/Dtuple.h"


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
    std::string     fJSONFilename;
    Int_t fGoodRuns[4000];
    Int_t fGoodLumiStart[4000];
    Int_t fGoodLumiEnd[4000];
    Int_t fNGoodRuns;

    // Object handles, object, and clean objects
    edm::Handle< std::vector<pat::Electron> > PatElectrons; 
    edm::Handle<reco::TrackCollection> recoTracks;
    edm::Handle< std::vector<pat::Photon> > PatPhotons; 
    edm::Handle< std::vector<pat::Muon> > PatMuons; 
    edm::Handle< std::vector<pat::Jet> > PatJets;
    edm::Handle< std::vector<pat::MET> > MetColl;  
    edm::Handle<reco::VertexCollection> recVtxs;

    std::vector<const pat::Electron*> fGoodElectrons;
    std::vector<const pat::Photon*>   fGoodPhotons;
    std::vector<const pat::Muon*>     fGoodMuons;
    std::vector<const pat::Jet*>      fGoodJets;

    std::vector<const pat::Muon*>     fCleanMuons;
    std::vector<const pat::Electron*> fCleanElectrons;
    std::vector<const pat::Photon*>   fCleanPhotons;
    std::vector<const pat::Jet*>      fCleanJets;
};

