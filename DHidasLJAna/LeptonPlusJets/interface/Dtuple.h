#ifndef GUARD_Dtuple_h
#define GUARD_Dtuple_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Apr 26 03:47:29 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"


class Dtuple
{
  public:
    Dtuple ();
    Dtuple (std::vector<TString> const&);
    ~Dtuple ();

    enum LepType {
      kElectron = 0,
      kMuon,
      kNLeptonTypes
    };

    void SetBranchAddresses (TTree*);
    void SetBranches (TTree*);
    void ClearDtuple ();
    int  GetEntry (long long&);


  private:
    TTree* fTree;

  public:
    struct SimpleEvent {
      int Run;
      int LumiSection;
      int Event;


      std::vector<float>* LeptonPx;
      std::vector<float>* LeptonPy;
      std::vector<float>* LeptonPz;
      std::vector<float>* LeptonPt;
      std::vector<int>*   LeptonType;
      int   NLeptons;

      std::vector<float>* JetPx;
      std::vector<float>* JetPy;
      std::vector<float>* JetPz;
      std::vector<float>* JetPt;
      int   NJets;

      float SumPtJets;
      std::vector<float>* TriJetSumPt;
      std::vector<float>* TriJetMasses;

      float METMag;
      float METPhi;


      // Not stored in Dtuple...but useful
      std::vector<TLorentzVector> Lep;
      std::vector<TLorentzVector> Jet;
      TVector2 MET;
    };
    SimpleEvent& GetEvt () {
      return fEvt;
    }


  protected:
    SimpleEvent fEvt;
};







#endif
