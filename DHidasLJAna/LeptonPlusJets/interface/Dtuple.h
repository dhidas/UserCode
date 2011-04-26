////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Apr 26 03:47:29 EDT 2011
//
////////////////////////////////////////////////////////////////////

#include <iostream>

#include "TString.h"
#include "TTree.h"


class Dtuple
{
  public:
    Dtuple ();
    ~Dtuple ();

    static int const kMaxLeptons      =  4;
    static int const kMaxJets         = 10;
    static int const kMaxCombinations = 26;

    enum LepType {
      kElectron = 0,
      kMuon,
      kNLeptonTypes
    };

    void SetBranches (TTree*);
    void ClearDtuple ();


  protected:
    struct SimpleEvent {
      int Run;
      int LumiSection;
      int Event;

      float LeptonPx[kMaxLeptons];
      float LeptonPy[kMaxLeptons];
      float LeptonPz[kMaxLeptons];
      float LeptonPt[kMaxLeptons];
      int   LeptonType[kMaxLeptons];
      int   NLeptons;

      float JetPx[kMaxJets];
      float JetPy[kMaxJets];
      float JetPz[kMaxJets];
      float JetPt[kMaxJets];
      int   NJets;

      float SumPtJets;
      float TriJetSumPt[kMaxCombinations];
      float TriJetMasses[kMaxCombinations];

      float MET;
    };

    SimpleEvent fEvt;
};
