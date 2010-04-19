#ifndef GUARD_Dtuple_h
#define GUARD_Dtuple_h

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TVector2.h"


#include "DHidasAna/Dtuple/interface/TLepton.h"
#include "DHidasAna/Dtuple/interface/TJet.h"
#include "DHidasAna/Dtuple/interface/TPhoton.h"
#include "DHidasAna/Dtuple/interface/TDUtility.h"


class TDtuple : public TDUtility
{

  public:
    TDtuple ();
    TDtuple (std::string const);
    TDtuple (std::string const, std::string const);
    TDtuple (TFile*);
    TDtuple (TChain*);
    virtual ~TDtuple ();

    void SetMaxTreeSize (int const);
    void SetBranches ();
    void SetBranchAddresses ();
    void Fill ();
    void Write ();
    int  GetEntry (int const);

  private:
    int   fRun;              // Run number
    int   fEvent;            // Event number
    int   fRunSection;       // Run section
    int   fEventFlags;       // Bits for different event flags
    int   fTriggerBits;      // Bits for triggers that passed
    float fMetX;             // X-component of Missing Et
    float fMetY;             // Y-component of Missing Et
    float fRawMetX;          // X-component of Raw Missing Et
    float fRawMetY;          // Y-component of Raw Missing Et
    float fSumEt;            // Sum Et
    float fRawSumEt;         // Raw Sum Et
    float fLum;              // Instantaneous Luminosity

  public:
    void  AddLepton (TLepton const&);
    void  AddLeptons (std::vector<TLepton> const&);
    void  AddLeptons (std::vector<TLepton>::iterator, std::vector<TLepton>::iterator);
    void  AddJet (TJet const&);
    void  AddJets (std::vector<TJet> const&);
    void  AddJets (std::vector<TJet>::iterator, std::vector<TJet>::iterator);
    void  AddPhoton (TPhoton const&);
    void  AddPhotons (std::vector<TPhoton> const&);
    void  AddPhotons (std::vector<TPhoton>::iterator, std::vector<TPhoton>::iterator);
    void  Clear ();
    void  SetRun (int const);
    void  SetEvent (int const);
    void  SetRunSection (int const);
    void  SetMetXY (float const, float const);
    void  SetMetX (float const);
    void  SetMetY (float const);
    void  SetRawMetX (float const);
    void  SetRawMetY (float const);
    void  SetRawMetXY (float const, float const);
    void  SetSumEt (float const);
    void  SetRawSumEt (float const);
    void  SetLum (float const);
    void  SetTriggerBit (std::string const);
    void  SetEventFlag (std::string const);


    int   GetRun ();
    int   GetEvent ();
    int   GetRunSection ();
    float GetMetX ();
    float GetMetY ();
    float GetMet ();
    float GetMetPhi ();
    float GetRawMetX ();
    float GetRawMetY ();
    float GetRawMet ();
    float GetSumEt ();
    float GetRawSumEt ();
    float GetLum ();
    bool  GetTriggerBit (std::string const);
    bool  GetEventFlag (std::string const);
    int   GetEventFlags ();

    std::vector<TLepton>*  GetLeptons ();
    std::vector<TJet>*     GetJets ();
    std::vector<TPhoton>*  GetPhotons ();

    void SortLeptons (std::vector<TLepton>::iterator, std::vector<TLepton>::iterator);
    void SortJets (std::vector<TJet>::iterator, std::vector<TJet>::iterator);
    void SortPhotons (std::vector<TPhoton>::iterator, std::vector<TPhoton>::iterator);

    TTree* GetTree ();
    void  DefaultValues ();


  private:
    void  SetRootFile (std::string const);
    void  SetRootFile (TFile*);

    void AddLeptonsToArray ();


    TFile *fRootFile;  // root file for output
    TTree *fDtupleTree;    // dtuple tree


    TClonesArray *fLepton;      // TLepton object array
    TClonesArray *fJet;         // TJet object array
    TClonesArray *fPhoton;         // TJet object array



  protected:
    std::vector<TLepton> Leptons;     // lepton vector filled in GetEntry
    std::vector<TJet> Jets;           // jet vector filled in GetEntry
    std::vector<TPhoton> Photons;     // photon vector filled in GetEntry

  public:
    // Root Likes ClassDef and ClassImp.
    // Comment them out if you don't need them.
    // There should NOT be a ; since this is a root macro and not a function
    // ClassDef must be the last line of the class before the };
    ClassDef(TDtuple,3) // Dtuple class
};


#endif
