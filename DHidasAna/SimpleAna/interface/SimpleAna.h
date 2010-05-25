#ifndef GUARD_SimpleAna_h
#define GUARD_SimpleAna_h

#include "DHidasAna/Dtuple/interface/TDtupleReader.h"
#include "DHidasAna/Dtuple/interface/TFakeRate.h"

class SimpleAna : public TDtupleReader
{
  public:
    SimpleAna (TString const, TChain*);
    ~SimpleAna ();

    void Analyze (long unsigned int const);
    void BeginJob ();
    void EndJob ();

    void SetFakeRateFile(TString const);
    void RunFakes (bool const);
    bool RunFakes ();

    void PlotEventQuantities ();
    void PlotElectronId ();
    void PlotLeptons ();
    void PlotPhotons ();
    void PlotJets ();

    TGenP* FindClosestGenP (TLepton&, int const Type = 0);

    void Selection ();
    void SelectionLepton ();
    bool PassSelectionElectron (TLepton&);
    bool PassSelectionMuon (TLepton&);
    void SelectionPhoton ();
    void SelectionJet ();

    void GetElectronFromJet (TJet&, TLepton&);




  private:
    TString fProcName;
    TFile* fOutFile;
    bool fRunFakes;
    TFakeRate* fFakeRate;
    TDtuple* fFakeDtuple;

};







#endif
