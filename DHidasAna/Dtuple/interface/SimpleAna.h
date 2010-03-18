#ifndef GUARD_SimpleAna_h
#define GUARD_SimpleAna_h

#include "DHidasAna/Dtuple/interface/TDtupleReader.h"

class SimpleAna : public TDtupleReader
{
  public:
    SimpleAna (TString const, TChain*);
    ~SimpleAna ();

    void Analyze (long unsigned int const);
    void PlotLeptons ();
    void PlotPhotons ();
    void PlotJets ();
    void PlotTriLeptons ();
    void PlotZllE ();
    void PlotDileptonMass ();
    void Plot6JetEvents ();
    void PlotNLepSumEtJetsMet ();

    TGenP* FindClosestGenP (TLepton&, int const Type = 0);

    void Selection ();
    void SelectionLepton ();
    void SelectionPhoton ();
    void SelectionJet ();

    std::vector<TLepton> ClosestZMatch (std::vector<TLepton>&, bool const RequireOS = true, bool const RequireSF = true, bool const AddOtherLeptonsAtEnd = false);



  private:
    TString fProcName;
    TFile* fOutFile;

};







#endif