#ifndef GUARD_SimpleAna_h
#define GUARD_SimpleAna_h

#include "DHidasAna/Dtuple/interface/TDtupleReader.h"

class SimpleAna : public TDtupleReader
{
  public:
    SimpleAna (TString const, TChain*);
    ~SimpleAna ();

    void Analyze (long unsigned int const);
    void BeginJob ();
    void EndJob ();

    void PlotEventQuantities ();
    void PlotLeptons ();
    void PlotPhotons ();
    void PlotJets ();
    void PlotTriLeptons ();
    void PlotlGamma ();
    void PlotllGamma ();
    void PlotZllE ();
    void PlotDileptonMass ();
    void Plot6JetEvents ();
    void PlotNLepSumEtJetsMet ();

    TGenP* FindClosestGenP (TLepton&, int const Type = 0);

    void Selection ();
    void SelectionLepton ();
    void SelectionPhoton ();
    void SelectionJet ();

    void ElectronJetTest ();

    std::vector<TLepton> ClosestZMatch (std::vector<TLepton>&, bool const RequireOS = true, bool const RequireSF = true, bool const AddOtherLeptonsAtEnd = false);



  private:
    TString fProcName;
    TFile* fOutFile;
    std::map<int, int> fPlotTriLeptons_ElectronGenPMap;

};







#endif
