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
    void PlotLeptons ();
    void PlotPhotons ();
    void PlotJets ();
    void PlotMonoLeptons ();
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
    bool PassSelectionElectron (TLepton&);
    bool PassSelectionMuon (TLepton&);
    bool IsDenominatorObject (TLepton&);
    bool IsDenominatorObject (TJet&);
    void PlotElectronId ();
    void SelectionPhoton ();
    void SelectionJet ();

    void PlotFakes ();
    void GetElectronFromJet (TJet&, TLepton&);
    void AddFakesToEvent (int const, int const NFakesToAdd = 1);
    void ElectronJetTest ();

    std::vector<TLepton> ClosestZMatch (std::vector<TLepton>&, bool const RequireOS = true, bool const RequireSF = true, bool const AddOtherLeptonsAtEnd = false);



  private:
    TString fProcName;
    TFile* fOutFile;
    std::map<int, std::pair<int, int> > fPlotTriLeptons_ElectronGenPMap;
    std::map< std::pair<int, int>, int> fPlotTriLeptons_ElectronGenPMotherMap;
    std::map<TString, int> fPlotTriLeptons_Counter;
    bool fRunFakes;
    TFakeRate* fFakeRate;
    TDtuple* fFakeDtuple;

};







#endif
