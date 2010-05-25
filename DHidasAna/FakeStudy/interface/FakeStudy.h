#ifndef GUARD_FakeStudy_h
#define GUARD_FakeStudy_h

#include "DHidasAna/Dtuple/interface/TDtupleReader.h"
#include "DHidasAna/Dtuple/interface/TFakeRate.h"

class FakeStudy : public TDtupleReader
{
  public:
    FakeStudy (TString const, TChain*);
    ~FakeStudy ();

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




  private:
    TString fProcName;
    TFile* fOutFile;
    std::map<int, std::pair<int, int> > fPlotTriLeptons_ElectronGenPMap;
    std::map< std::pair<int, int>, int> fPlotTriLeptons_ElectronGenPMotherMap;
    std::map<TString, int> fPlotTriLeptons_Counter;
    std::map< std::pair<int, int>, std::pair<int, int> > fPlotLeptons_ElectronGenPMotherMap;
    int fPlotLeptons_Electrons;
    int fPlotLeptons_ElectronsAfterConvVeto;
    bool fRunFakes;
    TFakeRate* fFakeRate;
    TDtuple* fFakeDtuple;

};







#endif
