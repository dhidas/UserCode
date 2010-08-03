#include "VgammaAna/TaiwanAna/interface/VgNtuple.h"

#include "TFile.h"

class VgAna : public VgNtuple
{
  public:
    VgAna ();
    ~VgAna ();

    void Loop ();
    void SetOutFile (TString const&);

    void PlotElectrons ();
    void PlotMETs ();
    void PlotJets ();
    void PlotElectronPhoton ();

    TFile* OutFile_;
};
