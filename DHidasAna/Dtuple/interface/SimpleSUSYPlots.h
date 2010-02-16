#ifndef GUARD_SimpleSUSYPlots_h
#define GUARD_SimpleSUSYPlots_h

#include "DHidasAna/Dtuple/interface/DtupleReader.h"

#include <map>

#include "TH1D.h"


class SimpleSUSYPlots : public DtupleReader
{
  public:
    SimpleSUSYPlots (TString const&);
    SimpleSUSYPlots (std::vector<TString> const&);
    ~SimpleSUSYPlots ();

    void InitOutFile (TString const& OutFileName = "Output.root");

    void InitializeHists ();
    void Loop ();

    void SortOutOverlaps (DtupleReader::Event_Struct&);
    void CopyILeptonFromTo (int const, Dtuple::Event_Struct&, Dtuple::Event_Struct&);

  private:
    TFile* fOutFile;
    std::map<TString, TH1D*> Hist1D;

};





















#endif
