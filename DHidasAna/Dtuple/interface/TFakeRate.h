
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TH2D.h"




class TFakeRate
{
  public:
    TFakeRate (TString const);
    ~TFakeRate ();

    void ReadFakeHists ();
    float GetFakeRate (int const, int const, float const, float const);

  private:
    TFile* fFakeFile;
    std::map< std::pair<int, int>, TH2D*> FakeMap2D;
};
