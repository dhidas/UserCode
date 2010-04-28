#include "DHidasAna/Dtuple/interface/TFakeRate.h"

#include <iostream>
#include <map>

TFakeRate::TFakeRate (TString const FileName)
{
  fFakeFile = new TFile(FileName, "read");
  if (!fFakeFile) {
    std::cerr << "ERROR: TFakeRate::TFakeRate() cannot open fake file." << std::endl;
  }
}


TFakeRate::~TFakeRate ()
{
  if (fFakeFile) {
    fFakeFile->Close();
  }
}


void TFakeRate::ReadFakeHists ()
{
  FakeMap2D[ std::make_pair<int, int>(0, 0) ] = (TH2D*) fFakeFile->Get("FakeRate");
  return;
}


float TFakeRate::GetFakeRate(int const Flavor, int const Type, float const Pt, float const Eta)
{
  TH2D* Hist = FakeMap2D[ std::make_pair<int, int>(Flavor, Type) ];
  int const Bin = Hist->FindBin(Pt, Eta);
  return Hist->GetBinContent(Bin);
}
