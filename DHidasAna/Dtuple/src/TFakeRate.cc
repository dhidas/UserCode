#include "DHidasAna/Dtuple/interface/TFakeRate.h"

#include <iostream>
#include <map>

TFakeRate::TFakeRate (TString const FileName)
{
  std::cout << "TFakeRate: Reading fake file: " << FileName << std::endl;
  fFakeFile = new TFile(FileName, "read");
  if (!fFakeFile) {
    std::cerr << "ERROR: TFakeRate::TFakeRate() cannot open fake file." << std::endl;
  }

  ReadFakeHists();
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
  FakeMap1D[ std::make_pair<int, int>(0, 0) ] = (TH1D*) fFakeFile->Get("FakeRateJ2E");
  //FakeMap1D[ std::make_pair<int, int>(0, 0) ] = (TH1D*) fFakeFile->Get("FakeRateE");
  return;
}


float TFakeRate::GetFakeRate(int const Flavor, int const Type, float const Pt)
{
  TH1D* Hist = FakeMap1D[ std::make_pair<int, int>(Flavor, Type) ];
  if (Hist == 0x0) {
    std::cerr << "ERROR: TFakeRate::GetFakeRate does not see 1D hist" << std::endl;
    exit(1);
  }
  int Bin = Hist->FindBin(Pt);
  if (Bin == Hist->GetNbinsX() + 1) {
    --Bin;
  } else if (Bin == 0) {
    Bin = 1;
  }
  return Hist->GetBinContent(Bin);
}


float TFakeRate::GetFakeRate(int const Flavor, int const Type, float const Pt, float const Eta)
{
  TH2D* Hist = FakeMap2D[ std::make_pair<int, int>(Flavor, Type) ];
  int Bin = Hist->FindBin(Pt, Eta);
  //if (Bin == Hist->GetNBins() + 1) {
  //  --Bin;
  //}
  return Hist->GetBinContent(Bin);
}
