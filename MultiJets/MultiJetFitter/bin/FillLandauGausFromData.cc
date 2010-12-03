////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Oct 12 12:18:33 EDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"



int FillLandauGaus (TString const InFileName)
{
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open Input file" << std::endl;
    return -1;
  }

  TFile f("DeanDataAndPE_3.root", "recreate");
  if (!f.IsOpen()) {
    std::cerr << "ERROR: cannot open output file" << std::endl;
    return -1;
  }

  TRandom r(143);




  std::vector<TString> HistNames;
  HistNames.push_back("20_20_200");
  HistNames.push_back("45_20_120");


  for (size_t iDiag = 0; iDiag != HistNames.size(); ++iDiag) {
    TH1* Hist = (TH1*) InFile.Get("Mjjj_"+HistNames[iDiag]);
    if (!Hist) {
      std::cerr << "ERROR: cannot get base hist: " << HistNames[iDiag] << std::endl;
      exit(1);
    }

    // Copy to new file
    TH1F* HistClone = (TH1F*) Hist->Clone();
    HistClone->SetDirectory(&f);
    HistClone->Write();

    TF1 Land("land", "[0]*TMath::Landau(x, [1], [2], 1)", Hist->GetXaxis()->GetXmin(), Hist->GetXaxis()->GetXmax());
    Land.SetParameter(0, Hist->GetEntries());
    Land.SetParameter(1, Hist->GetMean());
    Land.SetParameter(2, Hist->GetRMS());
    Hist->Fit("land");

    float const nland   = Land.GetParameter(0);
    float const enland  = Land.GetParError(0);
    float const lmpv    = Land.GetParameter(1);
    float const elmpv   = Land.GetParError(1);
    float const lsigma  = Land.GetParameter(2);
    float const elsigma = Land.GetParError(2);

    int   const NData = (int) Hist->Integral();

    int   const NPE = 10000;

    char BUFF[200];
    float ThisLMPV;
    float ThisLSigma;
    bool const DoSyst = false;
    for (int ipe = 0; ipe != NPE; ++ipe) {
      if (ipe % 1000 == 0) {
        std::cout << HistNames[iDiag] << " Filling PE: " << ipe << std::endl;
      }
      sprintf(BUFF, "PE_%s_%i", HistNames[iDiag].Data(), ipe);
      TH1F h(BUFF, BUFF, Hist->GetNbinsX(), Hist->GetXaxis()->GetXmin(), Hist->GetXaxis()->GetXmax());

      ThisLMPV   = DoSyst ?    lmpv + r.Gaus(0, elmpv) : lmpv;
      ThisLSigma = DoSyst ? lsigma + r.Gaus(0, elsigma) : lsigma;

      int NDataThisPE = NData;//r.Poisson(NData);
      for (int i = 0; i != NDataThisPE; ++i) {
        h.Fill( r.Landau(ThisLMPV, ThisLSigma) );
      }
      h.Write();
    }

  }

  f.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFile]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];

  return FillLandauGaus(InFileName);
}
