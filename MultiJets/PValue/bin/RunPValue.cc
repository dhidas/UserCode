////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Fri Mar 11 17:41:48 CET 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <map>

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMinuit.h"






TH1F* GetPEExpo (int const N, float const cexp, TString const label = "blah")
{
  TH1F* h = new TH1F("PE", "PE", 65, 170, 800);

  TF1 f("myexp", "[0]*TMath::Exp([1]*x)", 170, 800);
  f.SetParameter(0, 1);
  f.SetParameter(1, cexp);
  f.SetParameter(0, 1.0/f.Integral(170,800));
  std::cout << cexp << "  " << f.Integral(170,800) << std::endl;

  for (int i = 0; i < N; ++i) {
    h->Fill(f.GetRandom());
  }

  TF1 ff("myexp1", "[0]*TMath::Exp([1]*x)", 170, 800);
  ff.SetParameter(0, 1);
  ff.SetParameter(1, cexp);
  ff.SetParameter(0, 10*N/ff.Integral(170,800));
  std::cout << N << "  " << ff.Integral(170,800) << std::endl;

  if (false) {
    TCanvas c;
    c.cd();
    h->Draw();
    ff.Draw("same");
    c.SaveAs(TString("MyPE_")+label+".eps");
  }

  return h;
}





float GetAcceptanceForMjjj (float const Mjjj)
{
  // Get the acceptance for a given mass

  // This is equivalent to the above, but a hell of a lot faster
  return -5.56179969055969892e-03 - 4.01623842089755741e-06 * Mjjj + 2.49580149009780901e-07 * Mjjj * Mjjj;

}



std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  // This will eventually be a parametrization

  if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  if (Mjjj < 350) return std::make_pair<float, float>(10 + 10 * (Mjjj - 250.) / 100., 15 + 10 * (Mjjj - 250.) / 100.);
  return std::make_pair<float, float>(20, 25);
}


//TMinuit* MyMinuit;
TH1F* hToFit;
TF1* fToFit;

void NegativeLogLikelihood2D (int& NParameters, double* gin, double& f, double* Par, int iflag)
{
  // Function we want to minimize!


  // Number of bins
  int const NBinsX = hToFit->GetNbinsX();

  double LogLikelihood = 0.0;

  // Loop over all bins in histogram
  for (int ibinX=1; ibinX <= NBinsX; ++ibinX) {

      double mu = 0.0;

      //std::cout << "Likelihood: " << ibinX << " " << ibinY << " " << -LogLikelihood << std::endl;
  }

  f = 0;

  return;
}








float DoFit (int const Section, int const ipe, float const SignalMass, TH1F* hPE, TF1* fFunc)
{

  std::pair<float, float> SignalMassRange = GetGausWidthRange(SignalMass);

  TF1 fSigBG("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + expo(3)",  170, 800);
  fSigBG.SetParLimits(0, -1000, 1000);
  fSigBG.SetParLimits(1, SignalMass, SignalMass);
  fSigBG.SetParLimits(2, SignalMassRange.first, SignalMassRange.second);
  //fSigBG.SetParLimits(2, SignalMassRange.first, SignalMassRange.first);
  fSigBG.SetParLimits(3, fFunc->GetParameter(3), fFunc->GetParameter(3));
  fSigBG.SetParLimits(4, fFunc->GetParameter(4), fFunc->GetParameter(4));
  fSigBG.SetParameter(0, 0);
  fSigBG.SetParameter(1, SignalMass);
  fSigBG.SetParameter(2, (SignalMassRange.first + SignalMassRange.second)/2);
  //fSigBG.SetParameter(2, SignalMassRange.first);
  fSigBG.SetParameter(3, fFunc->GetParameter(3));
  fSigBG.SetParameter(4, fFunc->GetParameter(4));

  hPE->Fit(&fSigBG, "LPQ", "", 170, 800);

  // return -9999 if the fit did not converge
  if (!gMinuit->fCstatu.Contains("CONVERGED")) {
    return -9999;
  }



  if (Section == -1) {
    char BUFF[200];
    sprintf(BUFF, "Fit_Data_%i.eps", (int) SignalMass);
    TCanvas c;
    c.cd();
    hPE->Draw();
    fSigBG.Draw("same");
    c.SaveAs(BUFF);

    printf("Gaus Params: %8i %12.3f %12.3f\n", (int) SignalMass, fSigBG.GetParameter(0), fSigBG.GetParameter(2));
  }


  float const CrossSection = (fSigBG.GetParameter(0) / 10.) * GetAcceptanceForMjjj(SignalMass) * 35.1;
  return CrossSection;
}




















int RunPValue (TString const InFileName, int const Section)
{
  // Some basic parameters
  float const StepSize    =  10;
  float const BeginMass   = 200;
  float const EndMass     = 500;
  int   const NPerSection = 1000;


  // Set the randome seed based on section number
  gRandom->SetSeed(12332 * (Section + 2));


  // Setup output file
  char OutName[150];
  if (Section == -1) {
    sprintf(OutName, "Data_CrossSection.dat");
  } else {
    sprintf(OutName, "CrossSection_%i.dat", Section);
  }
  FILE* Out = fopen(OutName, "w");
  if (Out == NULL) {
    std::cerr << "ERROR: cannot open output file: " << OutName << std::endl;
    exit(1);
  }
  for (float ThisMass = BeginMass; ThisMass <= EndMass; ThisMass += StepSize) {
    fprintf(Out, "%10.3f ", ThisMass);
  }
  fprintf(Out, "\n");




  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1F* DataTH1F = (TH1F*) InFile.Get("Mjjj_45_20_130_6jet");
  TF1* fFunc = DataTH1F->GetFunction("total");

  int const NFromData = DataTH1F->Integral(17, 80);


  if (Section == -1) {
    // Do some Data Fit

    // In order to do this correctly you need a clone...ya I know..
    TH1F* DataClone = (TH1F*) DataTH1F->Clone("DataClone");

    // Loop over all masses
    for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {
      float const CrossSection = DoFit(Section, 0, SignalMass, DataClone, fFunc);
      printf("ipe: %12i Mass: %5i  xs:%15.2f\n", 0, (int) SignalMass, CrossSection);
      fprintf(Out, "%12E ", CrossSection);
    }
    fprintf(Out, "\n");
    fflush(Out);

    delete DataClone;
  } else {

    // Loop over all PEs
    for (int ipe = 0; ipe != NPerSection; ++ipe) {

      // Number in PE
      int const NPE = gRandom->Poisson(NFromData);

      // Get PE
      TH1F* hPE = GetPEExpo(NPE, fFunc->GetParameter(1));

      // Loop over signal masses
      for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {
        float const CrossSection = DoFit(Section, ipe, SignalMass, hPE, fFunc);
        printf("ipe: %12i Mass: %5i  xs:%15.2f\n", ipe, (int) SignalMass, CrossSection);
        fprintf(Out, "%12E ", CrossSection);
      }
      fprintf(Out, "\n");
      fflush(Out);

      // Better delete that PE
      delete hPE;
    }


  }



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [Section]" << std::endl;
    return 1;
  }

  //TString const InFileName = "/Users/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  TString const InFileName = "/users/h2/dhidas/Data35pb/ExpFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  int const Section = atoi(argv[1]);

  RunPValue(InFileName, Section);

  return 0;
}
