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
#include "TMath.h"






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

template <class T> double KahanSummation(T begin, T end)
{
  // You should be careful how you sum a lot of things...

  double result = 0.f;

  double c = 0.f;
  double y, t;
  for ( ; begin != end; ++begin) {
    y = *begin - c;
    t = result + y;
    c = (t - result) - y;
    result = t;
  }
  return result;
}





long double LogFactorial (int const in)
{
  // Get the log of the factorial.
  // The reason for this function is that we're computationally limited to
  // i <= 170 for the machine (and TMath::Factorial)

  std::vector<float> Logs;
  Logs.reserve(in);
  long double LogSum = 0.0;
  for (int i=1; i <= in; ++i) {
    Logs.push_back(TMath::Log(i));
  }

  // Are you really asking for the log factorial of such a large number?
  // Shame on you.
  LogSum = KahanSummation(Logs.begin(), Logs.end());

  // Fine, I'll give it to you anyway.
  return LogSum;
}





// I'm using these as globals just for ease.
TMinuit* MyMinuit;
TH1F* hToFit;
TF1* fToFit;
float ParC[5];
float ParE[5];

void NegativeLogLikelihood (int& NParameters, double* gin, double& f, double* Par, int iflag)
{
  // Function we want to minimize!  This function is given as is and no attention has been paid
  // to the normalization.  ie if you care about what the likelihood actually is, you'll have to
  // add the correct factors.  Makes no difference for minimization though...so to save
  // a few cycles they have been ignored.
  // Also, take care when using bins with a very large number of data points.  I've taken a bit
  // of care here using a Kalman Summation, but one should know what that is and how useful it
  // is if one should go beyond this

  // You should not print anything in this function save for debug time
  //for (int i = 0; i != 5; ++i) {
  //  printf("Par: %i %7E %7E %7E\n", i, Par[i], ParC[i], ParE[i]);
  //}


  // These are the parameters of the fit...
  // ("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + [3]*TMath::Exp([4]*x)",  170, 800);
  // Systematics: start on the 5.  should be gaussian centered at 0 width one
  fToFit->SetParameter(0, Par[0]);
  fToFit->SetParameter(1, Par[1]);
  fToFit->SetParameter(2, Par[2]);
  fToFit->SetParameter(3, Par[3]);
  fToFit->SetParameter(4, Par[4]);

  // Number of bins
  static int const NBinsX = hToFit->GetNbinsX();

  static long double LogLikelihood;
  LogLikelihood = 0.0;
  static long double mu = 0.0;

  // Loop over all bins in histogram
  static int const iStart = hToFit->FindBin(170);
  static int const iStop  = hToFit->FindBin(800);
  for (int ibinX = iStart; ibinX <= iStop; ++ibinX) {

    mu = fToFit->Integral( hToFit->GetBinLowEdge(ibinX), hToFit->GetBinLowEdge(ibinX) + hToFit->GetBinWidth(ibinX)) / hToFit->GetBinWidth(ibinX);
    //printf("Bin X mu: %4i %9.1f %9E\n", ibinX, hToFit->GetBinLowEdge(ibinX), mu);

    if (mu > 0) {
      LogLikelihood += (hToFit->GetBinContent(ibinX) * TMath::Log(mu)
          - mu - LogFactorial((int) hToFit->GetBinContent(ibinX)  ) );
    }
  }

  for (int i = 0; i < 5; ++i) {
    if (ParE[i] != -999) {
      LogLikelihood -= TMath::Power((Par[i] - ParC[i])/ParE[i], 2);
    }
  }


  f = -1.0*LogLikelihood;
  //printf("%15E\n", f);

  return;
}






float MinimizeNLL (int const Section, int const ipe, float const SignalMass, TF1* fFunc)
{
  int const N = 5;
  double ArgList[10];
  int ErrorFlag;

  ArgList[0] = 1;

  MyMinuit = new TMinuit(N);
  MyMinuit->SetPrintLevel(-2);


  std::pair<float, float> gWidth = GetGausWidthRange(SignalMass);
  float const NBG = TMath::Exp(fFunc->GetParameter(0));
  float const NBG_min = 0;
  float const NBG_max = 1000;
  float const BGExp = fFunc->GetParameter(1);
  float const BGExp_min = -0.1;
  float const BGExp_max =  0;

  float const eNBG = NBG * 0.03;
  float const eBGExp = fFunc->GetParError(1);

  // TF1 fSigBG("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + [3]*TMath::Exp([4]*x)",  170, 800);
  MyMinuit->DefineParameter(0, "GausNorm", 0, 0.01, 0, 1000);
  MyMinuit->DefineParameter(1, "GausMean",  SignalMass, 0.01, SignalMass-1, SignalMass+1);
  MyMinuit->DefineParameter(2, "GausWidth", (gWidth.first+gWidth.second)/2., 0.01, gWidth.first, gWidth.second);
  MyMinuit->DefineParameter(3, "BGNorm", NBG, 0.01, NBG_min, NBG_max);
  MyMinuit->DefineParameter(4, "BGExp", BGExp, 0.0001, BGExp_min, BGExp_max);
  MyMinuit->FixParameter(1);
  //MyMinuit->FixParameter(3);
  //MyMinuit->FixParameter(4);

  ParC[0] = -999;
  ParC[1] = SignalMass;
  ParC[2] = (gWidth.first+gWidth.second)/2;
  ParC[3] = NBG;
  ParC[4] = BGExp;

  ParE[0] = -999;
  ParE[1] = -999;
  ParE[2] = -999;
  ParE[3] = eNBG;
  ParE[4] = eBGExp;;

  // Set error definition
  // 1 for Chi squared
  // 0.5 for nagative log likelihood
  MyMinuit->SetErrorDef(0.5);

  // Set Minimization strategy
  // 1 standard 
  // 2 try to improve minimum (slower) 
  ArgList[0]=2;
  MyMinuit->SetFCN(NegativeLogLikelihood);
  MyMinuit->mnexcm("SET STR", ArgList, 1, ErrorFlag);

  // Set the maximum number of iterations
  MyMinuit->SetMaxIterations(1000);

  // Actual call to do minimimazation
  MyMinuit->Migrad();

  TString FitStatus = MyMinuit->fCstatu;
  if (!FitStatus.Contains("CONVERGED")) {
    std::cerr << "WARNING: Fit did not converge" << std::endl;
  }
  std::cout << MyMinuit->fCstatu << " " << MyMinuit->GetStatus() << std::endl;
  if (Section == -1 || false) {
    char BUFF[100];
    sprintf(BUFF, "Fit_Data_%i.eps", (int) SignalMass);
    TCanvas c;
    hToFit->Draw("hist");
    fToFit->Draw("same");
    c.SaveAs(BUFF);
  }

  // Call for minos error analysis
  MyMinuit->mnmnos();

  // Grab parameters from minuit
  double OPar[N], OEPar[N];
  for (int i = 0; i < N; ++i) {
    MyMinuit->GetParameter(i, OPar[i], OEPar[i]);
    printf("Fitted Param: %i %9E  +/-  %9E\n", i, OPar[i], OEPar[i]);
  }


  printf("Mass / CrossSection: %4i  %9E\n", (int) SignalMass, (OPar[0] / 10.) * GetAcceptanceForMjjj(SignalMass) * 35.1);
  return (OPar[0] / 10.) * GetAcceptanceForMjjj(SignalMass) * 35.1;
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

  fToFit = new TF1("FitFunction", "[0]*TMath::Gaus(x, [1], [2], 1) + [3]*TMath::Exp([4]*x)", 170, 800);

  if (Section == -1) {
    // Do some Data Fit

    // In order to do this correctly you need a clone...ya I know..
    TH1F* DataClone = (TH1F*) DataTH1F->Clone("DataClone");
    hToFit = DataClone;

    // Loop over all masses
    for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {
      //float const CrossSection = DoFit(Section, 0, SignalMass, DataClone, fFunc);
      float const CrossSection = MinimizeNLL(Section, 0, SignalMass, fFunc);
      //float const CrossSection = DoFitNLL(Section, 0, SignalMass, DataClone, fFunc);
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
      hToFit = hPE;

      // Loop over signal masses
      for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {
        //float const CrossSection = DoFit(Section, ipe, SignalMass, hPE, fFunc);
        float const CrossSection = MinimizeNLL(Section, ipe, SignalMass, fFunc);
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

  TString const InFileName = "/Users/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  //TString const InFileName = "/users/h2/dhidas/Data35pb/ExpFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  int const Section = atoi(argv[1]);

  RunPValue(InFileName, Section);

  return 0;
}
