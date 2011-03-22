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



// I'm using these as globals just for ease.
TH1F* hToFit;  // Histogram to fit
TF1* fToFit;   // Function used in fitting
float ParC[5]; // Parameter central value
float ParE[5]; // Parameter error




TH1F* GetPEExpo (float const inN, float const incexp, int const ipe,  bool const DoSyst)
{
  // This function gets a PE in any way you tell it to!  ha!

  // Systematics in PEs or not??
  float const cexp = DoSyst ? incexp * (1. + 0.05 * gRandom->Gaus(0, 1)) : incexp;
  int   const N    = DoSyst ? gRandom->Poisson(inN * (1. + 0.03 * gRandom->Gaus(0, 1))) : gRandom->Poisson(inN);
  printf("PE cexp/N: %12E  %9i\n", cexp, N);


  // New histogram we'll return later on..  you will own this, not me!!
  TH1F* h = new TH1F("PE", "PE", 65, 170, 800);

  // Define normalized BG function in the range
  TF1 f("myexp", "[0]*TMath::Exp([1]*x)", 170, 800);
  f.SetParameter(0, 1);
  f.SetParameter(1, cexp);
  f.SetParameter(0, 1.0/f.Integral(170,800));

  // Fill the new histogram based on the BG pdf defined above
  for (int i = 0; i < N; ++i) {
    h->Fill(f.GetRandom());
  }

  // This is a test so that I know that I have the normalizaton correct
  // In general, don't run this..
  if (true) {
    // Define the function..
    TF1 ff("myexp1", "[0]*TMath::Exp([1]*x)", 170, 800);
    ff.SetParameter(0, 1);
    ff.SetParameter(1, cexp);
    ff.SetParameter(0, 10*N/ff.Integral(170,800));

    // Plot it and the PE
    TCanvas c;
    c.cd();
    h->Draw("ep");
    ff.Draw("same");
    char BUFF[100];
    sprintf(BUFF, "PE_%i.eps", ipe);
    c.SaveAs(BUFF);
  }

  // YOU own h, not me!!  remember to delete it
  return h;
}





float GetAcceptanceForMjjj (float const Mjjj)
{
  // Get the acceptance for a given mass.  Thanks Dan!

  // Hell of a lot faster
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
  double LogSum = 0.0;
  for (int i=1; i <= in; ++i) {
    Logs.push_back(TMath::Log(i));
  }

  // Are you really asking for the log factorial of such a large number?
  // Shame on you.
  LogSum = KahanSummation(Logs.begin(), Logs.end());

  // Fine, I'll give it to you anyway.
  return LogSum;
}





void NegativeLogLikelihood (int& NParameters, double* gin, double& f, double* Par, int iflag)
{
  // Function we want to minimize!  This function is given as is and no attention has been paid
  // to the normalization.  ie if you care about what the likelihood actually is, you'll have to
  // add the correct factors.  Makes no difference for minimization though...so to save
  // a few cycles they have been ignored.
  // Also, take care when using bins with a very large number of data points.  I've taken a bit
  // of care here using a Kahan Summation, but one should know what that is and how useful it
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

  // LL variable
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






float MinimizeNLL (int const Section, int const ipe, float const SignalMass, TF1* fFunc, bool const IsBGOnly, bool const DoSyst)
{
  // This function will setup Minuit and call the minimization routine.

  // Number of minuit parameters
  int const NParams = 5;

  // These are for Minuit
  double ArgList[10];
  int ErrorFlag;

  // Grab a new instance of Minuit with NParams
  TMinuit MyMinuit(NParams);

  // Set the printlevel -2 minimum, higher for more junk
  MyMinuit.SetPrintLevel(-2);

  // Get the Expected numbers from the input function and other params and limits
  // NOTE: You should always stay AWAY from the limits... ie the fit should NOT
  // hit the limit
  std::pair<float, float> gWidth = GetGausWidthRange(SignalMass);
  float const NBG = TMath::Exp(fFunc->GetParameter(0));
  float const NBG_min = 0;
  float const NBG_max = 1000;
  float const BGExp = fFunc->GetParameter(1);
  float const BGExp_min = -0.1;
  float const BGExp_max =  0;
  // Errors...
  float const eNBG = NBG * 0.03;
  float const eBGExp = fFunc->GetParError(1);


  // This is the function we're working with
  // TF1 ("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + [3]*TMath::Exp([4]*x)",  170, 800);
  if (IsBGOnly) {
    // This is for BG only fits
    ParC[0] = -999;
    ParC[1] = SignalMass;
    ParC[2] = (gWidth.first+gWidth.second)/2;
    ParC[3] = NBG;
    ParC[4] = BGExp;

    ParE[0] = -999;
    ParE[1] = -999;
    ParE[2] = -999;
    ParE[3] = eNBG;
    ParE[4] = eBGExp;
    if (!DoSyst) {
      ParE[3] = -999;
      ParE[4] = -999;

      int NPar = NParams;
      double Pars[5] = {ParC[0], ParC[1], ParC[2], ParC[3], ParC[4]};
      double *Blank;
      double NLL = 0;

      NegativeLogLikelihood(NPar, Blank, NLL, Pars, 0);
      return NLL;
    }

    MyMinuit.DefineParameter(0, "GausNorm", 0, 0.001, -1000, 1000);
    MyMinuit.DefineParameter(1, "GausMean",  SignalMass, 0.01, SignalMass-1, SignalMass+1);
    MyMinuit.DefineParameter(2, "GausWidth", (gWidth.first+gWidth.second)/2., 0.001, gWidth.first, gWidth.second);
    MyMinuit.DefineParameter(3, "BGNorm", NBG, 0.001, NBG_min, NBG_max);
    MyMinuit.DefineParameter(4, "BGExp", BGExp, 0.000001, BGExp_min, BGExp_max);
    MyMinuit.FixParameter(0);
    MyMinuit.FixParameter(1);
    MyMinuit.FixParameter(2);
    if (!DoSyst) {
      MyMinuit.FixParameter(3);
      MyMinuit.FixParameter(4);
    }
  } else {
    // This is for Signal+BG fits
    ParC[0] = -999;
    ParC[1] = SignalMass;
    ParC[2] = (gWidth.first+gWidth.second)/2;
    ParC[3] = NBG;
    ParC[4] = BGExp;

    ParE[0] = -999;
    ParE[1] = -999;
    ParE[2] = -999;
    ParE[3] = eNBG;
    ParE[4] = eBGExp;
    if (!DoSyst) {
      ParE[3] = -999;
      ParE[4] = -999;
    }

    MyMinuit.DefineParameter(0, "GausNorm", 0, 0.001, -1000, 1000);
    MyMinuit.DefineParameter(1, "GausMean",  SignalMass, 0.01, SignalMass-1, SignalMass+1);
    MyMinuit.DefineParameter(2, "GausWidth", (gWidth.first+gWidth.second)/2., 0.001, gWidth.first, gWidth.second);
    MyMinuit.DefineParameter(3, "BGNorm", NBG, 0.001, NBG_min, NBG_max);
    MyMinuit.DefineParameter(4, "BGExp", BGExp, 0.000001, BGExp_min, BGExp_max);
    MyMinuit.FixParameter(1);
    if (!DoSyst) {
      MyMinuit.FixParameter(3);
      MyMinuit.FixParameter(4);
    }

  }

  // Set error definition
  // 1 for Chi squared
  // 0.5 for nagative log likelihood
  MyMinuit.SetErrorDef(0.5);

  // Set Minimization strategy
  // 1 standard 
  // 2 try to improve minimum (slower) 
  ArgList[0] = 1;
  MyMinuit.mnexcm("SET STR", ArgList, 1, ErrorFlag);

  // Set the function to minimize and maximum number of iterations
  MyMinuit.SetFCN(NegativeLogLikelihood);
  MyMinuit.SetMaxIterations(1000);

  // Actual call to do minimimazation
  MyMinuit.Migrad();

  // Grab the status of the fit.  If the fit does not converge we better return
  // some positive value...
  TString FitStatus = MyMinuit.fCstatu;
  if (!FitStatus.Contains("CONVERGED")) {
    std::cerr << "ERROR: Fit did not converge" << std::endl;
    return 9999;
  }

  // Save a plot if you like!
  if (Section == -1 || false) {
    char BUFF[100];
    if (SignalMass > 0) {
      sprintf(BUFF, "Fit_Data_%i.eps", (int) SignalMass);
    } else {
      sprintf(BUFF, "Fit_Data_BG.eps");
    }
    TCanvas c;
    hToFit->SetAxisRange(170, 800, "X");
    hToFit->Draw("hist");
    fToFit->Draw("same");
    c.SaveAs(BUFF);
  }

  // Call for minos error analysis
  // nah, screw that
  // MyMinuit.mnmnos();

  // Grab parameters from minuit
  double OPar[NParams], OEPar[NParams];
  for (int i = 0; i < NParams; ++i) {
    MyMinuit.GetParameter(i, OPar[i], OEPar[i]);
    printf("Fitted Param: %i %9E  +/-  %9E\n", i, OPar[i], OEPar[i]);
  }


  // For the parameters which minimize the NLL, lets get the NLL value and return it
  double NLL = 0;
  double* Blank;
  int NPar = NParams;
  NegativeLogLikelihood(NPar, Blank, NLL, OPar, 0);

  return NLL;
}




















int RunPValue (TString const InFileName, int const Section, bool const DoSyst)
{
  // Some basic parameters
  float const StepSize    =  10;
  float const BeginMass   = 200;
  float const EndMass     = 500;
  int   const NPerSection =  50;


  // Set the randome seed based on section number
  gRandom->SetSeed(943437 * (Section + 2));


  // Setup output file
  char OutName[150];
  if (Section == -1) {
    sprintf(OutName, "Data_TestStatistic.dat");
  } else {
    sprintf(OutName, "TestStatistic_%i.dat", Section);
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




  // Get the data hist and fit function
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1F* DataTH1F = (TH1F*) InFile.Get("Mjjj_45_20_130_6jet");
  TF1* fFunc = DataTH1F->GetFunction("total");

  // Number of events from data
  int const NFromData = (int) DataTH1F->Integral(17, 80);

  // Set the functino to fit.  Be VERY careful if oyu change this... you need to change a lot
  // of other things too...
  fToFit = new TF1("FitFunction", "[0]*TMath::Gaus(x, [1], [2], 1) + [3]*TMath::Exp([4]*x)", 170, 800);

  if (Section == -1) {
    // Do some Data Fit

    // In order to do this correctly you need a clone...ya I know..
    TH1F* DataClone = (TH1F*) DataTH1F->Clone("DataClone");
    hToFit = DataClone;

    // Grab the background only LL
    float const BGLL = -1.0 * MinimizeNLL(Section, -1, -999, fFunc, true, DoSyst);
    if (BGLL == -9999) {
      std::cerr << "ERROR: bad BG only fit for data.  This is a real problem..." << std::endl;
      throw;
    }

    // Loop over all masses
    for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {

      // Get LL for this mass..
      float const MyLL = -1.0 * MinimizeNLL(Section, -1, SignalMass, fFunc, false, DoSyst);

      // Get the test statistic
      float const TestStatistic = MyLL != -9999 ? -2.0 * (BGLL - MyLL) : -9999.;

      // Print out and file!
      printf("ipe: %12i Mass:%5i  BGLL:%12E  MyLL:%12E  D:%12E\n", -1, (int) SignalMass, BGLL, MyLL, TestStatistic);
        fprintf(Out, "%12E ", TestStatistic);
    }
    fprintf(Out, "\n");
    fflush(Out);

    delete DataClone;
  } else {

    // Loop over all PEs
    for (int ipe = 0; ipe != NPerSection; ++ipe) {

      // Number in PE
      //float const NPE = (float) NFromData;
      float const NPE = fFunc->Integral(170, 800) / 10.;
      printf("iep: %i NPE: %12.2f\n", ipe, NPE);

      // Get PE
      hToFit = GetPEExpo(NPE, fFunc->GetParameter(1), ipe, DoSyst);

      // Grab the background only LL
      float const BGLL = -1.0 * MinimizeNLL(Section, ipe, -999, fFunc, true, DoSyst);
      if (BGLL == -9999) {
        std::cerr << "skipping this PE... it's bad!" << std::endl;
        --ipe;
        continue;
      }

      // Loop over signal masses
      for (float SignalMass = BeginMass; SignalMass <= EndMass; SignalMass += StepSize) {

        // Get LL for this mass..
        float const MyLL = -1.0 * MinimizeNLL(Section, ipe, SignalMass, fFunc, false, DoSyst);

        // Get the test statistic
        float const TestStatistic = MyLL != -9999 ? -2.0 * (BGLL - MyLL) : -9999.;

      // Print out and file!
        printf("ipe: %12i Mass:%5i  BGLL:%12E  MyLL:%12E  D:%12E\n", ipe, (int) SignalMass, BGLL, MyLL, TestStatistic);
        fprintf(Out, "%12E ", TestStatistic);
      }
      fprintf(Out, "\n");
      fflush(Out);

      // Better delete that PE
      delete hToFit;
    }


  }

  delete fToFit;



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [Section]" << std::endl;
    return 1;
  }

  TString const InFileName = "/uscms/home/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  //TString const InFileName = "/Users/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  //TString const InFileName = "/users/h2/dhidas/Data35pb/ExpFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  int const Section = atoi(argv[1]);

  // Run Systematics?
  bool const DoSyst = false;

  RunPValue(InFileName, Section, DoSyst);

  return 0;
}
