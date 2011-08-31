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

int const NPAR = 6;
float ParC[NPAR]; // Parameter central value
float ParE[NPAR]; // Parameter error

float const XMIN     =  230;
float const XMAX     = 1520;
float const BINWIDTH = 10;
int   const NBINS    = (XMAX - XMIN) / BINWIDTH;



TH1F* GetPE (float const inN, std::vector< std::pair<float, float> >& InPar, int const ipe,  bool const DoSyst)
{
  // This function gets a PE in any way you tell it to!  ha!

  // Parameters for PE, with or without syst
  std::vector<float> Par;
  for (size_t i = 0; i != InPar.size(); ++i) {
    if (DoSyst) {
      Par.push_back( InPar[i].first * (1. + InPar[i].second * gRandom->Gaus(0, 1)) );
    } else {
      Par.push_back(InPar[i].first);
    }
  }


  // Systematics in PEs or not??
  int   const N    = DoSyst ? gRandom->Poisson(inN * (1. + 0.03 * gRandom->Gaus(0, 1))) : gRandom->Poisson(inN);

  // New histogram we'll return later on..  you will own this, not me!!
  TH1F* h = new TH1F("PE", "PE", NBINS, XMIN, XMAX);

  // Define normalized BG function in the range
  TF1 f("myexp", "[0]*(((1.-x/7000.)^[1])/((x/7000.)^[2]))", XMIN, XMAX);
  for (size_t i = 0; i != Par.size(); ++i) {
    f.SetParameter(i, Par[i]);
  }
  f.SetParameter(0, 1);
  f.SetParameter(0, 1.0 / f.Integral(XMIN, XMAX));

  // Fill the new histogram based on the BG pdf defined above
  for (int i = 0; i < N; ++i) {
    h->Fill(f.GetRandom());
  }

  // This is a test so that I know that I have the normalizaton correct
  // In general, don't run this..
  if (false) {
    // Define the function..
    TF1 ff("myexp1", "[0]*(((1.-x/7000.)^[1])/((x/7000.)^[2]))", XMIN, XMAX);
    for (size_t i = 0; i != Par.size(); ++i) {
      ff.SetParameter(i, Par[i]);
    }
    ff.SetParameter(0, 1);
    ff.SetParameter(0, BINWIDTH * N / ff.Integral(XMIN, XMAX));

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


float GetACCERROR (float const Mjjj)
{
  // 250    0.203786
  // 300    0.175059
  // 350    0.127826
  // 400    0.136186
  // 450    0.140467
  // 500    0.171287
  // 750    0.212378
  //1250    0.335977
  //1500    0.565487

  // Fixed 250-500 to the 250 val and fit pol3. might undershoot around 1000, but ain't too bad
  //p0                        =     0.140834   +/-   0.0160314   
  //p1                        =  0.000392742   +/-   7.3481e-05  
  //p2                        = -7.49412e-07   +/-   9.60659e-08 
  //p3                        =  4.50564e-10   +/-   3.69923e-11 

  return 0.140834 + 0.000392742*Mjjj - -7.49412e-07*Mjjj*Mjjj + 4.50564e-10*Mjjj*Mjjj*Mjjj;
}



float GetAcceptanceForMjjj (float const Mjjj)
{
  // Get the acceptance for a given mass
  //TF1* f = Accept->GetFunction("fit1");
  //f->GetParameter(0)
  //f->GetParameter(1)
  //f->GetParameter(2)
  //2nd degree polynomial numbers:
  //p0                        =   -0.0173967   +/-   0.000576407 
  //p1                        =  8.54121e-05   +/-   2.28888e-06 
  //p2                        = -2.44194e-08   +/-   1.63146e-09 
  //
  //3rd degree polynomial numbers:
  //p0                        =     -0.01027   +/-   0.00171542  
  //p1                        =  4.38331e-05   +/-   1.08208e-05 
  //p2                        =  4.43791e-08   +/-   1.96505e-08 
  //p3                        = -3.69426e-11   +/-   9.06989e-12


  return -0.01027 + 4.38331e-05*Mjjj + 4.43791e-08*Mjjj*Mjjj - 3.69426e-11*Mjjj*Mjjj*Mjjj;
  //return -0.0173967 + 8.54121e-05 * Mjjj + -2.44194e-08 * Mjjj * Mjjj;

}



std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  return std::make_pair<float, float>(Mjjj * 0.065, Mjjj * 0.075);
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

  if (in < 1) {
    return 0;
  }

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
  //for (int i = 0; i != 6; ++i) {
  //  printf("Par: %i %7E %7E %7E\n", i, Par[i], ParC[i], ParE[i]);
  //}


  // These are the parameters of the fit...
  // ("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + [3]*(((1.-x/7000.)^[4])/((x/7000.)^[5]))",  170, 800);
  // Systematics: start on the 5.  should be gaussian centered at 0 width one
  for (int i = 0; i != NPAR; ++i) {
    fToFit->SetParameter(i, Par[i]);
  }

  // LL variable
  static long double LogLikelihood;
  LogLikelihood = 0.0;
  static long double mu = 0.0;

  // Loop over all bins in histogram
  static int const iStart = hToFit->FindBin(XMIN);
  static int const iStop  = hToFit->FindBin(XMAX);
  for (int ibinX = iStart; ibinX <= iStop; ++ibinX) {

    mu = fToFit->Integral( hToFit->GetBinLowEdge(ibinX), hToFit->GetBinLowEdge(ibinX) + hToFit->GetBinWidth(ibinX)) / hToFit->GetBinWidth(ibinX);
    //printf("Bin X-X mu: %4i %9.1f-%9.1f %12.3E\n", ibinX, hToFit->GetBinLowEdge(ibinX), hToFit->GetBinLowEdge(ibinX) + hToFit->GetBinWidth(ibinX), (float) mu);

    if (mu > 0.1) {
      LogLikelihood += (hToFit->GetBinContent(ibinX) * TMath::Log(mu)
          - mu - LogFactorial((int) hToFit->GetBinContent(ibinX)  ) );
    }
  }

  for (int i = 0; i < NParameters; ++i) {
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

  // These are for Minuit
  double ArgList[10];
  int ErrorFlag;

  // Grab a new instance of Minuit with NParams
  TMinuit MyMinuit(NPAR);

  // Set the printlevel -2 minimum, higher for more junk
  MyMinuit.SetPrintLevel(-2);

  // Get the Expected numbers from the input function and other params and limits
  // NOTE: You should always stay AWAY from the limits... ie the fit should NOT
  // hit the limit
  std::pair<float, float> gWidth = GetGausWidthRange(SignalMass);
  float const p3     = fFunc->GetParameter(0);
  float const p3_min = 0;
  float const p3_max = 6000;
  float const p4     = fFunc->GetParameter(1);
  float const p4_min =    0;
  float const p4_max =  100;
  float const p5     = fFunc->GetParameter(2);
  float const p5_min = -2;
  float const p5_max =  0;

  // Errors...
  float const ep3 = ep3 * 0.03;
  //float const ep4 = fFunc->GetParError(4);
  //float const ep5 = fFunc->GetParError(5);


  // This is the function we're working with
  // ("fSigBG","[0]*TMath::Gaus(x, [1], [2], 1) + [3]*(((1.-x/7000.)^[4])/((x/7000.)^[5]))",  170, 800);
  if (IsBGOnly) {
    // This is for BG only fits
    ParC[0] = -999;
    ParC[1] = SignalMass;
    ParC[2] = (gWidth.first+gWidth.second)/2;
    ParC[3] = fFunc->GetParameter(0);
    ParC[4] = fFunc->GetParameter(1);
    ParC[5] = fFunc->GetParameter(2);

    ParE[0] = -999;
    ParE[1] = -999;
    ParE[2] = -999;
    ParE[3] = fFunc->GetParError(0);
    ParE[4] = fFunc->GetParError(1);
    ParE[5] = fFunc->GetParError(2);
    if (!DoSyst) {
      ParE[3] = -999;
      ParE[4] = -999;
      ParE[5] = -999;

      double Pars[NPAR] = {ParC[0], ParC[1], ParC[2], ParC[3], ParC[4], ParC[5]};
      double *Blank;
      double NLL = 0;

      int NPar = NPAR;
      NegativeLogLikelihood(NPar, Blank, NLL, Pars, 0);
      return NLL;
    }

    MyMinuit.DefineParameter(0, "GausNorm", 0, 0.001, 0, 1000);
    MyMinuit.DefineParameter(1, "GausMean",  SignalMass, 0.01, SignalMass-1, SignalMass+1);
    MyMinuit.DefineParameter(2, "GausWidth", (gWidth.first+gWidth.second)/2., 0.001, gWidth.first, gWidth.second);
    MyMinuit.DefineParameter(3, "BGp3", p3, 0.001, p3_min, p3_max);
    MyMinuit.DefineParameter(4, "BGp4", p4, 0.000001, p4_min, p4_max);
    MyMinuit.DefineParameter(5, "BGp5", p5, 0.000001, p5_min, p5_max);
    MyMinuit.FixParameter(0);
    MyMinuit.FixParameter(1);
    MyMinuit.FixParameter(2);
    if (!DoSyst) {
      MyMinuit.FixParameter(3);
      MyMinuit.FixParameter(4);
      MyMinuit.FixParameter(5);
    }
  } else {
    // This is for Signal+BG fits
    ParC[0] = -999;
    ParC[1] = SignalMass;
    ParC[2] = (gWidth.first+gWidth.second)/2;
    ParC[3] = fFunc->GetParameter(0);
    ParC[4] = fFunc->GetParameter(1);
    ParC[5] = fFunc->GetParameter(2);

    ParE[0] = -999;
    ParE[1] = -999;
    ParE[2] = -999;
    ParE[3] = fFunc->GetParError(0);
    ParE[4] = fFunc->GetParError(1);
    ParE[5] = fFunc->GetParError(2);
    if (!DoSyst) {
      ParE[3] = -999;
      ParE[4] = -999;
      ParE[5] = -999;
    }

    MyMinuit.DefineParameter(0, "GausNorm", 0, 0.001, 0, 1000);
    MyMinuit.DefineParameter(1, "GausMean",  SignalMass, 0.1, SignalMass-1, SignalMass+1);
    MyMinuit.DefineParameter(2, "GausWidth", (gWidth.first+gWidth.second)/2., 0.01, gWidth.first, gWidth.second);
    MyMinuit.DefineParameter(3, "BGp3", p3, 0.1, p3_min, p3_max);
    MyMinuit.DefineParameter(4, "BGp4", p4, 0.01, p4_min, p4_max);
    MyMinuit.DefineParameter(5, "BGp5", p5, 0.00001, p5_min, p5_max);
    MyMinuit.FixParameter(1);
    if (!DoSyst) {
      MyMinuit.FixParameter(3);
      MyMinuit.FixParameter(4);
      MyMinuit.FixParameter(5);
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
  MyMinuit.SetMaxIterations(10000);

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
    hToFit->SetAxisRange(XMIN, XMAX, "X");
    hToFit->Draw("hist");

    fToFit->Draw("same");
    c.SaveAs(BUFF);
  }

  // Call for minos error analysis
  // nah, screw that
  // MyMinuit.mnmnos();

  // Grab parameters from minuit
  double OPar[NPAR], OEPar[NPAR];
  for (int i = 0; i < NPAR; ++i) {
    MyMinuit.GetParameter(i, OPar[i], OEPar[i]);
    printf("Fitted Param: %i %9E  +/-  %9E\n", i, OPar[i], OEPar[i]);
  }


  // For the parameters which minimize the NLL, lets get the NLL value and return it
  double NLL = 0;
  double* Blank;
  int NPar = NPAR;
  NegativeLogLikelihood(NPar, Blank, NLL, OPar, 0);

  return NLL;
}




















int RunPValue (TString const InFileName, int const Section, bool const DoSyst)
{
  // Some basic parameters
  float const StepSize    =  10;
  float const BeginMass   =  250;
  float const EndMass     = 1500;
  int   const NPerSection =   50;


  // Set the randome seed based on section number
  gRandom->SetSeed(12332 * (Section + 2));


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
  TH1F* DataTH1F = (TH1F*) InFile.Get("Mjjj_70_20_160_6jet");
  TF1* fFunc = DataTH1F->GetFunction("g3");

  // Number of events from data
  //int const NFromData = (int) DataTH1F->Integral(23, 152);

  // Set the functino to fit.  Be VERY careful if oyu change this... you need to change a lot
  // of other things too...
  fToFit = new TF1("FitFunction", "[0]*TMath::Gaus(x, [1], [2], 1) + [3]*(((1.-x/7000.)^[4])/((x/7000.)^[5]))", XMIN, XMAX);

  if (Section == -1) {
    // Do some Data Fit

    // In order to do this correctly you need a clone...ya I know..
    TH1F* DataClone = (TH1F*) DataTH1F->Clone("DataClone");
    hToFit = DataClone;

    // Grab the background only LL
    float const BGLL = 123;//-1.0 * MinimizeNLL(Section, -1, -999, fFunc, true, DoSyst);
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
      float const NPE = fFunc->Integral(XMIN, XMAX) / BINWIDTH;
      printf("iep: %i NPE: %12.2f\n", ipe, NPE);

      // Set parameters and errors
      std::vector< std::pair<float, float> > ParAndE;
      for (int ipar = 0; ipar != NPAR; ++ipar) {
        ParAndE.push_back( std::make_pair<float, float>(fFunc->GetParameter(ipar), fFunc->GetParError(ipar)) );
      }

      // Get PE
      hToFit = GetPE(NPE, ParAndE, ipe, DoSyst);

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
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Section = atoi(argv[2]);

  // Run Systematics?
  bool const DoSyst = false;

  RunPValue(InFileName, Section, DoSyst);

  return 0;
}
