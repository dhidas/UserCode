////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Fri Nov 12 11:52:11 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TStyle.h"

#include "MultiJets/MultiJetFitter/interface/TAnaHist.h"


// Let's maks this self contained and put all function and class declaritions here

// Object for fit stuff
class FitObj;

// In the order they appear in.
std::pair<float, float> GetGausWidthRange (float const);
int   GetDiagForMjjj (float const);
float GetAcceptanceForMass (float const);
float Luminosity ();
TH1D* GetPE (FitObj const&);
TH1D* GetHistForMjjj (float const, TFile*);
void  MakeGraph (std::vector< std::pair<float, float> > const&, TString const, TString const, TString const, TString const);

template <class T> double KahanSummation(T, T);
long double LogFactorial (int const);

void  SetFitObjParams (FitObj&);

std::pair<float, float> BestFitSigBG (FitObj const&, bool const ReturnTotal = false);

float LimitAtMass (FitObj const&);
int   RunMultiJetLimits (int const, TString const);


enum kFitType { kFitLandau, kFitExp, kFitDijet };

// Some global variables just to make my life easier
kFitType gFitType;
TFile* gOutFile;







class FitObj : public TObject
{
  // This is just a simple object to use in the fit.
  // Just easier to pass things that way I think...

  public:
    FitObj () {};
    ~FitObj () {};

    TFile* File;
    TH1D* Hist;
    float gmean;
    float gsigma;
    std::pair<float, float> gsigmaRange;
    bool  IsData;
    int   Section;
    int   ipe;
    int   nbins;
    float xmin;
    float xmax;
    bool  DoSyst;
    bool  DoAccSmear;
    TString LimitsFileName;

    TF1* Function;
};



std::vector< std::pair<float, float> > ReadLimitsFile (TString const FileName)
{
  // Vector we'll return
  std::vector< std::pair<float, float> > Vec;

  // string
  std::string word;

  // Open the input file
  std::ifstream In(FileName.Data());
  if (!In.is_open()) {
    std::cerr << "ERROR: cannot open limits file for reading" << std::endl;
    exit(1);
  }

  // Skip first line which is title
  std::getline(In, word);

  // Get the first line
  std::string line;
  std::getline(In, line);
  std::istringstream linestream;
  linestream.str(line);

  for (int i; linestream >> word; ++i) {
    if (linestream.eof()) {
      break;
    }
    Vec.push_back( std::make_pair<float, float>(atof(word.c_str()), -1) );
  }

  // Skip two lines
  std::getline(In, word);
  std::getline(In, word);

  // Sure, read the values..
  for (size_t i = 0; i != Vec.size(); ++i) {
    In >> Vec[i].second;
  }

  // Print just to check the are correct
  if (true) {
    for (size_t i = 0; i != Vec.size(); ++i) {
      printf("Read Meadian Limits  M = %8.1f    Limits = %7.3f\n", Vec[i].first, Vec[i].second);
    }
  }

  // Close the input file
  In.close();

  return Vec;
}



float GetPESignalXSForMass (FitObj const& MyFitObj)
{
  static std::vector< std::pair<float, float> > Vec = ReadLimitsFile(MyFitObj.LimitsFileName);

  for (size_t i = 0; i != Vec.size(); ++i) {
    if (Vec[i].first == MyFitObj.gmean) {
      return Vec[i].second;
    }
  }

  std::cerr << "ERROR: did not find mass match for M=" << MyFitObj.gmean << std::endl;

  return 9999;
}




std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  // This will eventually be a parametrization

  return std::make_pair<float, float>(15, 15);
  if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  if (Mjjj < 350) return std::make_pair<float, float>(10 + 10 * (Mjjj - 250.) / 100., 15 + 10 * (Mjjj - 250.) / 100.);
  //if (Mjjj < 350) return std::make_pair<float, float>(15, 20);
  return std::make_pair<float, float>(20, 25);
}



int GetDiagForMjjj (float const Mjjj)
{
  // Let's say you want to get the diag cut for a given mass..
  // I was using this to choose which hist to use..

  if (Mjjj < 150)      return 140;
  else if (Mjjj < 200) return 160;
  else if (Mjjj < 250) return 180;
  else if (Mjjj < 300) return 200;
  else if (Mjjj < 350) return 220;
  else if (Mjjj < 400) return 240;
  else                 return 260;
}



float GetAcceptanceForMass (float const Mjjj)
{
  // Get the acceptance for a given mass

  // Acceptance numbers from dan
  // Pol0:  1.93282628706596682e-03
  // Pol1: -3.61649088607910918e-05
  // Pol2:  1.84254776396677824e-07
  //TF1 AcceptanceFunc("AcceptanceFunc", "pol2(0)", BeginMass, EndMass);
  //AcceptanceFunc.SetParameter(0,  1.93282628706596682e-03);
  //AcceptanceFunc.SetParameter(1, -3.61649088607910918e-05);
  //AcceptanceFunc.SetParameter(2,  1.84254776396677824e-07);
  //return AcceptanceFunc.Eval(Mjjj);

  // This is equivalent to the above, but a hell of a lot faster
  if (Mjjj <= 250) {
    return -5.56179969055969892e-03 - 4.01623842089755741e-06 * Mjjj + 2.49580149009780901e-07 * Mjjj * Mjjj;
  } else {
    return 1.93282628706596682e-03 - 3.61649088607910918e-05 * Mjjj + 1.84254776396677824e-07 * Mjjj * Mjjj;
  }

}



float Luminosity ()
{
  return 35.0;
}



TH1D* GetPE (FitObj const& Obj)
{
  // This function shall return a PE given the parameters in the FitObj
  // YOU are the owner of this PE so please delete it.

  // Landau function!
  TF1* Func = (TF1*) Obj.Function->Clone();

  // If we're doing systematics let's add some randomness.  These numbers taken from
  // the CDF code
  if (Obj.DoSyst) {
    switch (gFitType) {
      case kFitLandau:
        Func->FixParameter(0, Func->GetParameter(0) * (1.0 + gRandom->Gaus(0, 0.10)));
        Func->FixParameter(1, Func->GetParameter(1) * (1.0 + gRandom->Gaus(0, 0.01)));
        Func->FixParameter(2, Func->GetParameter(2) * (1.0 + gRandom->Gaus(0, 0.10)));
        break;
      case kFitExp:
        Func->FixParameter(0, Func->GetParameter(0) * (1.0 + gRandom->Gaus(0, 0.10)));
        Func->FixParameter(1, Func->GetParameter(1) * (1.0 + gRandom->Gaus(0, 0.01)));
        break;
      case kFitDijet:
        Func->FixParameter(0, Func->GetParameter(0) * (1.0 + gRandom->Gaus(0, 0.10)));
        Func->FixParameter(1, Func->GetParameter(1) * (1.0 + gRandom->Gaus(0, 0.01)));
        Func->FixParameter(2, Func->GetParameter(2) * (1.0 + gRandom->Gaus(0, 0.10)));
        break;
      default:
        std::cerr << "Uhm, dude, you should not be here and now I am going to barf on you" << std::endl;
        exit(1);
        break;
    }
  } else {
    // Do nothing!
  }

  // Create a histogram
  char BUFF[40];
  sprintf(BUFF, "PE_%i_%i", (int) Obj.gmean, Obj.ipe);
  TH1D* hPE = new TH1D(BUFF, BUFF, Obj.nbins, Obj.xmin, Obj.xmax);

  // Fill the histogram
  float xval;
  for (int ibin = 1; ibin <= Obj.nbins; ++ibin) {
    xval = Obj.xmin + (Obj.xmax - Obj.xmin) / ((float) Obj.nbins) * (ibin - 0.5);
    hPE->SetBinContent(ibin, gRandom->Poisson( Func->Eval(xval) ) );
  }

  // If we're adding signal to this..
  if (Obj.LimitsFileName != "") {
    float const Mass = Obj.gmean;
    int const NSignal = (int) (GetPESignalXSForMass(Obj) / 0.683 * Luminosity() * GetAcceptanceForMass(Mass));
    std::cout << "MM " << Mass << "  " << NSignal << std::endl;
    std::pair<float, float> const WidthRange = Obj.gsigmaRange;
    float const Width = gRandom->Uniform(WidthRange.first, WidthRange.second);
    TF1 Gaus("Gaus", "[0] * TMath::Gaus(x, [1], [2], 1)", Mass - Width, Mass + Width);
    Gaus.FixParameter(0, NSignal);
    Gaus.FixParameter(1, Mass);
    Gaus.FixParameter(2, Width);

    hPE->FillRandom("Gaus", NSignal);
  }

  TCanvas Can(hPE->GetName(), hPE->GetName());
  Can.cd();
  hPE->Draw();
  Func->Draw("same");

  if (gOutFile) {
    gOutFile->cd();
    Can.Write();
  }

  delete Func;

  return hPE;
}



TF1* GetFitForMjjj (float const Mjjj, TFile* File, TString const FitName)
{
  // Get the fit from the hist for this
  TH1D* Hist = GetHistForMjjj(Mjjj, File);
  TF1* Func = Hist->GetFunction(FitName);


  if (false) {
    TF1* Func2 = new TF1("land", "[0] * TMath::Landau(x, [1], [2], 1)", Func->GetXmin(), Func->GetXmax());
    Func2->FixParameter(0, Func->GetParameter(0) * Func->GetParameter(2));
    Func2->FixParameter(1, Func->GetParameter(1));
    Func2->FixParameter(2, Func->GetParameter(2));

    //delete Func;
    return Func2;
  }

  return Func;
}



TH1D* GetHistForMjjj (float const Mjjj, TFile* File)
{
  // This function will return a histogram
  // Before, I had put all PEs in a file and ran over them.. now generate
  // on the go.. so this may or may not be used for more than data..

  char BUFF[200];
  if (Mjjj <= 250)        sprintf(BUFF,  "Mjjj_45_20_130_6jet");
  else                    sprintf(BUFF,  "Mjjj_45_20_170_6jet");
  //sprintf(BUFF,  "Mjjj_45_20_170_6jet");


  std::cout << "Getting Hist: " << BUFF << "  for mass " << Mjjj << std::endl;
  TH1D* Hist = (TH1D*) File->Get(BUFF);
  if (!Hist) {
    std::cerr << "ERROR: cannot get hist " << BUFF << std::endl;
    exit(1);
  }

  return Hist;
}


void MakeGraph (std::vector< std::pair<float, float> > const& P, TString const Title, TString const XTitle, TString const YTitle, TString const OutFileName)
{
  // This function will make a TGraph and save it as a file given any vector of float pairs.
  // Give it some titles, and an output name..

  // arrays for the x and y points
  float grX[P.size()];
  float grY[P.size()];

  // Fill the arrays
  for (size_t ifit = 0; ifit != P.size(); ++ifit) {
    printf("Summary * GausMean = %12.4f  ==>  %12.4f\n", P[ifit].first, P[ifit].second);
    grX[ifit] = P[ifit].first;
    grY[ifit] = P[ifit].second;
  }


  // Why not make a grap of that..?
  TGraph Graph(P.size(), grX, grY);
  TCanvas Can;
  Graph.SetTitle(Title);
  Graph.GetXaxis()->SetTitle(XTitle);
  Graph.GetYaxis()->SetTitle(YTitle);
  Graph.Draw("AL*");
  Can.SaveAs(OutFileName);

  return;
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



void SetFitObjParams (FitObj& MyFitObj)
{
  // Set the fig object parameters however you want

  // Some default parameters.  I could change how this is done someday because it shoudl be automatic
  // ie done from the data hist.. I'll get to that later I guess

  TH1D* DataHist = GetHistForMjjj(MyFitObj.gmean, MyFitObj.File);
  MyFitObj.nbins      = DataHist->GetNbinsX();
  MyFitObj.xmin       = DataHist->GetXaxis()->GetXmin();
  MyFitObj.xmax       = DataHist->GetXaxis()->GetXmax();
  MyFitObj.DoSyst     = true;
  MyFitObj.DoAccSmear = true;

  MyFitObj.Function   = GetFitForMjjj(MyFitObj.gmean, MyFitObj.File, "total");


  //TCanvas Can;
  //MyFitObj.Hist->Draw("hist");
  //Land.Draw("same");
  //Can.SaveAs("Fit.eps");

  //printf("LandauFit: %8.2f %8.2f %8.2f\n", MyFitObj.nland, MyFitObj.lmpv, MyFitObj.lsigma);


  return;
}




std::pair<float, float> BestFitSigBG (FitObj const& Obj, bool const ReturnTotal)
{

  TString BGFuncStr;
  switch (gFitType) {
    case kFitLandau:
      BGFuncStr = "[0] * TMath::Landau(x, [1], [2], 1)";
      break;
    case kFitExp:
      BGFuncStr = "TMath::Exp([0] + [1] * x) + 0*[2]";
      break;
    case kFitDijet:
      BGFuncStr = "[0] * ( ((1.0 - x/7000.0)^[1]) / (x^[2]) )";
      break;
    default:
      std::cerr << "Uhm, dude, you should not be here and now I am going to barf on you" << std::endl;
      exit(1);
      break;
  }

  // Get out function
  TF1 FuncGaus("FuncGaus", BGFuncStr + " + [3] * TMath::Gaus(x, [4], [5], 1)", Obj.Function->GetXmin(), Obj.Function->GetXmax());

  // This error taken from CDF code
  float err=0.000000001;

  // Set the parameters depending on which fit it is
  switch (gFitType) {
    case kFitLandau:
      FuncGaus.SetParameter(0, Obj.Function->GetParameter(0) * Obj.Function->GetParameter(2));
      FuncGaus.SetParLimits(0, FuncGaus.GetParameter(0) - err, FuncGaus.GetParameter(0) + err);
      FuncGaus.SetParameter(1, Obj.Function->GetParameter(1));
      FuncGaus.SetParLimits(1, Obj.Function->GetParameter(1) - err, Obj.Function->GetParameter(1) + err);
      FuncGaus.SetParameter(2, Obj.Function->GetParameter(2));
      FuncGaus.SetParLimits(2, Obj.Function->GetParameter(2) - err, Obj.Function->GetParameter(2) + err);
      break;
    case kFitExp:
      FuncGaus.SetParameter(0, Obj.Function->GetParameter(0));
      FuncGaus.SetParLimits(0, FuncGaus.GetParameter(0) - err, FuncGaus.GetParameter(0) + err);
      FuncGaus.SetParameter(1, Obj.Function->GetParameter(1));
      FuncGaus.SetParLimits(1, Obj.Function->GetParameter(1) - err, Obj.Function->GetParameter(1) + err);
      FuncGaus.FixParameter(2, 0);
      break;
    case kFitDijet:
      FuncGaus.SetParameter(0, Obj.Function->GetParameter(0));
      FuncGaus.SetParLimits(0, FuncGaus.GetParameter(0) - err, FuncGaus.GetParameter(0) + err);
      FuncGaus.SetParameter(1, Obj.Function->GetParameter(1));
      FuncGaus.SetParLimits(1, Obj.Function->GetParameter(1) - err, Obj.Function->GetParameter(1) + err);
      FuncGaus.SetParameter(2, Obj.Function->GetParameter(2));
      FuncGaus.SetParLimits(2, FuncGaus.GetParameter(2) - err, FuncGaus.GetParameter(2) + err);
      break;
    default:
      std::cerr << "Uhm, dude, you should not be here and now I am going to barf on you" << std::endl;
      exit(1);
      break;
  }

  // Set the gaussian parameters
  FuncGaus.SetParameter(3, 0);
  FuncGaus.SetParLimits(3, -1000, 1000);
  FuncGaus.SetParameter(4, Obj.gmean);
  FuncGaus.SetParLimits(4, Obj.gmean - 0.1, Obj.gmean + 0.1);
  FuncGaus.SetParameter(5, (Obj.gsigmaRange.first + Obj.gsigmaRange.second)/2);
  FuncGaus.SetParLimits(5, Obj.gsigmaRange.first, Obj.gsigmaRange.second);


  // Fit the function to the histogram given
  Obj.Hist->Fit("FuncGaus", "QL", "", FuncGaus.GetXmin(), FuncGaus.GetXmax());


  // Right here if you want the total integral no more is needed since the gaussian is normalized
  if (ReturnTotal) {
    // Note that the landau is an approximation since the function is cut off at some x value..
    // this returns the integral to infinity for both...
    return std::make_pair<float, float>(FuncGaus.GetParameter(3) / Obj.Hist->GetBinWidth(0),
        FuncGaus.GetParameter(0) / Obj.Hist->GetBinWidth(0));
  }


  // Define landau and gaus given the parameters fo the fit
  TF1 BGFunc("BGFunc", BGFuncStr, Obj.Function->GetXmin(), Obj.Function->GetXmax());
  TF1 Gaus("Gaus", "[0]*TMath::Gaus(x, [1], [2], 1)", Obj.Function->GetXmin(), Obj.Function->GetXmax());
  BGFunc.FixParameter(0, FuncGaus.GetParameter(0));
  BGFunc.FixParameter(1, FuncGaus.GetParameter(1));
  BGFunc.FixParameter(2, FuncGaus.GetParameter(2));
  Gaus.FixParameter  (0, FuncGaus.GetParameter(3));
  Gaus.FixParameter  (1, FuncGaus.GetParameter(4));
  Gaus.FixParameter  (2, FuncGaus.GetParameter(5));

  // Get the signal in +/- 1 sigma
  float const NGausTotal = FuncGaus.GetParameter(3);
  float const Sig = 0.683 * NGausTotal / Obj.Hist->GetBinWidth(0);

  // Need to know where to integrate the landau from and to
  float const begin = Obj.gmean - Gaus.GetParameter(2);
  float const end   = Obj.gmean + Gaus.GetParameter(2);

  // Get the background
  float const BG = BGFunc.Integral(begin, end) / Obj.Hist->GetBinWidth(0);


  // If you want to see the signal and BG totals uncomment the following line
  printf("SigBG Total: %12.3f %12.3f %12.3f\n", Sig, BG, BGFunc.Integral(Obj.xmin, Obj.xmax) / Obj.Hist->GetBinWidth(0));


  bool const DoAllPlots = false;
  if (Obj.IsData || DoAllPlots) {
    // This just saves the plot for data, and optionally for all...if you reall want

    char FuncGausHistName[100];
    if (Obj.Section < 0) {
      sprintf(FuncGausHistName, "FuncGaus_%i_Data", (int) Obj.gmean);
    } else {
      sprintf(FuncGausHistName, "FuncGaus_%i_%i", (int) Obj.gmean, Obj.ipe);
    }

    TCanvas Can(FuncGausHistName, FuncGausHistName);
    Can.cd();
    Obj.Hist->Draw();
    FuncGaus.Draw("same");
    BGFunc.SetLineStyle(2);
    BGFunc.Draw("same");
    Gaus.SetLineColor(2);
    Gaus.Draw("same");

    if (Obj.Section < 0) {
      sprintf(FuncGausHistName, "FuncGaus_%i_Data.eps", (int) Obj.gmean);
    } else {
      sprintf(FuncGausHistName, "FuncGaus_%i_%i.eps", (int) Obj.gmean, Obj.ipe);
    }
    if (Obj.IsData) {
      Can.SaveAs(FuncGausHistName);
    }
    if (gOutFile) {
      gOutFile->cd();
      Can.Write();
    }

    std::cout << "MYLimit norms: " << Obj.Hist->Integral() << "  "
      << FuncGaus.GetParameter(0) << "  " 
      << FuncGaus.GetParameter(1) << "  " 
      << FuncGaus.GetParameter(2) << "  " 
      //<< Gaus.Integral(0, 3000) << "  "
      //<< LandGaus.Integral(0, 3000)
      << std::endl;
  }

  // return the signal and background pair
  return std::make_pair<float, float>(Sig, BG);

}



float LimitAtMass (FitObj const& Obj)
{
  // Get the 95% limit given the best fit signal and BG.  This uses
  // a poisson likelihood, which is integrated to 95%..

  // You have to define where to integrate to and from.  The important part
  // is the end or cutooff.  It must be practical to trust this.
  float const begin = 0;
  float const end = 60;

  // Get the best fir signal and background
  std::pair<float, float> SigBG = BestFitSigBG(Obj);

  // Define a poisson likelhood given the best fit parameters
  TF1 Poisson("land", "[0]*TMath::Poisson([1], [2]+x)", begin, end);
  Poisson.FixParameter(0, 1);
  Poisson.FixParameter(1, SigBG.first + SigBG.second);
  Poisson.FixParameter(2, SigBG.second);

  // This is to normalize the poisson
  float const Total = Poisson.Integral(begin, end);
  Poisson.FixParameter(0, 1/Total);

  // Integrate to 95% =)
  float mean = 0;
  for ( ; Poisson.Integral(begin, mean) < 0.95; mean += 0.01) {}



  if (Obj.IsData) {
    // Let's draw these for data

    float const Max = Poisson.GetMaximum();
    float const Mass = Obj.gmean;
    char HistName[200];
    sprintf(HistName, "CL95_%i.eps", (int) Mass);
    TCanvas Can;
    Poisson.Draw();
    TLine Line(mean, 0, mean, Max);
    Line.SetLineColor(2);
    Line.SetLineWidth(2);
    Line.Draw("same");
    Can.SaveAs(HistName);
  }

  // Teturn the 95% limit
  return mean;
}








int RunMultiJetLimits (int const Section, TString const InFileName, TString const LimitsFileName)
{
  // Set the default style to plain to not waste ink
  gROOT->SetStyle("Plain");

  // Open the inpu file
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }


  // Open the output file
  char OutFileName[200];
  if (Section < 0) {
    sprintf(OutFileName, "LimitPlots_Data.root");
  } else {
    sprintf(OutFileName, "LimitPlots_%i.root", Section);
  }

  // If you want to save all the plots to a root file
  // make tis true...
  gOutFile = 0x0;
  if (false) {
    gOutFile = new TFile(OutFileName, "recreate");
    if (!gOutFile->IsOpen()) {
      std::cerr << "ERROR: cannot open Output file: " << OutFileName << std::endl;
      exit(1);
    }
  }

  // Set the random seed.
  //gRandom->SetSeed(5 * Section);
  gRandom->SetSeed( (Section + 2) * (int) fmod(time(NULL),100000));
  std::cout << "Setting seed to: " << gRandom->GetSeed() << std::endl;



  // Setup the output text file
  // Output file name
  if (Section < 0) {
    sprintf(OutFileName, "DataLimits.dat");
  } else {
    sprintf(OutFileName, "PELimits_%i.dat", Section);
  }

  // Open the output text file for this section
  FILE*  OutFile = fopen(OutFileName, "w");
  if (OutFile == NULL) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
    exit(1);
  }

  // These numbers are what range you want to calculate limits for.
  float const BeginMass = 200;
  float const EndMass   = 500;
  float const StepSize  =  10;


  float const AcceptErr  = 0.30;

  // Print the masses to the text file
  for (float ThisMass = BeginMass; ThisMass <= EndMass; ThisMass += StepSize) {
    fprintf(OutFile, "%10.3f ", ThisMass);
  }
  fprintf(OutFile, "\n");


  // Setup FitObj
  FitObj MyFitObj;
  MyFitObj.File = &InFile;
  MyFitObj.LimitsFileName = LimitsFileName;


  // Which section is this... data or PE?
  if (Section < 0) {
    // This is data

    // Something to keep track of the limits
    std::vector< std::pair<float, float> > GausMeanNGaus;

    // Run through the masses
    for (float ThisMass = BeginMass; ThisMass <= EndMass; ThisMass += StepSize) {

      // Grab the data hist
      TH1D* DataTH1 = GetHistForMjjj(ThisMass, &InFile);
      if (!DataTH1) {
        std::cerr << "ERROR: cannot get data histogram: " << std::endl;
        exit(1);
      }

      // Setup the fit object
      MyFitObj.Hist = (TH1D*) DataTH1->Clone(0);
      MyFitObj.gmean = ThisMass;
      SetFitObjParams(MyFitObj);
      MyFitObj.gsigmaRange = GetGausWidthRange(ThisMass);
      MyFitObj.IsData = true;
      MyFitObj.Section = Section;


      // Calcuate the limit for this mass
      float const Limit = LimitAtMass(MyFitObj);

      // Technically I guess a cross section//
      float const Acceptance = GetAcceptanceForMass(ThisMass);
      float const XSecLimit = Limit / (Luminosity() * Acceptance);

      // Print that out
      printf("MYLimit %9i %12.4f\n", (int) ThisMass, XSecLimit);
      fprintf(OutFile, "%10E ", XSecLimit);

      // Save it to our little friendly vector
      GausMeanNGaus.push_back(std::make_pair<float, float>(ThisMass, XSecLimit));

      // delete our new hist
      delete MyFitObj.Hist;
    }

    // May as well plot the limit as a functino of mass..
    MakeGraph (GausMeanNGaus, "95% C.L. Number of signal events", "Gaus mean [GeV]", " Number of signal events", "Limit95NEvents_Data.eps");

  } else {
    // This is for PEs

    // This just defines how many to run per section
    int const NPerSection = 100;
    int const Start = Section * NPerSection;
    int const End   = Start   + NPerSection;

    // Set the basic fit obj parameters
    MyFitObj.IsData = false;
    MyFitObj.Section = Section;

    // Let's run some PEs
    for (int ipe = Start; ipe != End; ++ipe) {

      // Keep track of the best fit by mass for each PE just in case
      std::vector< std::pair<float, float> > GausMeanBestFit;

      // Set the PE number and get a PE
      MyFitObj.ipe = ipe;

      // Chec to see that we have a hist
      if (!MyFitObj.Hist) {
        std::cerr << "ERROR: cannot get PE: " << ipe << std::endl;
        exit(1);
      }



      // Do the masses
      for (float ThisMass = BeginMass; ThisMass <= EndMass; ThisMass += StepSize) {

        // Set the mass and width range
        MyFitObj.gmean = ThisMass;
        SetFitObjParams(MyFitObj);
        MyFitObj.gsigmaRange = GetGausWidthRange(ThisMass);

        // Get PE for this ipe and mass (this need to be here for when we add signal...)
        MyFitObj.Hist = GetPE(MyFitObj);


        // Grab the best fit given this PE and mass
        std::pair<float, float> SigBG = BestFitSigBG(MyFitObj);

        // Technically a cross section is okay too I guess
        // If you're doing acceptance smear do it here.
        float const Acceptance = GetAcceptanceForMass(ThisMass);
        float const ThisAcceptance = MyFitObj.DoAccSmear ? Acceptance * (1.0 + gRandom->Gaus(0, AcceptErr)) : Acceptance;
        float const XSec = SigBG.first / (Luminosity() * ThisAcceptance);

        // Print that out and save it to file
        printf("MYLimit XSec %9i %9i %12.4f\n", ipe, (int) ThisMass, XSec);
        fprintf(OutFile, "%10E ", XSec);

        // Why not put that in a vector
        GausMeanBestFit.push_back( std::make_pair<float, float>(ThisMass, XSec) );

        // I told you we should delete this
        delete MyFitObj.Hist;
      }
      fprintf(OutFile, "\n");
      fflush(OutFile);

      // Every so often output a graph just because
      if (ipe % 1000 == 0) {
        char FitHistName[100];
        sprintf(FitHistName, "Limit95NEvents_%i.eps", ipe);
        MakeGraph (GausMeanBestFit, "Best Fit Number of signal events", "Gaus mean [GeV]", "Number of signal events", FitHistName);
      }

    }
  }

  // Nicely close the text output file
  fclose(OutFile);

  // Close the out root file if it exists
  if (gOutFile) {
    gOutFile->Close();
    delete gOutFile;
  }

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 4 && argc != 5) {
    std::cerr << "Usage: " << argv[0] << " [FitType] [Section] [InFileName] [Optional Limits file]" << std::endl;
    return 1;
  }

  TString const FitTypeStr = argv[1];
  int const Section = atoi(argv[2]);
  TString const InFileName = argv[3];
  TString const LimitsFileName = argc == 5 ? argv[4] : "";

  if (Section < -1) {
    std::cerr << "Well, I really intended data to be -1 and PEs to be >= 0.  Be careful" << std::endl;
    return 1;
  }

  // Set the fit Ty[e
  if      (FitTypeStr == "Landau") gFitType = kFitLandau;
  else if (FitTypeStr == "Exp")    gFitType = kFitExp;
  else if (FitTypeStr == "Dijet")  gFitType = kFitDijet;
  else {
    std::cerr << "ERROR: I don't know what fit type you want" << std::endl;
    exit(1);
  }

  RunMultiJetLimits(Section, InFileName, LimitsFileName);

  return 0;
}
