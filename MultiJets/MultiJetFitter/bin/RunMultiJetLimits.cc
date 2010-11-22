////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Fri Nov 12 11:52:11 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TLine.h"





class FitObj : public TObject
{
  public:
    FitObj () {};
    ~FitObj () {};

    TH1D* Hist;

    float lmpv;
    float lsigma;
    float nland;
    float gmean;
    float gsigma;
    std::pair<float, float> gsigmaRange;
    bool IsData;
    int Section;
    int ipe;
};




std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj

  return std::make_pair<float, float>(8, 25);
  if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  if (Mjjj < 350) return std::make_pair<float, float>(15, 20);
  return std::make_pair<float, float>(20, 25);
}



int GetDiagForMjjj (float const Mjjj)
{
  if (Mjjj < 150)      return 140;
  else if (Mjjj < 200) return 160;
  else if (Mjjj < 250) return 180;
  else if (Mjjj < 300) return 200;
  else if (Mjjj < 350) return 220;
  else if (Mjjj < 400) return 240;
  else                 return 260;
}


TH1D* GetHistForMjjj (float const Mjjj, TFile* File, int const iPE)
{
  char BUFF[200];
  //if (iPE < 0)          sprintf(BUFF,  "Mjjj_20_20_%i", GetDiagForMjjj(Mjjj));
  //else                  sprintf(BUFF, "PE_20_20_%i_%i", GetDiagForMjjj(Mjjj), iPE);
  if (iPE < 0)          sprintf(BUFF,  "Mjjj_20_20_200");
  else                  sprintf(BUFF, "PE_20_20_200_%i", iPE);

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


template <class T>
double KahanSummation(T begin, T end) {
  double result = 0.f;

  double c = 0.f;
  double y, t;
  for(;begin != end; ++begin) {
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
    //LogSum += TMath::Log(i);
    Logs.push_back(TMath::Log(i));
  }

  LogSum = KahanSummation(Logs.begin(), Logs.end());

  return LogSum;
}



void SetFitObjParams (FitObj& MyFitObj)
{
  // Set the fig object parameters however you want

  if (false) {
    TF1 Land("land", "[0]*TMath::Landau(x, [1], [2], 1)", MyFitObj.Hist->GetXaxis()->GetXmin(), MyFitObj.Hist->GetXaxis()->GetXmax());
    Land.SetParameter(0, MyFitObj.Hist->GetEntries());
    Land.SetParameter(1, MyFitObj.Hist->GetMean());
    Land.SetParameter(2, MyFitObj.Hist->GetRMS());

    MyFitObj.Hist->Fit("land");

    MyFitObj.nland  = Land.GetParameter(0);
    MyFitObj.lmpv   = Land.GetParameter(1);
    MyFitObj.lsigma = Land.GetParameter(2);
  } else {
    MyFitObj.nland  = 209.245;
    MyFitObj.lmpv   = 176.320;
    MyFitObj.lsigma = 37.1101;
  }

  //TCanvas Can;
  //MyFitObj.Hist->Draw("hist");
  //Land.Draw("same");
  //Can.SaveAs("Fit.eps");

  //printf("LandauFit: %8.2f %8.2f %8.2f\n", MyFitObj.nland, MyFitObj.lmpv, MyFitObj.lsigma);


  return;
}




std::pair<float, float> BestFitSigBG (FitObj const& Obj)
{
  TF1 LandGaus("LandGaus", "[0]*TMath::Landau(x, [1], [2], 1) + [3]*TMath::Gaus(x, [4], [5], 1)", Obj.Hist->GetXaxis()->GetXmin(), Obj.Hist->GetXaxis()->GetXmax());
  //LandGaus.FixParameter(0, Obj.Hist->Integral() * Obj.Hist->GetBinWidth(0));
  LandGaus.SetParameter(0, Obj.nland);
  LandGaus.SetParLimits(0, Obj.nland * 0.97, Obj.nland * 1.03);
  LandGaus.FixParameter(1, Obj.lmpv);
  LandGaus.SetParameter(2, Obj.lsigma);
  LandGaus.SetParLimits(2, Obj.lsigma*0.97, Obj.lsigma*1.03);
  LandGaus.SetParameter(3, 200);
  LandGaus.SetParLimits(3, 0, 1000);
  LandGaus.FixParameter(4, Obj.gmean);
  LandGaus.SetParameter(5, (Obj.gsigmaRange.first + Obj.gsigmaRange.second)/2);
  LandGaus.SetParLimits(5, Obj.gsigmaRange.first, Obj.gsigmaRange.second);

  //Obj->Hist->Fit("land");
  TH1D* FitHist = (TH1D*) Obj.Hist->Clone();
  FitHist->Fit("LandGaus");

  TF1 Land("Land", "[0]*TMath::Landau(x, [1], [2], 1)", Obj.Hist->GetXaxis()->GetXmin(), Obj.Hist->GetXaxis()->GetXmax());
  TF1 Gaus("Gaus", "[3]*TMath::Gaus(x, [4], [5], 1)",   Obj.Hist->GetXaxis()->GetXmin(), Obj.Hist->GetXaxis()->GetXmax());
  Land.FixParameter(0, LandGaus.GetParameter(0));
  Land.FixParameter(1, LandGaus.GetParameter(1));
  Land.FixParameter(2, LandGaus.GetParameter(2));
  Gaus.FixParameter(3, LandGaus.GetParameter(3));
  Gaus.FixParameter(4, LandGaus.GetParameter(4));
  Gaus.FixParameter(5, LandGaus.GetParameter(5));

  float const begin = Obj.gmean - Gaus.GetParameter(5);
  float const end   = Obj.gmean + Gaus.GetParameter(5);

  float Sig = Gaus.Integral(begin, end) / Obj.Hist->GetBinWidth(0);
  float BG  = Land.Integral(begin, end) / Obj.Hist->GetBinWidth(0);

  if (Obj.IsData) {
    TCanvas Can;
    Land.Draw();
    Gaus.Draw("same");
    Obj.Hist->Draw("same");
    char LandGausHistName[100];
    if (Obj.Section < 0) {
      sprintf(LandGausHistName, "LandGaus_%i_Data.eps", (int) Obj.gmean);
    } else {
      sprintf(LandGausHistName, "LandGaus_%i_%i.eps", (int) Obj.gmean, Obj.ipe);
    }
    Can.SaveAs(LandGausHistName);
  }

  printf("SigBG LandTotal: %12.3f %12.3f %12.3f\n", Sig, BG, Land.Integral(0, 3000) / Obj.Hist->GetBinWidth(0));

  return std::make_pair<float, float>(Sig, BG);

}



float LimitAtMass (FitObj const& Obj)
{
  float const begin = 0;
  float const end = 100;

  std::pair<float, float> SigBG = BestFitSigBG(Obj);
  TF1 Poisson("land", "[0]*TMath::Poisson([1], [2]+x)", begin, end);
  Poisson.FixParameter(0, 1);
  Poisson.FixParameter(1, SigBG.first + SigBG.second);
  Poisson.FixParameter(2, SigBG.first);

  float const Total = Poisson.Integral(begin, end);
  Poisson.FixParameter(0, 1/Total);

  float mean = 0;
  for ( ; Poisson.Integral(begin, mean) < 0.95; mean += 0.1) {}


  float const Max = Poisson.GetMaximum();

  float const Mass = Obj.gmean;

  if (Obj.IsData) {
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

  return mean;
}








int RunMultiJetLimits (int const Section, TString const InFileName)
{
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }


  // Setup the output text file
  // Output file name
  char OutFileName[200];
  if (Section < 0) {
    sprintf(OutFileName, "DataLimits.dat");
  } else {
    sprintf(OutFileName, "Limits_%i.dat", Section);
  }

  FILE*  OutFile = fopen(OutFileName, "w");
  if (OutFile == NULL) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
    exit(1);
  }


  float const BeginMass = 150;
  float const EndMass   = 450;
  float const StepSize  =   10;

  for (float ThisMass = BeginMass; ThisMass <= EndMass; ThisMass += StepSize) {
    fprintf(OutFile, "%10.3f ", ThisMass);
  }
  fprintf(OutFile, "\n");




  // Which section is this... data or PE?
  float Limit;
  if (Section < 0) {
    std::vector< std::pair<float, float> > GausMeanNGaus;

    for (float ThisMass = BeginMass; ThisMass < EndMass; ThisMass += StepSize) {
      TH1D* DataTH1 = GetHistForMjjj(ThisMass, &InFile, -1);
      if (!DataTH1) {
        std::cerr << "ERROR: cannot get data histogram: " << std::endl;
        exit(1);
      }

      FitObj MyFitObj;
      MyFitObj.Hist = DataTH1;
      SetFitObjParams(MyFitObj);
      MyFitObj.gmean = ThisMass;
      MyFitObj.gsigmaRange = GetGausWidthRange(ThisMass);
      MyFitObj.IsData = true;
      MyFitObj.Section = Section;

      Limit = LimitAtMass(MyFitObj);
      printf("DataLimit %9i %12.4f\n", (int) ThisMass, Limit);
      fprintf(OutFile, "%10E ", Limit);
      GausMeanNGaus.push_back(std::make_pair<float, float>(ThisMass, Limit));
    }
    MakeGraph (GausMeanNGaus, "95% C.L. Number of signal events", "Gaus mean [GeV]", " Number of signal events", "LimitNEvents.eps");

  } else {

    int const NPerSection = 250; // For 70 sections
    //int const NPerSection = 500; // For 200
    int const Start = Section * NPerSection;
    int const End   = Start   + NPerSection;

    for (int ipe = Start; ipe != End; ++ipe) {

      std::vector< std::pair<float, float> > GausMeanNGaus;

      for (float ThisMass = BeginMass; ThisMass < EndMass; ThisMass += StepSize) {
        TH1D* HistPE = GetHistForMjjj(ThisMass, &InFile, ipe);
        if (!HistPE) {
          std::cerr << "ERROR: cannot get PE: " << ipe << std::endl;
          exit(1);
        }

        FitObj MyFitObj;
        MyFitObj.Hist = HistPE;
        SetFitObjParams(MyFitObj);
        MyFitObj.gmean = ThisMass;
        MyFitObj.gsigmaRange = GetGausWidthRange(ThisMass);
        MyFitObj.IsData = false;
        MyFitObj.Section = Section;
        MyFitObj.ipe = ipe;

        Limit = LimitAtMass(MyFitObj);
        printf("PELimit %9i %9i %12.4f\n", ipe, (int) ThisMass, Limit);
        fprintf(OutFile, "%10E ", Limit);
        GausMeanNGaus.push_back(std::make_pair<float, float>(ThisMass, Limit));
      }
      fprintf(OutFile, "\n");
      fflush(OutFile);
      char FitHistName[100];
      sprintf(FitHistName, "Limit95NEvents_%i.eps", ipe);
      MakeGraph (GausMeanNGaus, "95\% CL Number of signal events", "Gaus mean [GeV]", "95\% CL Number of signal events", FitHistName);

    }
  }

  fclose(OutFile);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [Section] [InFileName]" << std::endl;
    return 1;
  }

  int const Section = atoi(argv[1]);
  TString const InFileName = argv[2];

  RunMultiJetLimits(Section, InFileName);

  return 0;
}
