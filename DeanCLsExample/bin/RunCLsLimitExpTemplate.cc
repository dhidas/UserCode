////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Sep 26 11:11:46 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "StandardHypoTestInvDemo.h"

#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TTree.h"
#include "TRandom.h"


#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooRandom.h"


float const BGSlope = -0.003;
float const XMin =    0;
float const XMax = 2000;
int   const NBins = 100;

TH1F* GetFakeData (int const NBackground, int const NSignal, float const SignalMass, TString const Name)
{
  // Make a fake data hist with or without some signal
  // Data is expo dist and signal gaussian

  // Expo Function - BG
  TF1 fMyExpo("g2", "[0]*TMath::Exp([1]*x)", XMin, XMax);
  fMyExpo.SetParameter(0, 1);
  fMyExpo.SetParameter(1, BGSlope);

  // Gaus function - Signal
  TF1 fMyGaus("MyGaus", "TMath::Gaus(x, [0], [0]*0.10)");
  fMyGaus.SetParameter(0, SignalMass);

  // Data histogram
  TH1F hData("MyFakeData", "MyFakeData", NBins, XMin, XMax);

  // Fill wish signal and background
  hData.FillRandom("g2", NBackground);
  hData.FillRandom("MyGaus", NSignal);

  // Save plot as eps
  TCanvas C;
  C.cd();
  hData.Draw("ep");
  C.SaveAs(Name + TString::Format("_%i.eps", (int) SignalMass));

  return (TH1F*) hData.Clone(Name);
}




int RunCLsLimit (float const TestMass, TH1F* DataHist, TH1F* MCHist_Background, TH1F* MCHist_Signal)
{
  // Grab data histogram
  if ( !(DataHist && MCHist_Signal && MCHist_Background) ) {
    std::cerr << "ERROR: cannot get data histogram" << std::endl;
    throw;
  }


  // Start a workspace
  RooWorkspace ws("ws");

  // Obseravable
  ws.factory("x[0]");
  ws.var("x")->setRange(XMin, XMax);

  // Import the data to Roo
  RooDataHist Data("Data", "dataset with x", *ws.var("x"), DataHist);
  RooDataHist Signal("signal", "dataset with x", *ws.var("x"), MCHist_Signal);
  RooDataHist Background("background", "dataset with x", *ws.var("x"), MCHist_Background);

  ws.import(Data);
  ws.import(Signal);
  ws.import(Background);


  // Parameter of Interest
  ws.factory("nsig[0, 2000]");
  ws.defineSet("POI","nsig");
  ws.factory("nbkg[0, 60000]");

  RooHistPdf SignalPdf("SignalPdf", "SignalPdf", *ws.var("x"), Signal, 2);
  RooHistPdf BackgroundPdf("BackgroundPdf", "BackgroundPdf", *ws.var("x"), Background, 2);
  ws.import(SignalPdf);
  ws.import(BackgroundPdf);



  // Define models (sig+bg, and bg only)
  ws.factory("SUM::sbmodel_noprior(nsig*SignalPdf, nbkg*BackgroundPdf)");
  ws.factory("SUM::bgmodel_noprior(nbkg*BackgroundPdf)");


  // This is here for studying systematics
  // Pick which constraints and nuisance params you want to use
  switch (0) {
    case 0:
      ws.factory("RooUniform::constraints(x)");
      ws.defineSet("nuisance","");
      break;
    default:
      std::cerr << "ERROR: Incorrect choice for systematicsi" << std::endl;
      exit(1);
      break;
  }



  // Build modelconfig for sbmodel and set current workspace
  RooStats::ModelConfig ModelConfigSB("ModelConfigSB");
  ModelConfigSB.SetWorkspace(ws);

  // Setup this model
  ModelConfigSB.SetPdf(*ws.pdf("sbmodel_noprior"));
  ModelConfigSB.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigSB.SetObservables(*ws.var("x"));
  ModelConfigSB.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigSB.SetNuisanceParameters(*ws.set("nuisance"));

  // Fit model and add POI snapsnot, import this modelconfig to the workspace
  ws.pdf("sbmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(XMin, XMax), RooFit::Extended(kTRUE));
  RooArgSet POIAndNuisSB("POIAndNuisSB");
  POIAndNuisSB.add(*ModelConfigSB.GetParametersOfInterest());
  ModelConfigSB.SetSnapshot(POIAndNuisSB);
  ModelConfigSB.SetGlobalObservables( RooArgSet() );
  ws.import(ModelConfigSB);


  // Build modelconfig for bgmodel and set current workspace
  RooStats::ModelConfig ModelConfigBG("ModelConfigBG");
  ModelConfigBG.SetWorkspace(ws);

  // Setup this model
  ModelConfigBG.SetPdf(*ws.pdf("bgmodel_noprior"));
  ModelConfigBG.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigBG.SetObservables(*ws.var("x"));
  ModelConfigBG.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigBG.SetNuisanceParameters(*ws.set("nuisance"));

  // Set the observable to zero and add POI snapsnot, import this modelconfig to the workspace
  ws.var("nsig")->setVal(0);
  RooArgSet POIAndNuisBG("POIAndNuisBG");
  POIAndNuisBG.add(*ModelConfigBG.GetParametersOfInterest());
  ModelConfigBG.SetSnapshot(POIAndNuisBG);
  ModelConfigBG.SetGlobalObservables( RooArgSet() );
  ws.import(ModelConfigBG);


  // Parameters of the CLs method we'll call
  int   const calculatorType    = 2;
  int   const testStatType      = 3;
  bool  const useCls            = true;
  int   const npoints           = 500;
  float const poimin            = 0;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
  float const poimax            = 400;
  int   const ntoys             = 50;
  bool  const useNumberCounting = false;
  const char* nuisPriorName     = "";

  // For debugging, why not save the workspace
  ws.SaveAs( TString::Format("Workspace_Exp_%i.root", (int) TestMass));

  // Run the actual CLs
  RooStats::HypoTestInvTool HTIT;
  HTIT.SetParameter("GenerateBinned", true);
  RooStats::HypoTestInverterResult* MyResult = HTIT.RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, useCls, npoints, poimin, poimax, ntoys, useNumberCounting, nuisPriorName);

  // Number of entries in result
  const int NEntries = MyResult->ArraySize();

  // Just some names
  const char* TypeName = "Hybrid";
  const char* ResultName = MyResult->GetName();
  TString PlotTitle = TString::Format("%s CL Scan for workspace %s", TypeName, ResultName);

  // Grab the result plot
  RooStats::HypoTestInverterPlot *Plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", PlotTitle, MyResult);
  TCanvas CanCLb;
  CanCLb.cd();
  Plot->Draw("CLb 2CL");  // plot all and Clb
  CanCLb.SaveAs( TString::Format("CLb2L_Exp_%i.eps", (int) TestMass));

  // Does not work with newer root version.  To be figured out later.
  if (false) {
    TCanvas Can;
    Can.Divide(2, (int) TMath::Ceil(NEntries/2));
    for (int i = 0; i < NEntries; ++i) {
      Can.cd(i + 1);
      RooStats::SamplingDistPlot * SamplingPlot = Plot->MakeTestStatPlot(i);
      SamplingPlot->SetLogYaxis(true);
      SamplingPlot->Draw();
      delete SamplingPlot;
    }
    Can.SaveAs( TString::Format("HTI_Result_Exp_%i.eps", (int) TestMass));
  }

  printf(" expected limit (-2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-2));
  printf(" expected limit (-1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-1));
  printf(" expected limit (median) %12.3E\n", MyResult->GetExpectedUpperLimit(0) );
  printf(" expected limit (+1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(1) );
  printf(" expected limit (+2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(2) );
  printf(" observed limit          %12.3E +/- %12.3E\n", MyResult->UpperLimit(), MyResult->UpperLimitEstimatedError()); 


  TString const FileName = TString::Format("Limits_%i.dat", (int) TestMass);
  FILE* f = fopen(FileName.Data(), "w");
  if (!f) {
    std::cerr << "ERROR: cannot open file for writing: " << FileName << std::endl;
    exit(1);
  }
  fprintf(f, "%15.3E\n", TestMass);
  fprintf(f, "%15.3E\n", MyResult->GetExpectedUpperLimit(-2));
  fprintf(f, "%15.3E\n", MyResult->GetExpectedUpperLimit(-1));
  fprintf(f, "%15.3E\n", MyResult->GetExpectedUpperLimit( 0));
  fprintf(f, "%15.3E\n", MyResult->GetExpectedUpperLimit( 1));
  fprintf(f, "%15.3E\n", MyResult->GetExpectedUpperLimit( 2));
  fprintf(f, "%15.3E\n", MyResult->UpperLimit());
  fclose(f);


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << "" << std::endl;
    return 1;
  }

  // Just in case you want to set the random seed
  RooRandom::randomGenerator()->SetSeed(27);
  gRandom->SetSeed(29);

  // What mass should we put in the data?
  float const SignalMass = 900;


  // You only get one data hist and one BG hist!
  TH1F* DataHist = GetFakeData(60000, 150, SignalMass, "Data");
  TH1F* MCHist_Background = GetFakeData(10000000, 0, 0, "Background");

  for (float mass = 200; mass <= 1800; mass += 50.) {
    // Generate the signal model for this mass
    TH1F* MCHist_Signal = GetFakeData(0, 10000000, mass, "Signal");

    // Run the limit code for this data, bg model and signal model
    RunCLsLimit(mass, DataHist, MCHist_Background, MCHist_Signal);
  }

  return 0;
}
