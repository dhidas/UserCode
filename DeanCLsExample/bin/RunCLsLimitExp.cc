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


#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"


float const BGSlope = -0.003;
float const XMin =    0;
float const XMax = 2000;
int   const NBins = 200;

TH1F* GetFakeData (int const NBackground, int const NSignal, float const SignalMass)
{
  // Make a fake data hist with or without some signal
  // Data is expo dist and signal gaussian

  // Open output root file
  TFile RootFile("FakeDataExp.root", "recreate");

  // Expo Function - BG
  TF1 fMyExpo("g2", "[0]*TMath::Exp([1]*x)", XMin, XMax);
  fMyExpo.SetParameter(0, 1);
  fMyExpo.SetParameter(1, BGSlope);

  // Gaus function - Signal
  TF1 fMyGaus("MyGaus", "TMath::Gaus(x, [0], [0]*0.065)");
  fMyGaus.SetParameter(0, SignalMass);

  // Data histogram
  TH1F hData("MyFakeData", "MyFakeData", NBins, XMin, XMax);

  // Fill wish signal and background
  hData.FillRandom("g2", NBackground);
  hData.FillRandom("MyGaus", NSignal);
  hData.Fit("g2");

  // Save plot as eps
  TCanvas C;
  C.cd();
  hData.Draw();
  fMyExpo.Draw("same");
  C.SaveAs("FakeDataExp.eps");

  TTree Tree("FakeDataTree", "FakeDataTree");
  Tree.SetDirectory(&RootFile);
  float Value;
  Tree.Branch("x", &Value);
  for (int i = 0; i != NBackground; ++i) {
    Value = fMyExpo.GetRandom();
    Tree.Fill();
  }
  for (int i = 0; i != NSignal; ++i) {
    Value = fMyGaus.GetRandom();
    Tree.Fill();
  }

  // write and close file
  RootFile.Write();
  RootFile.Close();


  return (TH1F*) hData.Clone(0);
}




int RunCLsLimit ()
{
  // What mass for signal to consider?
  float const SignalMass = 900;

  // Grab data histogram
  TH1F* DataHist = GetFakeData(15000, 0, SignalMass);
  if (!DataHist) {
    std::cerr << "ERROR: cannot get data histogram" << std::endl;
    throw;
  }

  float const NBackground = DataHist->Integral( DataHist->FindBin(XMin), DataHist->FindBin(XMax) );
  std::cout << "NBackground: " << NBackground << std::endl;


  // Start a workspace
  RooWorkspace ws("ws");

  // Obseravable
  ws.factory("x[0]");
  ws.var("x")->setRange(XMin, XMax);

  // Import the data to Roo
  RooDataHist Data("Data", "dataset with x", *ws.var("x"), DataHist);

  // If instead you want to use the unbinned data you can uncomment below and comment out above
  //TFile InRootFile("FakeDataExp.root", "read");
  //TTree* Tree = (TTree*) InRootFile.Get("FakeDataTree");
  //std::cout << Tree->GetEntries() << std::endl;
  //RooDataSet Data("Data", "dataset with x", *ws.var("x"), RooFit::Import(*Tree));
  ws.import(Data);




  // Define background and parameters
  ws.factory("p1[0,1]");
  ws.var("p1")->setVal(-BGSlope);
  ws.var("p1")->setConstant(kTRUE);
  ws.factory("RooGenericPdf::background('exp(-p1*x)', {x, p1})");


  // Parameter of Interest
  ws.factory("nsig[0, 2000]");
  ws.defineSet("POI","nsig");


  // Signal shape
  ws.factory("RooGaussian::signal(x, sigMean[0], sigWidth[0,500])");
  ws.var("sigWidth")->setVal(SignalMass*0.065);
  ws.var("sigWidth")->setRange( (SignalMass*0.065)*(1.0-0.10), (SignalMass*0.065)*(1.0+0.10));
  ws.var("sigWidth")->setConstant(true);
  ws.var("sigMean")->setVal(SignalMass);
  ws.var("sigMean")->setConstant(true);

  // prior on the signal width
  ws.factory("RooUniform::sigWidth_prior(sigWidth)");


  // nbkg prior
  ws.factory("RooLognormal::nbkg_prior(nbkg[0,60000], nbkgM0[0], nbkgS0[0])");
  ws.var("nbkgS0")->setVal(1. + 0.03);  // 3% error
  ws.var("nbkgS0")->setConstant(true);
  ws.var("nbkgM0")->setVal(NBackground);
  ws.var("nbkgM0")->setConstant(true);

  ws.var("nbkg")->setVal(NBackground);
  ws.var("nbkg")->setRange(NBackground * (1 - 3*0.03), NBackground * (1 + 3*0.03));
  ws.var("nbkg")->setConstant(true);

  // Define a prior for this p1 parameter
  ws.factory("RooLognormal::p1_prior(p1, p1M0[0], p1S0[1])");
  ws.var("p1S0")->setVal(1.05);  // 5% error
  ws.var("p1S0")->setConstant(true);
  ws.var("p1M0")->setVal(-BGSlope); // mean val
  ws.var("p1M0")->setConstant(true);


  // Define models (sig+bg, and bg only)
  ws.factory("SUM::sbmodel_noprior(nsig*signal, nbkg*background)");
  ws.factory("SUM::bgmodel_noprior(nbkg*background)");


  // Pick which constraints and nuisance params you want to use
  switch (1) {
    case 0:
      ws.factory("RooUniform::constraints(x)");
      ws.defineSet("nuisance","");
      break;
    case 1:
      ws.factory("PROD::constraints(nbkg_prior)");
      ws.defineSet("nuisance", "nbkg");
      ws.var("nbkg")->setConstant(false);
      break;
    case 2:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior)");
      ws.defineSet("nuisance", "nbkg,p1");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      break;
    case 3:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,sigWidth_prior)");
      ws.defineSet("nuisance", "nbkg,p1,sigWidth");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("sigWidth")->setConstant(false);
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
  int   const npoints           = 20;
  float const poimin            = 0;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
  float const poimax            = 100;
  int   const ntoys             = 400;
  bool  const useNumberCounting = false;
  const char* nuisPriorName     = "";

  // For debugging, why not save the workspace
  ws.SaveAs("Workspace_Exp.root");

  // Run the actual CLs
  RooStats::HypoTestInvTool HTIT;
  HTIT.SetParameter("GenerateBinned", true);
  RooStats::HypoTestInverterResult* MyResult = HTIT.RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, useCls, npoints, poimin, poimax, ntoys, useNumberCounting, nuisPriorName);
  //RooStats::HypoTestInverterResult* MyResult = RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, npoints, poimin, poimax, ntoys, useCls, useNumberCounting, nuisPriorName);

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
  CanCLb.SaveAs("CLb2L_Exp.eps");

  if (true) {
    TCanvas Can;
    Can.Divide(2, TMath::Ceil(NEntries/2));
    for (int i = 0; i < NEntries; ++i) {
      Can.cd(i + 1);
      RooStats::SamplingDistPlot * SamplingPlot = Plot->MakeTestStatPlot(i);
      SamplingPlot->SetLogYaxis(true);
      SamplingPlot->Draw();
      delete SamplingPlot;
    }
    Can.SaveAs("HTI_Result_Exp.eps");
  }

  printf(" expected limit (-2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-2));
  printf(" expected limit (-1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-1));
  printf(" expected limit (median) %12.3E\n", MyResult->GetExpectedUpperLimit(0) );
  printf(" expected limit (+1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(1) );
  printf(" expected limit (+2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(2) );

  printf(" observed limit          %12.3E +/- %12.3E\n", MyResult->UpperLimit(), MyResult->UpperLimitEstimatedError()); 



  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << "" << std::endl;
    return 1;
  }

  RunCLsLimit();

  return 0;
}
