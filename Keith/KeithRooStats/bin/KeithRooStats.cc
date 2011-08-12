////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Tue Jul 19 16:41:16 CEST 2011
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
#include "TLine.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TStyle.h"

#include "RooWorkspace.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "RooBinning.h"
#include "RooLognormal.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooRandom.h"

#include "RooStats/BayesianCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/UpperLimitMCSModule.h"


int KeithRooStats (TString const InFileName)
{
  // Setup workspace
  RooWorkspace ws("ws");

  // Observable

  //ws.factory("JetPt[507.,1684.]");
  ws.factory("JetPt[0.,3000.]");
  RooRealVar* JetPt = ws.var("JetPt");
  //ws.factory("RooGenericPdf::CutOffPt(-568. + 370*Scale)");  
  //ws.factory("JetPt[707.,1684.]");
  //RooRealVar* CutOffPt = ws.var("CutOffPt");
  //RooRealVar* JetPt = new RooRealVar("JetPt", "JetPt", CutOffPt, 1684.);

  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  //TH1F* DataTH1F = (TH1F*) InFile.Get("PtSpectrum_Ratio");
  TH1F* DataTH1F = (TH1F*) InFile.Get("BSignal");
  //TH1F* DataTH1F = (TH1F*) InFile.Get("Extinction2");
  // TH1F* DataTH1FDen = (TH1F*) InFile.Get("Pt_Ybin_0");
  if (!DataTH1F) {
    std::cerr << "ERROR: cannot get histogram!!" << std::endl;
    throw;
  }
  TH1* DataTH1 = (TH1*) DataTH1F;
  RooDataHist binnedData("data", "dataset with x", *ws.var("JetPt"), DataTH1F);
  // RooDataHist binnedDataDen("dataden", "dataset with x", *ws.var("JetPt"), DataTH1FDen);
  RooDataHist fitData("DataToFit", "dataset with x", *ws.var("JetPt"), RooFit::Import(*DataTH1, false));
  ws.import(fitData);

  // Just draw what we see in the data..
  TCanvas Can;
  Can.cd();
  RooPlot* DataPlot = ws.var("JetPt")->frame();
  ws.data("DataToFit")->plotOn(DataPlot);
  DataPlot->Draw();
  Can.SaveAs("DataPlot.eps");


  // standardmodel is a line
  ws.factory("RooUniform::standardmodel(JetPt)");

  // Signal
  ws.factory("CutOffPt[100.,3000.]");
  ws.factory("CutOffSpeed[10., 200]");
  //ws.factory("NLO[0.3, 0.5]");
  //ws.factory("Height[0.5]");
  ws.factory("MyMult[1]");
  ws.factory("MyOffset[0]");
  ws.factory("RooGenericPdf::signalmodel('MyMult / (1.0 + TMath::Exp((JetPt - CutOffPt) / CutOffSpeed)) + MyOffset', {MyMult, JetPt, CutOffPt, CutOffSpeed, MyOffset})");

  //ws.factory("RooGaussian::signalmodel(JetPt, CutOffPt, CutOffSpeed)");
  //ws.import(signalmodel);
  // model = line + cutoff function thing

  // signalmodel fit
  RooPlot* datafit = ws.var("JetPt")->frame();
  ws.data("DataToFit")->plotOn(datafit);
  //JetPt->setRange("ptrange", 507.,1684.);
  ws.pdf("signalmodel")->fitTo(*ws.data("DataToFit"),RooFit::SumW2Error(kTRUE));
  ws.pdf("signalmodel")->plotOn(datafit);
  datafit->Draw();
  Can.SaveAs("SignalModelFit.eps");
  Can.Clear();
  delete datafit;
  ws.var("CutOffPt")->Print();
  ws.var("CutOffSpeed")->Print();

  // standard model fit
  datafit = ws.var("JetPt")->frame();
  ws.data("DataToFit")->plotOn(datafit);
  ws.pdf("standardmodel")->fitTo(*ws.data("DataToFit"));
  ws.pdf("standardmodel")->plotOn(datafit);
  datafit->Draw();
  Can.SaveAs("StandardModelFit.eps");


  // Make the model configuration
  RooStats::ModelConfig modelConfig("modelConfig");
  modelConfig.SetWorkspace(ws);
  modelConfig.SetPdf(*ws.pdf("signalmodel"));
  ws.defineSet("POI","CutOffPt");
  modelConfig.SetParametersOfInterest(*ws.set("POI"));
  ws.factory("RooUniform::prior_CutOffPt(CutOffPt)");
  ws.factory("RooUniform::prior_CutOffSpeed(CutOffSpeed)");
  
  ws.factory("RooPolynomial::SystS0(JetPt, {-0.003677, .000005819, -.00000000404, .000000000001082})");
  ws.factory("RooLognormal::Syst_prior(Syst[1], SystM0[0], SystS0)");
  ws.var("Syst")->setRange(0., 2.);
  ws.var("Syst")->setVal(1.);
  ws.var("SystM0")->setVal(1.);

  ws.factory("RooPolynomial::TheoS0(JetPt, {-.004072, .000006684, -0.000000004592, .000000000001226})");
  ws.factory("RooLognormal::Theo_prior(Theo[1], TheoM0[0], TheoS0)");
  ws.var("Theo")->setRange(0., 2.);
  ws.var("Theo")->setVal(1.);
  ws.var("TheoM0")->setVal(1.);

  ws.factory("PROD::priors(prior_CutOffPt, prior_CutOffSpeed)");
  ws.factory("PROD::priors_withSys(priors, Syst_prior)");
  ws.factory("PROD::priors_withAll(priors_withSys, Theo_prior)");
  ws.factory("PROD::signalmodelWprior(signalmodel,priors_withAll)");
  ws.defineSet("nuisSet","CutOffSpeed,Syst_prior,Theo_prior");
  ws.import(modelConfig);

  ws.Print();
  modelConfig.Print();


  RooStats::MCMCCalculator mc(*ws.data("DataToFit"), modelConfig);
  mc.SetNumBins(50);
  mc.SetConfidenceLevel(0.95);
  mc.SetLeftSideTailFraction(1.0);
  //mc.SetNumIters(10000000);
  mc.SetNumIters(100000);
  mc.SetNumBurnInSteps(500);
  RooStats::MCMCInterval* interval = (RooStats::MCMCInterval*)mc.GetInterval();

  // Upper limit
  float LowerLimit = interval->LowerLimit(*ws.var("JetPt"));
  std::cout << "LowerLimit: " << LowerLimit << std::endl;

  TCanvas CanPost("MCMCpost", "MCMCpost");
  CanPost.cd();
  RooStats::MCMCIntervalPlot mcPlot(*interval);
  mcPlot.Draw();
  CanPost.SaveAs("MCMCPosterior.eps");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  KeithRooStats(InFileName);

  return 0;
}
