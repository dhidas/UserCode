///////////////////////////////////////////////////////////////////
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
#include "RooDataHist.h"
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

#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/FrequentistCalculator.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "StandardHypoTestInvDemo.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"



int KeithRooStats (TString const InFileName)
{

  using namespace RooStats;
  // Setup workspace
  RooWorkspace ws("ws");

  // Observable

  //ws.factory("JetPt[507.,1684.]");
  ws.factory("JetPt[507.,2000.]");
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
  
  TH1F* DataTH1F = (TH1F*) InFile.Get("PtcSpectrum_Ratio");
  //TH1F* DataTH1F = (TH1F*) InFile.Get("PtcSpectrum_ShiftDown_Ratio");
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

  // Signal
  ws.factory("CutOffPt[507.,2000.]");
  ws.factory("CutOffSpeed[10., 200.]");
  //ws.factory("NLO[0.3, 0.5]");
  //ws.factory("Height[0.5]");
  ws.factory("MyMult[0.95]");
  ws.factory("MyOffset[0.05]");
  ws.factory("RooGenericPdf::signalmodel('MyMult / (1.0 + TMath::Exp((JetPt - CutOffPt) / CutOffSpeed)) + MyOffset', {MyMult, JetPt, CutOffPt, CutOffSpeed, MyOffset})");
  ws.factory("RooGenericPdf::standardmodel('MyMult / (1.0 + TMath::Exp((JetPt - 2500.) / 100.)) + MyOffset', {MyMult, JetPt, MyOffset})");

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
  cout << "here" << endl;
  // standard model fit
  datafit = ws.var("JetPt")->frame();
  ws.data("DataToFit")->plotOn(datafit);
  ws.pdf("standardmodel")->fitTo(*ws.data("DataToFit"));
  ws.pdf("standardmodel")->plotOn(datafit);
  datafit->Draw();
  Can.SaveAs("StandardModelFit.eps");

	cout << "here two" << endl;
  // Make the model configuration
  RooStats::ModelConfig modelConfig("modelConfig");
  modelConfig.SetWorkspace(ws);
  modelConfig.SetPdf(*ws.pdf("signalmodel"));
  modelConfig.SetObservables(*ws.var("JetPt"));
  
  RooStats::ModelConfig BModel("modelConfigB");
  BModel.SetWorkspace(ws);
  BModel.SetPdf(*ws.pdf("standardmodel"));
  BModel.SetObservables(*ws.var("JetPt"));
  ws.defineSet("POI","CutOffPt");
  modelConfig.SetParametersOfInterest(*ws.set("POI"));
  BModel.SetParametersOfInterest(*ws.set("POI"));
  RooArgSet POIAndNuisS("POIAndNuisS");
  POIAndNuisS.add(*modelConfig.GetParametersOfInterest());

  RooArgSet POIAndNuisB("POIAndNuisB");
  POIAndNuisB.add(*BModel.GetParametersOfInterest());
  BModel.SetSnapshot(POIAndNuisB);

  ws.factory("RooUniform::prior_CutOffPt(CutOffPt)");
  ws.factory("RooUniform::prior_CutOffSpeed(CutOffSpeed)");
  //cout << "here" << endl;  
  ws.factory("RooPolynomial::SystS0(JetPt, {1.0,-0.003677, .000005819, -.00000000404, .000000000001082})");
  ws.factory("RooLognormal::Syst_prior(Syst[1], SystM0[0], SystS0)");
  ws.var("Syst")->setRange(0., 2.);
  ws.var("Syst")->setVal(1.);
  ws.var("SystM0")->setVal(1.);
  //cout << "here too" << endl;
  ws.factory("RooPolynomial::TheoS0(JetPt, {-.02508, 0.0003124, -0.00000003827, -0.000000000115, 0.00000000000006674})");
  ws.factory("RooLognormal::Theo_prior(Theo[1], TheoM0[0], TheoS0)");
  ws.var("Theo")->setRange(0., 2.);
  ws.var("Theo")->setVal(1.);
  ws.var("TheoM0")->setVal(1.);
  //cout << "here three" << endl;
  ws.factory("PROD::priors(prior_CutOffPt,prior_CutOffSpeed)");
  ws.factory("PROD::priors_withSys(priors,Syst_prior)");
  ws.factory("PROD::priors_withAll(priors_withSys,Theo_prior)");
  //cout << "here four" << endl;
  ws.factory("PROD::signalmodelWprior(signalmodel,priors_withAll)");
  ws.defineSet("nuisSet","CutOffSpeed,Theo_prior,Syst_prior");
  modelConfig.SetPriorPdf(*ws.pdf("priors_withAll"));
  BModel.SetPriorPdf(*ws.pdf("prior_CutOffPt"));
  //cout << "here five" << endl;


  RooArgSet nuisPar(*ws.var("CutOffSpeed"),*ws.var("Theo"),*ws.var("Syst"));
  //RooArgSet nuisPar(*ws.var("CutOffSpeed"));
	
  //cout << "here six" << endl;
  modelConfig.SetNuisanceParameters(nuisPar);
  POIAndNuisS.add(*modelConfig.GetNuisanceParameters());
  modelConfig.SetSnapshot(POIAndNuisS);
  BModel.SetNuisanceParameters(nuisPar);
  ws.import(modelConfig);
  ws.import(BModel);

  //nuisSet = &nuisPar;
  ws.Print();
  modelConfig.Print();
  BModel.Print();

  RooStats::MCMCCalculator mc(*ws.data("DataToFit"), modelConfig);
  mc.SetNumBins(50);
  mc.SetConfidenceLevel(0.95);
  mc.SetLeftSideTailFraction(1.0);
  //mc.SetNumIters(10000000);
  mc.SetNumIters(50000);
  mc.SetNumBurnInSteps(500);
  RooStats::MCMCInterval* interval = (RooStats::MCMCInterval*)mc.GetInterval();

  RooStats::BayesianCalculator bcalc(*ws.data("DataToFit"), modelConfig);
  double conflevel = 0.95; 
  bcalc.SetTestSize((1.-conflevel));

  RooStats::SimpleInterval* Bayes_interval = bcalc.GetInterval();
  
  // Upper limit

  //float LowerLimit = interval->LowerLimit(*ws.var("JetPt"));

  //std::cout << "LowerLimit: " << LowerLimit << std::endl;
  std::cout << "Bayesian Lower Limit:" << Bayes_interval->LowerLimit();
  TCanvas CanPost("MCMCpost", "MCMCpost");

  RooPlot *CLplot = bcalc.GetPosteriorPlot();
  CLplot->SetTitle("CMS Preliminary");
  CLplot->Draw();
  CanPost.SaveAs("BayesPosterior.eps");
  CanPost.SaveAs("BayesPosterior.pdf");



  HypoTestInverterResult * r = 0;
  const char * modelSBName = "modelConfig";
  const char * modelBName = "modelConfigB";
  const char * dataName = "DataToFit";
  int calculatorType = 1;
  int testStatType = 3;
  bool useCls = true;
  int npoints = 30;
  double poimin = 900;
  double poimax = 1940;
  int ntoys = 40;
  bool useNumberCounting = false;
  const char * nuisPriorName = "nuisSet";

  r = RunInverter(&ws, "modelConfig", "modelConfigB", "DataToFit", calculatorType, testStatType, npoints, poimin, poimax, ntoys, useCls, useNumberCounting, nuisPriorName );   
  //int nTestPlots = r->ArraySize();
  double CLsLowerLimit = r->LowerLimit();
  cout << "CLs Lower Limit:" << CLsLowerLimit << endl;
/*
   // compute expected limit
   std::cout << " expected limit (median) " <<  r->GetExpectedLowerLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedLowerLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedLowerLimit(1) << std::endl;
   std::cout << " expected limit (-2 sig) " << r->GetExpectedLowerLimit(-2) << std::endl;
   std::cout << " expected limit (+2 sig) " << r->GetExpectedLowerLimit(2) << std::endl;*/
  RooStats::MCMCIntervalPlot mcPlot(*interval);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.15); // X position of the title box from left
  gStyle->SetTitleY(.975); // Y position of the title box from bottom
  gStyle->SetLabelSize(0.035,"y");
  gStyle->SetStatX(.95);
  gStyle->SetStatY(.95);
  gStyle->SetStatW(0.20);
  gStyle->SetStatFontSize(0.044);
  gStyle->SetStatColor(0);
// For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

// For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

// For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  HypoTestInverterPlot CLS_plot("HTI_Result_Plot","CLs Result",r);
  CanPost.cd();
  //CLS_plot->Draw("CLs");
  TGraphErrors *CLS_graph = CLS_plot.MakePlot("CLsplusb");
  CLS_graph->Draw("ALP");
  CLS_graph->SetMinimum(.4);
  CLS_graph->SetMaximum(1.1);
  CLS_graph->GetXaxis()->SetTitle("Cutoff p_{T} (GeV)");
  TGraphErrors *CLS_graph2 = CLS_plot.MakePlot("CLb");
  CLS_graph2->Draw("ALP");
  CLS_graph2->SetMinimum(.4);
  CLS_graph2->SetMaximum(1.1);
  CLS_graph2->GetXaxis()->SetTitle("Cutoff p_{T} (GeV)");
  TGraphErrors *CLS_graph3 = CLS_plot.MakePlot("CLs");
  CLS_graph3->Draw("ALP");
  CLS_graph3->SetMinimum(.4);
  CLS_graph3->SetMaximum(1.1);
  CLS_graph3->GetXaxis()->SetTitle("Cutoff p_{T} (GeV)");
  TGraphErrors *CLS_graph4 = CLS_plot.MakePlot("2CL");
  CLS_graph4->Draw("ALP");
  CLS_graph4->SetMinimum(.4);
  CLS_graph4->SetMaximum(1.1); 
  CLS_graph4->GetXaxis()->SetTitle("Cutoff p_{T} (GeV)");
  TFile OutFile("LimitPlots.root", "RECREATE");
  OutFile.cd();

  CLS_graph->Write();
  CLS_graph2->Write();
  CLS_graph3->Write();
  CLS_graph4->Write();

  CanPost.SaveAs("CLs_limit.pdf");
  CanPost.SetBorderSize(0);
  CanPost.SetBorderMode(0);

/*
  for(int ti = 0; ti < nTestPlots; ti++){
	SamplingDistPlot *pl = CLS_plot.MakeTestStatPlot(ti);
	pl->SetLogYaxis(true);
	pl->Write();
  }


  mcPlot.SetTitle("CMS Preliminary");
  mcPlot.Draw();
  mcPlot.Write();  
  CanPost.SaveAs("MCMCPosterior.eps");
  CanPost.SaveAs("MCMCPosterior.pdf");
*/
  OutFile.Write(); 
  OutFile.Close();

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
