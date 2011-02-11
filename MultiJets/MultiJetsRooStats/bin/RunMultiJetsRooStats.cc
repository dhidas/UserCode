////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Jan 19 12:26:55 CET 2011
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

#include "TStopwatch.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"





// Utility functions
TH1D* GetHistForMjjj (float const Mjjj, TFile* File, TString const Suffix)
{
  // This function will return a histogram
  // Before, I had put all PEs in a file and ran over them.. now generate
  // on the go.. so this may or may not be used for more than data..

  char BUFF[200];
  if (Mjjj <= 250)        sprintf(BUFF,  "Mjjj_45_20_130_" + Suffix);
  else                    sprintf(BUFF,  "Mjjj_45_20_170_" + Suffix);
  //sprintf(BUFF,  "Mjjj_45_20_170_6jet");


  std::cout << "Getting Hist: " << BUFF << "  for mass " << Mjjj << std::endl;
  TH1D* Hist = (TH1D*) File->Get(BUFF);
  if (!Hist) {
    std::cerr << "ERROR: cannot get hist " << BUFF << std::endl;
    exit(1);
  }

  return Hist;
}



float GetAcceptanceForMjjj (float const Mjjj)
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



std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  // This will eventually be a parametrization

  if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  if (Mjjj < 350) return std::make_pair<float, float>(10 + 10 * (Mjjj - 250.) / 100., 15 + 10 * (Mjjj - 250.) / 100.);
  //if (Mjjj < 350) return std::make_pair<float, float>(15, 20);
  return std::make_pair<float, float>(20, 25);
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



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////






float RunMultiJetsRooStats (TString const InFileName, float const SignalMass, int const method, int const statLevel, int const Section)
{
  // Declare some constants
  float const MININVMASS =    0;
  float const MAXINVMASS = 1500;

  float const LUMINOSITY = 35.0;
  float const LUMIERROR  = 0.11;

  float const JESERROR   = 0.08;

  float const ACCERROR   = 0.13;

  float const MAXXS      = 1000;

  char label[100];
  sprintf(label, "Dean_%i_%i_%i", (int) SignalMass, method, statLevel);

  // Setup workspace
  RooWorkspace ws("ws");

  // Observable
  ws.factory("mjjj[0]");
  RooRealVar* mjjj = ws.var("mjjj");
  mjjj->setRange(MININVMASS, MAXINVMASS);

  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1* DataTH1 = (TH1*) GetHistForMjjj(SignalMass, &InFile, "6jet");
  RooDataHist binnedData("data", "dataset with x", *ws.var("mjjj"), DataTH1);
  ws.import(binnedData);
  TF1* Func = DataTH1->GetFunction("total");
  printf("Fit parameters: %12.3f %12.3f %12.3f %12.3f\n",
      Func->GetParameter(0),
      Func->GetParameter(1),
      Func->GetParameter(2),
      Func->GetParameter(0) * Func->GetParameter(2) / 10.0
      );
  float const NBGFromFit = Func->GetParameter(0) * Func->GetParameter(2) / 10.0;

  // Other fit, Used for syst numbers
  TH1* SystTH1 = (TH1*) GetHistForMjjj(SignalMass, &InFile, "4jetscaled");
  TF1* SystFunc = SystTH1->GetFunction("total");
  printf("Syst Numbers: %12.3f %12.3f %12.3f\n",
      SystFunc->GetParameter(0) / Func->GetParameter(0),
      SystFunc->GetParameter(1) / Func->GetParameter(1),
      SystFunc->GetParameter(2) / Func->GetParameter(2));
  TCanvas Can1;
  Can1.cd();
  DataTH1->Draw("ep");
  Func->Draw("samel");
  SystFunc->Draw("lsame");
  Can1.SaveAs(TString("DefaultFit_")+label+".gif");


  // background and background prior
  //ws.factory("CEXPR::background('TMath::Landau(mjjj, lmpv, lsigma) ', lmpv[200,0,1000], lsigma[40,0,500], mjjj)");
  ws.factory("RooLandau::background(mjjj, lmpv[200,0,1000], lsigma[40,0,500])");
  ws.factory("RooGaussian::nbkg_prior(nbkg[0,60000], nbkgM0[0], nbkgS0[0])");
  //ws.factory("RooUniform::nbkg_prior(nbkg[0,60000])");
  // Set the background parameters
  ws.var("nbkg")->setVal(NBGFromFit);
  ws.var("nbkg")->setRange(NBGFromFit * (1 - Func->GetParError(0)/Func->GetParameter(0)), NBGFromFit * (1 + Func->GetParError(0)/Func->GetParameter(0)));
  ws.var("nbkgM0")->setVal(NBGFromFit);
  ws.var("nbkgS0")->setVal(TMath::Sqrt(NBGFromFit));
  ws.var("lmpv")->setVal(Func->GetParameter(1));
  ws.var("lsigma")->setVal(Func->GetParameter(2));

  // This is for the signal
  ws.factory("RooGaussian::lumi_prior(lumi[0], lumiM0[0], lumiS0[0])");
  ws.var("lumi")->setVal(LUMINOSITY);
  ws.var("lumi")->setRange(LUMINOSITY * (1. - 3. * LUMIERROR), LUMINOSITY * (1. + 3. * LUMIERROR));
  ws.var("lumiM0")->setVal(LUMINOSITY);
  ws.var("lumiS0")->setVal(LUMINOSITY * LUMIERROR);
  ws.factory("RooUniform::xs_prior(xs[-5000,5000])");
  ws.factory("prod::nsig(xs, lumi, acceptance[0,1])");
  ws.var("acceptance")->setVal(GetAcceptanceForMjjj(SignalMass));
  ws.var("acceptance")->setConstant(true);
  ws.factory("RooGaussian::acceptance_prior(acceptanceDelta[1,0,2], acceptanceDeltaM0[1.0], acceptanceDeltaS0[0.1])");
  ws.var("acceptanceDelta")->setRange(1. - 3. * ACCERROR, 1. + 3. * ACCERROR);
  ws.var("acceptanceDeltaS0")->setVal(ACCERROR);

  ws.var("xs")->setRange(0, MAXXS);

  // Signal
  ws.factory("sigMass[1000]");
  ws.var("sigMass")->setRange(MININVMASS, MAXINVMASS);
  ws.var("sigMass")->setVal(SignalMass);
  ws.factory("RooGaussian::sigMassDelta_prior(sigMassDelta[1.0,0.0,2.0], sigMassDeltaM0[1.0], sigMassDeltaS0[0.1])");
  ws.var("sigMassDelta")->setRange(1. - 3. * JESERROR, 1. + 3. * JESERROR);
  ws.var("sigMassDelta")->setVal(1);
  ws.var("sigMassDeltaS0")->setVal(JESERROR);
  ws.factory("cexpr::invmassprime('sigMassDelta*mjjj',sigMassDelta,mjjj)");

  // Signal shape
  ws.factory("RooGaussian::signal(invmassprime, sigMass, sigWidth[0,50])");
  //ws.factory("RooGaussian::signal(mjjj, sigMass, sigWidth[0,50])");
  ws.var("sigWidth")->setVal( GetGausWidthRange(SignalMass).first/2. + GetGausWidthRange(SignalMass).second/2.);
  ws.var("sigWidth")->setRange( GetGausWidthRange(SignalMass).first, GetGausWidthRange(SignalMass).second);


  // Define model and parameters of interest and so on
  ws.factory("SUM::model(nsig*signal, nbkg*background)");
  ws.factory("SUM::background_model(nbkg*background)");
  ws.defineSet("observables","mjjj");
  ws.defineSet("POI","xs");
  ws.factory("PROD::prior0(xs_prior)");
  ws.factory("PROD::prior1(xs_prior,nbkg_prior)");
  ws.factory("PROD::prior2(xs_prior,lumi_prior)");
  ws.factory("PROD::prior3(xs_prior,acceptance_prior)");
  ws.factory("PROD::prior4(xs_prior,lumi_prior,nbkg_prior)");
  ws.factory("PROD::prior5(xs_prior,nbkg_prior,lumi_prior,sigMassDelta_prior)");
  ws.factory("PROD::prior6(xs_prior,nbkg_prior,lumi_prior,acceptance_prior)");
  ws.defineSet("nuisSet0","");
  ws.defineSet("nuisSet1","nbkg");
  ws.defineSet("nuisSet2","lumi");
  ws.defineSet("nuisSet3","acceptance");
  ws.defineSet("nuisSet4","lumi,nbkg");
  ws.defineSet("nuisSet5","nbkg,lumi,sigMassDelta");
  ws.defineSet("nuisSet6","nbkg,lumi,acceptance");

  ws.var("xs")->setConstant(false);

  ws.var("lumi")->setConstant(true);
  ws.var("acceptance")->setConstant(true);
  ws.var("lmpv")->setConstant(true);
  ws.var("lsigma")->setConstant(true);
  ws.var("nbkg")->setConstant(true);
  ws.var("sigWidth")->setConstant(true);
  ws.var("sigMassDelta")->setConstant(true);
  ws.var("sigMassDeltaM0")->setConstant(true);
  ws.var("sigMassDeltaS0")->setConstant(true);

  // setup model config
  RooStats::ModelConfig modelConfig("modelConfig");
  modelConfig.SetWorkspace(ws);
  modelConfig.SetPdf(*ws.pdf("model"));
  modelConfig.SetParametersOfInterest(*ws.set("POI"));

  if (statLevel == 0) {
    modelConfig.SetPriorPdf(*ws.pdf("prior0"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet0"));
  } else if (statLevel == 1) {
    ws.var("nbkg")->setConstant(false);
    modelConfig.SetPriorPdf(*ws.pdf("prior1"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet1"));
  } else if (statLevel == 2) {
    ws.var("lumi")->setConstant(false);
    modelConfig.SetPriorPdf(*ws.pdf("prior2"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet2"));
  } else if (statLevel == 3) {
    ws.var("acceptance")->setConstant(false);
    modelConfig.SetPriorPdf(*ws.pdf("prior3"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet3"));
  } else if (statLevel == 4) {
    ws.var("nbkg")->setConstant(false);
    ws.var("lumi")->setConstant(false);
    modelConfig.SetPriorPdf(*ws.pdf("prior4"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet4"));
  } else if (statLevel == 6) {
    ws.var("nbkg")->setConstant(false);
    ws.var("lumi")->setConstant(false);
    ws.var("acceptance")->setConstant(false);
    modelConfig.SetPriorPdf(*ws.pdf("prior6"));
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet6"));
  }
  ws.import(modelConfig);


  //RooFitResult* Fit = ws.pdf("model")->fitTo(binnedData);
  ws.var("nbkg")->Print();
  TCanvas Can;
  Can.cd();
  RooPlot* datafit = ws.var("mjjj")->frame();
  ws.data("data")->plotOn(datafit);
  ws.pdf("model")->fitTo(*ws.data("data"));
  ws.pdf("model")->plotOn(datafit);
  datafit->Draw();
  Can.SaveAs(TString("Fit_")+label+".gif");
  ws.var("sigMassDelta")->Print();
  //exit(0);

  ws.Print();
  modelConfig.Print();
  ws.var("nbkg")->Print();
  ws.var("nbkgM0")->Print();
  ws.var("nbkgS0")->Print();
  ws.var("lmpv")->Print();
  ws.var("lsigma")->Print();
  ws.var("lumi")->Print();
  ws.var("lumiM0")->Print();
  ws.var("lumiS0")->Print();
  ws.var("acceptance")->Print();
  ws.var("xs")->Print();
  ws.var("sigMass")->Print();
  ws.var("sigMassDeltaM0")->Print();
  ws.var("sigMassDeltaS0")->Print();
  ws.var("sigWidth")->Print();


  double lower=-1, upper=-1;
  if(method==0) {
    RooStats::BayesianCalculator bc(*ws.data("data"), modelConfig);
    bc.SetConfidenceLevel(0.95);
    bc.SetLeftSideTailFraction(0.0);
    RooStats::SimpleInterval* bInt = bc.GetInterval();
    TCanvas* canvas = new TCanvas("bcPlost","Posterior Distribution",500,500);
    RooPlot* plot = bc.GetPosteriorPlot();
    plot->Draw();
    canvas->SaveAs(TString("BCpost")+label+".gif");
    lower=bInt->LowerLimit();
    upper=bInt->UpperLimit();

  } else if(method==1) {
    //RooStats::ProposalHelper ph;
    //ph.SetVariables((RooArgSet&)fit->floatParsFinal());
    //ph.SetCovMatrix(fit->covarianceMatrix());
    //ph.SetUpdateProposalParameters(true);
    //ph.SetCacheSize(100);
    //RooStats::ProposalFunction* pdfProp = ph.GetProposalFunction();

    RooStats::MCMCCalculator mc(*ws.data("data"), modelConfig);
    mc.SetNumBins(100);
    mc.SetConfidenceLevel(0.95);
    mc.SetLeftSideTailFraction(0.0);
    mc.SetNumIters(500000);
    mc.SetNumBurnInSteps(500);
    //mc.SetProposalFunction(*pdfProp);
    RooStats::MCMCInterval* interval = (RooStats::MCMCInterval*)mc.GetInterval();

    // draw posterior plot
    TCanvas * c1 = new TCanvas("MCMCpost", "MCMCpost", 500, 500);
    RooStats::MCMCIntervalPlot mcPlot(*interval);
    mcPlot.Draw();
    c1->SaveAs(TString("mcmcPost")+label+".gif");
    lower=interval->LowerLimit(*ws.var("xs"));
    upper=interval->UpperLimit(*ws.var("xs"));

  } else if(method==2) {
    RooStats::ProfileLikelihoodCalculator plc(*ws.data("data"), modelConfig);
    plc.SetConfidenceLevel(0.95);
    RooStats::LikelihoodInterval* plInt=plc.GetInterval();
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    plInt->LowerLimit( *ws.var("xs") ); // get ugly print out of the way. Fix.
    RooMsgService::instance().setGlobalKillBelow(msglevel);

    // draw posterior plot
    TCanvas * c1 = new TCanvas("PLikePost", "PLikePost", 500, 500);
    RooStats::LikelihoodIntervalPlot* plotInt=new RooStats::LikelihoodIntervalPlot(plInt);
    plotInt->SetTitle("Profile Likelihood Ratio and Posterior for S");
    plotInt->Draw();
    c1->SaveAs(TString("PLikePost")+label+".gif");
    lower = plInt->LowerLimit(*ws.var("xs"));
    upper = plInt->UpperLimit(*ws.var("xs"));
  }

  printf("LOWER LIMIT: %12.3f\n", lower);
  printf("UPPER LIMIT: %12.3f\n", upper);

  std::cout << "BEFORE" << std::endl;
  // Setup some stuff for toy MC experiments TOY, everyone likes TOYS, but not to be toyed with.
  // I can get you a toy, believe me.  There are ways, Dude.
  RooMCStudy* MCS = 0x0;
  if (Section >= 0) {
    MCS = new RooMCStudy(*ws.pdf("background_model"), *ws.var("mjjj"), RooFit::Extended(kTRUE), RooFit::FitModel(*ws.pdf("model")), RooFit::FitOptions( RooFit::Extended(kTRUE), RooFit::PrintEvalErrors(1)), RooFit::Constrain(*modelConfig.GetNuisanceParameters()), RooFit::Silence());
    //MCS = new RooMCStudy(*ws.pdf("background"), *ws.var("mjjj"), RooFit::Extended(kTRUE), RooFit::FitModel(*ws.pdf("model")), RooFit::FitOptions( RooFit::Extended(kTRUE), RooFit::PrintEvalErrors(1)));
    //MCS = new RooMCStudy(*ws.pdf("model"), *ws.var("mjjj"), RooFit::Extended(kTRUE), RooFit::FitOptions( RooFit::Extended(kTRUE), RooFit::PrintEvalErrors(-1)));
    //MCS = new RooMCStudy(*ws.pdf("model"), *ws.var("mjjj"), RooFit::Extended(kTRUE), RooFit::FitOptions( RooFit::Extended(kTRUE), RooFit::PrintEvalErrors(-1)), RooFit::Constrain(*modelConfig.GetNuisanceParameters()));
    //MCS = new RooMCStudy(*ws.pdf("background"), *ws.var("mjjj"), RooFit::Extended(kTRUE), RooFit::FitModel(*ws.pdf("model")), RooFit::FitOptions( RooFit::Extended(kTRUE), RooFit::PrintEvalErrors(1)) );
  }

  //RooStats::UpperLimitMCSModule ULMCS(ws.set("POI"), 0.95);
  //RooDLLSignificanceMCSModule sigModule("xs", 0); // For significance test??
  //MCS->addModule(ULMCS);
  MCS->generateAndFit(2000);
  MCS->Print();
  std::cout << "AFTER" << std::endl;

  TString limitstr("xs");
  TH1* mcslimit_histo = (TH1F*) MCS->fitParDataSet().createHistogram(limitstr+ws.set("POI")->GetName());
  //Make a histogram of the upperlimit vs the number of generated events
  //TString limitstr2("ngen,xs");
  //TH1* mcslimitvsevt_histo = (TH1F*) MCS->fitParDataSet().createHistogram(limitstr2+ws.set("POI")->GetName());
  TCanvas* c2 =new TCanvas();
  c2->cd();
  c2->Divide(1,2);
  c2->cd(1);
  RooPlot* hist = MCS->plotParam(*ws.var("xs"));
  std::cout << "XXX " << hist->GetNbinsX() << std::endl;
  hist->Draw();
  c2->cd(2);
  mcslimit_histo->Draw();
  //mcslimitvsevt_histo->Draw();
  //c2->Draw();
  c2->SaveAs("plplpl.eps");
  exit(0);

  return upper;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  int const Section = 0;

  //TString const InFileName = "/Users/dhidas/Data35pb/LandauFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  TString const InFileName = "/home/dhidas/Data35pb/LandauFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";

  float const BeginMass = 200;
  float const EndMass   = 500;
  float const StepSize  =  10;

  int const Method      = 0;
  int const Systematics = 0;

  if (Section == -1) {
    std::vector< std::pair<float, float> > MassLimitVec;
    for (float Mass = BeginMass; Mass <= EndMass; Mass += StepSize) {
      MassLimitVec.push_back( std::make_pair<float, float>(Mass, RunMultiJetsRooStats(InFileName, Mass, Method, Systematics, -1)) );
    }
    MakeGraph(MassLimitVec, "95\% C.L.", "M_{jjj}", "95\% C.L. #sigma (pb)", "LimitPlot.gif");
  } else {
    for (float Mass = BeginMass; Mass <= EndMass; Mass += StepSize) {
      RunMultiJetsRooStats(InFileName, Mass, Method, Systematics, Section);
    }
  }

  return 0;
}