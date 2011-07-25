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

#include "TStopwatch.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"





// Utility functions


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



TH1F* GetPEExpoWithBump (int const NBG, float const cexp, int const NSig, float const gMean, float const gSigma, TString const label = "blah")
{
  TH1F* h = new TH1F("PE", "PE", 65, 170, 800);

  TF1 f("myexpgaus", "[0]*TMath::Exp([1]*x) + [2]*TMath::Gaus(x, [3], [4], 1)", 170, 800);
  // The order of this is very specific...
  f.SetParameter(0, 1);
  f.SetParameter(1, cexp);
  f.SetParameter(2, 0);
  f.SetParameter(3, gMean);
  f.SetParameter(4, gSigma);
  f.SetParameter(0, 10*NBG/f.Integral(170,800));
  f.SetParameter(2, NSig*10);

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    int const NThisBin = (int) f.Integral( h->GetBinLowEdge(i), h->GetBinLowEdge(i) + h->GetBinWidth(i)) / 10.;
    printf("PEExpoWithBump ibin/N %3i %5i\n", i, NThisBin);
    for (int ie = 0; ie < NThisBin; ++ie) {
      h->Fill( h->GetBinCenter(i) );
    }
    //printf("Error %3i %5.0f %8.3f  %8.3f\n", i, h->GetBinContent(i), sqrt(h->GetBinContent(i)), h->GetBinError(i));
  }

  if (true) {
    TF1 fE("myexp", "[0]*TMath::Exp([1]*x)", 170, 800);
    fE.SetParameter(0, f.GetParameter(0));
    fE.SetParameter(1, f.GetParameter(1));
    TCanvas c;
    c.cd();
    h->Draw("ep");
    fE.Draw("same");
    c.SaveAs("TestPEWithSig.eps");
  }

  return h;
}




TH1D* GetHistForMjjj (float const Mjjj, TFile* File, TString const Suffix)
{
  // This function will return a histogram
  // Before, I had put all PEs in a file and ran over them.. now generate
  // on the go.. so this may or may not be used for more than data..

  char BUFF[200];
  sprintf(BUFF,  "Mjjj_45_20_130_" + Suffix);
  //if (Mjjj <= 250)        sprintf(BUFF,  "Mjjj_45_20_130_" + Suffix);
  //else                    sprintf(BUFF,  "Mjjj_45_20_170_" + Suffix);


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
  return -5.56179969055969892e-03 - 4.01623842089755741e-06 * Mjjj + 2.49580149009780901e-07 * Mjjj * Mjjj;
  //if (Mjjj <= 250) {
  //  return -5.56179969055969892e-03 - 4.01623842089755741e-06 * Mjjj + 2.49580149009780901e-07 * Mjjj * Mjjj;
  //} else {
  //  return 1.93282628706596682e-03 - 3.61649088607910918e-05 * Mjjj + 1.84254776396677824e-07 * Mjjj * Mjjj;
  //}

}



std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  // This will eventually be a parametrization

  if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  if (Mjjj < 350) return std::make_pair<float, float>(10 + 10 * (Mjjj - 250.) / 100., 15 + 10 * (Mjjj - 250.) / 100.);
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







float DoFit (RooWorkspace& ws, RooStats::ModelConfig& modelConfig, TString const label, int const method, int const Section)
{
  // Used for the limits
  float upper = -1;

  // The different methods
  if(method==0) {
    RooStats::BayesianCalculator bc(*ws.data("DataToFit"), modelConfig);
    bc.SetConfidenceLevel(0.95);
    bc.SetLeftSideTailFraction(0.0);
    RooStats::SimpleInterval* bInt = bc.GetInterval();
    if (Section == -1) {
      TCanvas* canvas = new TCanvas("bcPlost","Posterior Distribution",500,500);
      RooPlot* plot = bc.GetPosteriorPlot();
      plot->Draw();
      canvas->SaveAs(TString("BCpost")+label+".eps");
    }
    upper=bInt->UpperLimit();

  } else if(method==1) {
    //RooStats::ProposalHelper ph;
    //ph.SetVariables((RooArgSet&)fit->floatParsFinal());
    //ph.SetCovMatrix(fit->covarianceMatrix());
    //ph.SetUpdateProposalParameters(true);
    //ph.SetCacheSize(100);
    //RooStats::ProposalFunction* pdfProp = ph.GetProposalFunction();

    RooStats::MCMCCalculator mc(*ws.data("DataToFit"), modelConfig);
    mc.SetNumBins(50);
    mc.SetConfidenceLevel(0.95);
    mc.SetLeftSideTailFraction(0.0);
    mc.SetNumIters(10000000);
    mc.SetNumBurnInSteps(500);
    //mc.SetProposalFunction(*pdfProp);
    RooStats::MCMCInterval* interval = (RooStats::MCMCInterval*)mc.GetInterval();

    // Upper limit
    upper=interval->UpperLimit(*ws.var("xs"));

    // draw posterior plot for data
    if (Section == -1) {
      TCanvas * c1 = new TCanvas("MCMCpost", "MCMCpost", 500, 500);
      RooStats::MCMCIntervalPlot mcPlot(*interval);
      mcPlot.Draw();
      c1->SaveAs(TString("mcmcPost")+label+".eps");
    }

  } else if(method==2) {
    RooStats::ProfileLikelihoodCalculator plc(*ws.data("DataToFit"), modelConfig);
    plc.SetConfidenceLevel(0.95);
    RooStats::LikelihoodInterval* plInt=plc.GetInterval();
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    plInt->LowerLimit( *ws.var("xs") ); // get ugly print out of the way. Fix.
    RooMsgService::instance().setGlobalKillBelow(msglevel);

    // draw posterior plot
    if (Section == -1) {
      TCanvas * c1 = new TCanvas("PLikePost", "PLikePost", 500, 500);
      RooStats::LikelihoodIntervalPlot* plotInt=new RooStats::LikelihoodIntervalPlot(plInt);
      plotInt->SetTitle("Profile Likelihood Ratio and Posterior for S");
      plotInt->Draw();
      c1->SaveAs(TString("PLikePost")+label+".eps");
    }

    upper = plInt->UpperLimit(*ws.var("xs"));
  }
  printf("UPPER LIMIT: %12.3f\n", upper);

  return upper;
}



TH1F* GetMeAHist (TH1F* h)
{
  TH1F* g = new TH1F("MyNewHist", "MyNewHist", 63, 170, 800);

  for (int i = 17; i != 80; ++i) {
    for (int j = 0; j != h->GetBinContent(i); ++j) {
      g->Fill(10*i - 0.5);
    }
  }

  if (false) {
    TCanvas* c = new TCanvas();
    c->cd();
    h->Draw();
    g->SetLineColor(6);
    g->Draw("same");
    c->SaveAs("NewHist.eps");
  }
  return g;
}






float RunMultiJetsRooStats (TString const InFileName, float const SignalMass, int const method, int const statLevel, int const Section)
{
  // Declare some constants
  float const MININVMASS =    0;
  float const MAXINVMASS = 1500;

  float const LUMINOSITY = 35.1;
  float const LUMIERROR  = 0.11;

  //float const ACCERROR   = 0.13; // Old
  float const ACCERROR   = 0.18; // Include MC stat pileup

  float const MINXS      =    0;
  float const MAXXS      = 1000;

  // Just a label
  char label[100];
  sprintf(label, "Dean_%i_%i_%i", (int) SignalMass, method, statLevel);


  // Setup workspace
  RooWorkspace ws("ws");

  // Observable
  ws.factory("mjjj[0]");
  RooRealVar* mjjj = ws.var("mjjj");

  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1F* DataTH1F = (TH1F*) GetHistForMjjj(SignalMass, &InFile, "6jet");
  TH1* DataTH1 = (TH1*) GetMeAHist(DataTH1F);
  //TH1* DataTH1 = (TH1*) GetHistForMjjj(SignalMass, &InFile, "6jet");
  RooDataHist binnedData("data", "dataset with x", *ws.var("mjjj"), DataTH1);
  RooDataHist fitData("DataToFit", "dataset with x", *ws.var("mjjj"), DataTH1);
  ws.import(binnedData);
  TF1* Func = DataTH1F->GetFunction("total");
  printf("Fit parameters: %12.3f %12.3f\n",
      Func->GetParameter(0),
      Func->GetParameter(1)
      );
  float const XMIN = 170;//Func->GetXmin();
  float const XMAX = Func->GetXmax();
  int   const BINSMJJJ   = (int) (XMAX - XMIN) / 10;
  ws.var("mjjj")->setBins(BINSMJJJ);
  printf("Fit Min/Max: %12.2f %12.2f\n", XMIN, XMAX);
  //float const NBGFromFit = Func->Integral(XMIN, XMAX) / DataTH1->GetBinWidth(0);
  float const NBGFromFit = DataTH1->Integral();
  printf("Data vs fit: %12.1f %12.1f\n", DataTH1->Integral(), Func->Integral(XMIN, XMAX) / DataTH1->GetBinWidth(0));

  mjjj->setRange(XMIN, XMAX);

  // Other fit, Used for syst numbers
  TH1* SystTH1 = (TH1*) GetHistForMjjj(SignalMass, &InFile, "4jetscaled");
  TF1* SystFunc = SystTH1->GetFunction("total");
  printf("Syst Numbers: %12.3f %12.3ff\n",
      SystFunc->GetParameter(0) / Func->GetParameter(0),
      SystFunc->GetParameter(1) / Func->GetParameter(1));

  if (Section == -1) {
    TCanvas Can1;
    Can1.cd();
    DataTH1->Draw("ep");
    Func->Draw("samel");
    SystFunc->Draw("lsame");
    Can1.SaveAs(TString("DefaultFit_")+label+".eps");
  }


  // background prior
  ws.factory("RooExponential::background(mjjj, cexp[-1,1])");
  ws.factory("RooGaussian::nbkg_prior(nbkg[0,60000], nbkgM0[0], nbkgS0[0])");
  ws.factory("RooGaussian::cexp_prior(cexp, cexpM0[0], cexpS0[0])");
  ws.var("nbkg")->setVal(NBGFromFit);
  ws.var("nbkg")->setRange(NBGFromFit * (1 - 3*0.03), NBGFromFit * (1 + 3*0.03));
  //ws.var("nbkg")->setRange(NBGFromFit * (1 - 3*Func->GetParError(0)/Func->GetParameter(0)), NBGFromFit * (1 + 3*Func->GetParError(0)/Func->GetParameter(0)));
  ws.var("nbkgM0")->setVal(NBGFromFit);
  //ws.var("nbkgS0")->setVal(NBGFromFit * Func->GetParError(0)/Func->GetParameter(0));
  ws.var("nbkgS0")->setVal(NBGFromFit * 0.03);
  ws.var("cexp")->setVal(Func->GetParameter(1));
  //ws.var("cexp")->setRange(Func->GetParameter(1)-3*Func->GetParError(1), Func->GetParameter(1)+3*Func->GetParError(1));
  ws.var("cexp")->setRange(Func->GetParameter(1)*(1+3*0.05), Func->GetParameter(1)*(1-3*0.05));
  ws.var("cexpM0")->setVal(Func->GetParameter(1));
  ws.var("cexpS0")->setVal(TMath::Abs(Func->GetParameter(1)*0.05));
  //ws.var("cexpS0")->setVal(Func->GetParError(1));


  // Lumi prior
  ws.factory("RooGaussian::lumi_prior(lumi[0], lumiM0[0], lumiS0[0])");
  ws.var("lumi")->setVal(LUMINOSITY);
  ws.var("lumi")->setRange(LUMINOSITY * (1. - 3. * LUMIERROR), LUMINOSITY * (1. + 3. * LUMIERROR));
  ws.var("lumiM0")->setVal(LUMINOSITY);
  ws.var("lumiS0")->setVal(LUMINOSITY * LUMIERROR);

  // cross section prior and set the allowed range for the cross section
  ws.factory("xs[0]");
  ws.var("xs")->setRange(MINXS, MAXXS);
  ws.factory("RooUniform::xs_prior(xs)");

  // define number of signal as xs*lumi*acceptance
  ws.factory("prod::nsig(xs, lumi, acceptance[0,1])");


  // Acceptance prior
  float const Acc = GetAcceptanceForMjjj(SignalMass);
  ws.factory("RooGaussian::acceptance_prior(acceptance, acceptanceM0[0], acceptanceS0[0])");
  ws.var("acceptance")->setVal(Acc);
  ws.var("acceptance")->setRange(Acc * (1. - 3. * ACCERROR), Acc * (1. + 3. * ACCERROR));
  ws.var("acceptanceM0")->setVal(Acc);
  ws.var("acceptanceS0")->setVal(Acc * ACCERROR);

  // set the signal mass
  ws.factory("sigMass[0]");
  ws.var("sigMass")->setRange(MININVMASS, MAXINVMASS);
  ws.var("sigMass")->setVal(SignalMass);

  // define signal gaussian
  ws.factory("RooGaussian::signal(mjjj, sigMass, sigWidth[0,50])");
  ws.var("sigWidth")->setVal( GetGausWidthRange(SignalMass).first/2. + GetGausWidthRange(SignalMass).second/2.);
  ws.var("sigWidth")->setRange( GetGausWidthRange(SignalMass).first, GetGausWidthRange(SignalMass).second);


  // Define models (sig+bg, and bg only)
  ws.factory("SUM::model_noprior(nsig*signal, nbkg*background)");
  ws.factory("SUM::background_noprior(nbkg*background)");

  // Observables and parameters or interest
  ws.defineSet("observables","mjjj");
  ws.defineSet("POI","xs");

  // Priors and nuis params
  // use this one
  ws.factory("PROD::prior0(xs_prior)");
  ws.defineSet("nuisSet0","sigWidth");
  ws.factory("PROD::prior5b(xs_prior,nbkg_prior,cexp_prior,lumi_prior,acceptance_prior)");
  ws.defineSet("nuisSet5b","nbkg,cexp,lumi,acceptance,sigWidth");
  

  // let's not fix the cross section..just to make sure
  ws.var("xs")->setConstant(false);

  //  Fix everything and turn off individually based on what we need later..
  ws.var("lumi")->setConstant(true);
  ws.var("lumiM0")->setConstant(true);
  ws.var("lumiS0")->setConstant(true);
  ws.var("acceptance")->setConstant(true);
  ws.var("acceptanceM0")->setConstant(true);
  ws.var("acceptanceS0")->setConstant(true);
  ws.var("cexp")->setConstant(true);
  ws.var("cexpM0")->setConstant(true);
  ws.var("cexpS0")->setConstant(true);
  ws.var("nbkg")->setConstant(true);
  ws.var("nbkgM0")->setConstant(true);
  ws.var("nbkgS0")->setConstant(true);
  ws.var("sigMass")->setConstant(true);

  // setup model config
  RooStats::ModelConfig modelConfig("modelConfig");
  modelConfig.SetWorkspace(ws);

  // Which prior and nuis are we using
  if (statLevel == 0) {
    ws.factory("PROD::model(model_noprior,prior0)");
    ws.factory("PROD::background_model(background_noprior,prior0)");
    ws.var("sigWidth")->setConstant(false);
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet0"));
  } else if (statLevel == 6) {
    ws.factory("PROD::model(model_noprior,prior5b)");
    ws.factory("PROD::background_model(background_noprior,prior5b)");
    ws.var("sigWidth")->setConstant(false);
    ws.var("nbkg")->setConstant(false);
    ws.var("cexp")->setConstant(false);
    ws.var("lumi")->setConstant(false);
    ws.var("acceptance")->setConstant(false);
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet5b"));
  } else {
    std::cerr << "HELP!" << std::endl;
    exit(1);
  }

  modelConfig.SetPdf(*ws.pdf("model"));
  modelConfig.SetParametersOfInterest(*ws.set("POI"));

  // import the modelconfig
  ws.import(modelConfig);



  // Just some print statements
  ws.Print();
  modelConfig.Print();
  ws.var("nbkg")->Print();
  ws.var("nbkgM0")->Print();
  ws.var("nbkgS0")->Print();
  ws.var("cexp")->Print();
  ws.var("cexpM0")->Print();
  ws.var("cexpS0")->Print();
  ws.var("lumi")->Print();
  ws.var("lumiM0")->Print();
  ws.var("lumiS0")->Print();
  ws.var("acceptance")->Print();
  ws.var("xs")->Print();
  ws.var("sigMass")->Print();
  ws.var("sigWidth")->Print();





  float upper = -1;
  if (Section == -2) {
    // Setup some stuff for toy MC experiments TOY, everyone likes TOYS, but not to be toyed with.
    // I can get you a toy, believe me.  There are ways, Dude.


    float const MyNBG  = Func->Integral(170, 800) / 10.;
    float const MyCEXP = Func->GetParameter(1);
    printf("ADASD: %12E %12E\n", DataTH1->Integral(), MyCEXP);
    ws.var("nbkg")->setVal(MyNBG);

    TH1F* hPEWithBump = GetPEExpoWithBump(DataTH1->Integral(), MyCEXP, 379.995/10., 380, 22.681);
    RooDataHist hPE("DataToFit", "dataset with x", *ws.var("mjjj"), hPEWithBump);

    //RooDataHist* hPE = PE->binnedClone("DataToFit");
    printf("Data vs PE: %12.1f %12.1f\n", fitData.sum(false), hPE.sum(false));

    // Just for debug
    if (false) {
      TCanvas CanPE;
      CanPE.cd();
      RooPlot* PEplot = ws.var("mjjj")->frame();
      hPE.plotOn(PEplot);
      PEplot->Draw();
      CanPE.SaveAs(TString("PE_")+label+".eps");
    }

    ws.import(hPE);
    if (false) {
      ws.var("nbkg")->Print();
      TCanvas Can;
      Can.cd();
      RooPlot* datafit = ws.var("mjjj")->frame();
      ws.data("DataToFit")->plotOn(datafit);
      ws.pdf("model")->fitTo(*ws.data("DataToFit"));
      ws.pdf("model")->plotOn(datafit);
      datafit->Draw();
      Can.SaveAs(TString("Fit_")+label+".eps");
      ws.var("nbkg")->Print();
    }
    upper = DoFit(ws, modelConfig, label, method, Section);

    delete hPEWithBump;
  } else if (Section == -1) {
    ws.import(fitData);
    ws.var("nbkg")->Print();

    TCanvas Can;
    Can.cd();
    RooPlot* datafit = ws.var("mjjj")->frame();
    ws.data("DataToFit")->plotOn(datafit);
    ws.pdf("model")->fitTo(*ws.data("DataToFit"), RooFit::Range(XMIN,XMAX), RooFit::Extended(kTRUE));
    ws.pdf("model")->plotOn(datafit);
    datafit->Draw();
    Can.SaveAs(TString("Fit_")+label+".eps");
    ws.var("nbkg")->Print();

    upper = DoFit(ws, modelConfig, label, method, Section);
  } else {
    // Setup some stuff for toy MC experiments TOY, everyone likes TOYS, but not to be toyed with.
    // I can get you a toy, believe me.  There are ways, Dude.


    float const rand_NBGThisPE = NBGFromFit + ws.var("nbkgS0")->getVal() * RooRandom::randomGenerator()->Gaus(0,1);

    // Grab a toy.  If you can't get a toy there is something wrong
    int const NBGThisPE = RooRandom::randomGenerator()->Poisson(rand_NBGThisPE);
    ws.var("nbkg")->setVal(NBGThisPE);

    float rand_cexp = 999;
    while (rand_cexp < ws.var("cexp")->getMin() || rand_cexp > ws.var("cexp")->getMax()) {
      rand_cexp = ws.var("cexp")->getVal() + ws.var("cexpS0")->getVal() * RooRandom::randomGenerator()->Gaus(0,1);
      std::cout << "rand_cexp:  " << rand_cexp << std::endl;
    }

    int NTest = RooRandom::randomGenerator()->Poisson(
      1.*(NBGFromFit + ws.var("nbkgS0")->getVal() * RooRandom::randomGenerator()->Gaus(0,1))
        );
    TH1F* hPEF = GetPEExpo(NTest, rand_cexp, label);
    RooDataHist hPE("DataToFit", "dataset with x", *ws.var("mjjj"), hPEF);

    //RooDataHist* hPE = PE->binnedClone("DataToFit");
    printf("Data vs PE: %12.1f %12.1f\n", fitData.sum(false), hPE.sum(false));

    // Just for debug
    if (false) {
      TCanvas CanPE;
      CanPE.cd();
      RooPlot* PEplot = ws.var("mjjj")->frame();
      hPE.plotOn(PEplot);
      PEplot->Draw();
      CanPE.SaveAs(TString("PE_")+label+".eps");
    }

    ws.import(hPE);
    if (false) {
      ws.var("nbkg")->Print();
      TCanvas Can;
      Can.cd();
      RooPlot* datafit = ws.var("mjjj")->frame();
      ws.data("DataToFit")->plotOn(datafit);
      ws.pdf("model")->fitTo(*ws.data("DataToFit"));
      ws.pdf("model")->plotOn(datafit);
      datafit->Draw();
      Can.SaveAs(TString("Fit_")+label+".eps");
      ws.var("nbkg")->Print();
    }
    upper = DoFit(ws, modelConfig, label, method, Section);

    delete hPEF;

  }



  return upper;
}


int main (int argc, char* argv[])
{
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [Section] (if -1, mass)" << std::endl;
    return 1;
  }

  float const StepSize  =  10;

  int const Section = atoi(argv[1]);
  float const BeginMass = argc == 3 ?  200 + atof(argv[2])*StepSize : 200;
  float const EndMass   = argc == 3 ?  200 + atof(argv[2])*StepSize : 500;
  std::cout << BeginMass << "  " << EndMass << std::endl;

  TString const InFileName = "/home/dhidas/Data35pb/ExpFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  //TString const InFileName = "/uscms/home/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";
  //TString const InFileName = "/Users/dhidas/Data35pb/ExpoFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root";

  //float const BeginMass = 200;
  //float const EndMass   = 500;

  int const Method      =  1;
  int const Systematics =  6;
  int const NPerSection =  2;

  // Set the roostats random seed based on secton number..fine..
  RooRandom::randomGenerator()->SetSeed(771723*(Section+2));
  gRandom->SetSeed(791723*(Section+2));


  // Setup output file
  char OutName[150];
  if (Section == -1 || Section == -2) {
    if (argc == 2) {
      sprintf(OutName, "Limits_Data.dat");
    } else if (argc == 3) {
      sprintf(OutName, "Limits_Data_%i.dat", (int) BeginMass);
    }
  } else {
    sprintf(OutName, "PELimits_%i.dat", Section);
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



  // Print header to outfile


  if (Section == -1 || Section == -2) {
    std::vector< std::pair<float, float> > MassLimitVec;
    for (float Mass = BeginMass; Mass <= EndMass; Mass += StepSize) {
      MassLimitVec.push_back( std::make_pair<float, float>(Mass, RunMultiJetsRooStats(InFileName, Mass, Method, Systematics, Section)) );
      fprintf(Out, "%10E ", MassLimitVec.back().second);
    }
    fprintf(Out, "\n");
    fflush(Out);
    MakeGraph(MassLimitVec, "95\% C.L.", "M_{jjj}", "95\% C.L. #sigma (pb)", "Limits_Data.eps");
  } else {
    for (int ipe = Section * NPerSection; ipe < (Section+1) * NPerSection; ++ipe) {
      std::vector< std::pair<float, float> > MassLimitVec;
      for (float Mass = BeginMass; Mass <= EndMass; Mass += StepSize) {
        MassLimitVec.push_back( std::make_pair<float, float>(Mass, RunMultiJetsRooStats(InFileName, Mass, Method, Systematics, Section)) );
        fprintf(Out, "%10E ", MassLimitVec.back().second);
      }
      fprintf(Out, "\n");
      fflush(Out);
    }
  }

  return 0;
}
