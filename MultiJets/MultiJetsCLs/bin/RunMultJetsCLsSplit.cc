////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Sat Oct  1 06:02:06 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>

#include "StandardHypoTestInvDemo.h"

#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooRandom.h"


// Min and max for observable
float const MJJJMIN =  260;
float const MJJJMAX = 1620;

// Luminosity
float const LUMINOSITY = 4980.0;
float const LUMIERROR  =   0.022;
//float const LUMINOSITY = 4632.0;
//float const LUMIERROR  =   0.045;




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


  // Update by dan on 29/2/2012
  return -0.0172876 + 8.94367e-05 * Mjjj - 4.05152e-08 * Mjjj*Mjjj;

  // One used for most 4.6/fb studies
  //return -2.12900022490592165e-02 + 1.08977199027775987e-04*Mjjj + -6.75844597940117972e-08*Mjjj*Mjjj + 1.03590838339000778e-11*Mjjj*Mjjj*Mjjj;

  // old 2.X
  //return -0.01027 + 4.38331e-05*Mjjj + 4.43791e-08*Mjjj*Mjjj - 3.69426e-11*Mjjj*Mjjj*Mjjj;
  //return -0.0173967 + 8.54121e-05 * Mjjj + -2.44194e-08 * Mjjj * Mjjj;
}


float GetAcceptanceError (float const m)
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

  // acceptance error test remove later on
  //if (m >= 800) {
  //  return 1.1653400;
  //}

  //return (0.140834 + 0.000392742*m - -7.49412e-07*m*m + 4.50564e-10*m*m*m);
  return 0.199;
  //return 0.25;
}









int RunMultJetsCLsSplit (TString const InFileName, int const Section, int const SeedOffset)
{
  // Get the mass for this section
  float const SignalMass = (int) MJJJMIN + (10 * Section);
  std::cout << "SignalMass = " << SignalMass << std::endl;

  // Set the roostats random seed based on secton number..fine..
  RooRandom::randomGenerator()->SetSeed(7 * (100 * Section + SeedOffset));


  // Min and max ax for roovar
  float const MINXS      =      0;
  float const MAXXS      =    200;

  // Min and max xs for CLs
  float const MINPOI     =      0;

  float const MAXPOI     =   
    SignalMass <=  260 ?   70   :
    SignalMass <=  280 ?   70   :

    SignalMass <=  300 ?   60   :
    SignalMass <=  320 ?   40   :
    SignalMass <=  340 ?   30   :
    SignalMass <=  360 ?   25   :
    SignalMass <=  380 ?   20   :

    SignalMass <=  400 ?   17   :
    SignalMass <=  420 ?   16   :
    SignalMass <=  440 ?   13   :
    SignalMass <=  460 ?   12   :
    SignalMass <=  480 ?   11   :

    SignalMass <=  500 ?    9   :
    SignalMass <=  520 ?    9   :
    SignalMass <=  540 ?    8   :
    SignalMass <=  560 ?    8   :
    SignalMass <=  580 ?    7   :

    SignalMass <=  600 ?    6   :
    SignalMass <=  620 ?    6   :
    SignalMass <=  640 ?    5   :
    SignalMass <=  660 ?    5   :
    SignalMass <=  680 ?    5   :

    SignalMass <=  700 ?    4   :
    SignalMass <=  720 ?    4   :
    SignalMass <=  740 ?    4   :
    SignalMass <=  760 ?    4   :
    SignalMass <=  780 ?    3   :

    SignalMass <=  800 ?    3   :
    SignalMass <=  820 ?    3   :
    SignalMass <=  840 ?    3   :
    SignalMass <=  860 ?    3   :
    SignalMass <=  880 ?    3   :

    SignalMass <=  900 ?    2.5 :
    SignalMass <=  920 ?    2   :
    SignalMass <=  940 ?    2   :
    SignalMass <=  960 ?    2   :
    SignalMass <=  980 ?    2   :

    SignalMass <= 1000 ?    2   :
    SignalMass <= 1020 ?    1.5 :
    SignalMass <= 1040 ?    1.5 :
    SignalMass <= 1060 ?    1.5 :
    SignalMass <= 1080 ?    1.5 :

    SignalMass <= 1100 ?    1.5 :
    SignalMass <= 1120 ?    1.0 :
    SignalMass <= 1140 ?    2.0 :
    SignalMass <= 1160 ?    1.0 :
    SignalMass <= 1180 ?    0.9 :

    SignalMass <= 1200 ?    0.8 :
    SignalMass <= 1220 ?    0.8 :
    SignalMass <= 1240 ?    0.7 :
    SignalMass <= 1260 ?    0.7 :
    SignalMass <= 1280 ?    0.6 :

    SignalMass <= 1300 ?    0.6 :
    SignalMass <= 1320 ?    0.6 :
    SignalMass <= 1340 ?    0.6 :
    SignalMass <= 1360 ?    0.5 :
    SignalMass <= 1380 ?    0.5 :

    SignalMass <= 1400 ?    0.5 :
    SignalMass <= 1420 ?    0.5 :
    SignalMass <= 1440 ?    0.4 :
    SignalMass <= 1460 ?    0.4 :
    SignalMass <= 1480 ?    0.4 :

    SignalMass <= 1500 ?    0.4 :
    SignalMass <= 1520 ?    0.4 :
    SignalMass <= 1540 ?    0.4 :
    SignalMass <= 1560 ?    0.4 :
    SignalMass <= 1580 ?    0.4 :

    SignalMass <= 1600 ?    0.35:
    SignalMass <= 1620 ?    0.35:
    SignalMass <= 1640 ?    0.35:
    SignalMass <= 1660 ?    0.35:
    SignalMass <= 1680 ?    0.35:
    0.35;
    /*
  float const MAXPOI     =   
    SignalMass <=  260 ?    5   :
    SignalMass <=  280 ?    5 :

    SignalMass <=  300 ?    5 :
    SignalMass <=  320 ?    5 :
    SignalMass <=  340 ?    1   :
    SignalMass <=  360 ?    1   :
    SignalMass <=  380 ?    1   :

    SignalMass <=  400 ?    0.5 :
    SignalMass <=  420 ?    0.5 :
    SignalMass <=  440 ?    0.5 :
    SignalMass <=  460 ?    0.4 :
    SignalMass <=  480 ?    0.4 :

    SignalMass <=  500 ?    0.3 :
    SignalMass <=  520 ?    0.3 :
    SignalMass <=  540 ?    0.3 :
    SignalMass <=  560 ?    0.3 :
    SignalMass <=  580 ?    0.3 :

    SignalMass <=  600 ?    0.25 :
    SignalMass <=  620 ?    0.25 :
    SignalMass <=  640 ?    0.25:
    SignalMass <=  660 ?    0.25 :
    SignalMass <=  680 ?    0.25 :

    SignalMass <=  700 ?    0.2 :
    SignalMass <=  720 ?    0.2 :
    SignalMass <=  740 ?    0.2 :
    SignalMass <=  760 ?    0.2 :
    SignalMass <=  780 ?    0.2 :

    SignalMass <=  800 ?    0.2 :
    SignalMass <=  820 ?    0.2 :
    SignalMass <=  840 ?    0.2 :
    SignalMass <=  860 ?    0.2 :
    SignalMass <=  880 ?    0.2 :

    SignalMass <=  900 ?    0.2 :
    SignalMass <=  920 ?    0.2 :
    SignalMass <=  940 ?    0.2 :
    SignalMass <=  960 ?    0.2 :
    SignalMass <=  980 ?    0.2 :

    SignalMass <= 1000 ?    0.2 :
    SignalMass <= 1020 ?    0.2 :
    SignalMass <= 1040 ?    0.2 :
    SignalMass <= 1060 ?    0.2 :
    SignalMass <= 1080 ?    0.2 :

    SignalMass <= 1100 ?    0.2 :
    SignalMass <= 1120 ?    0.2 :
    SignalMass <= 1140 ?    0.2 :
    SignalMass <= 1160 ?    0.2 :
    SignalMass <= 1180 ?    0.2 :

    SignalMass <= 1200 ?    0.2 :
    SignalMass <= 1220 ?    0.2 :
    SignalMass <= 1240 ?    0.2 :
    SignalMass <= 1260 ?    0.2 :
    SignalMass <= 1280 ?    0.2 :

    SignalMass <= 1300 ?    0.2 :
    SignalMass <= 1320 ?    0.2 :
    SignalMass <= 1340 ?    0.2 :
    SignalMass <= 1360 ?    0.2 :
    SignalMass <= 1380 ?    0.2 :

    SignalMass <= 1400 ?    0.2 :
    SignalMass <= 1420 ?    0.2 :
    SignalMass <= 1440 ?    0.2 :
    SignalMass <= 1460 ?    0.2 :
    SignalMass <= 1480 ?    0.2 :

    SignalMass <= 1500 ?    0.2 :
    SignalMass <= 1520 ?    0.2 :
    SignalMass <= 1540 ?    0.2 :
    SignalMass <= 1560 ?    0.2 :
    SignalMass <= 1580 ?    0.2 :

    SignalMass <= 1600 ?    0.2 :
    SignalMass <= 1620 ?    0.2 :
    SignalMass <= 1640 ?    0.2 :
    SignalMass <= 1660 ?    0.2 :
    SignalMass <= 1680 ?    0.2 :
    0.2;
    */


  // Setup output root file based on mass
  TFile OutRootFile(TString::Format("MultiJetsCLs_%i_%i.root", SeedOffset, (int) SignalMass), "recreate");
  if (!OutRootFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output root file for mass: " << SignalMass << std::endl;
    throw;
  }
  OutRootFile.cd();




  // Start a workspace
  RooWorkspace ws("ws");

  // Obseravable
  ws.factory("mjjj[0]");
  ws.var("mjjj")->setRange(MJJJMIN, MJJJMAX);

  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1* DataHist = (TH1*) InFile.Get("Mjjj_70_20_160_6jet");
  if (DataHist == 0x0) {
    std::cerr << "ERROR: cannot get data histogram" << std::endl;
    throw;
  }
  std::cout << "Got data hist with number of entries: " << DataHist->GetEntries() << std::endl;
  RooDataHist Data("Data", "dataset with x", *ws.var("mjjj"), DataHist);
  ws.import(Data);


  // Get the fit function from the hist
  //TF1* FitFunction = (TF1*) DataHist->GetFunction("g4");
  TF1* FitFunction = (TF1*) DataHist->GetFunction("Three-Jet Mass Spectrum");
  if (FitFunction == 0x0) {
    std::cout << "ERROR: cannot get fitted function." << std::endl;
    throw;
  }
  std::cout << "Found fitted function: " << std::endl;
  FitFunction->Print();

  // Compare number from fit and data
  float const NData = DataHist->Integral( DataHist->FindBin(MJJJMIN), DataHist->FindBin(MJJJMAX) );
  float const NBackground = FitFunction->Integral(MJJJMIN, MJJJMAX) / DataHist->GetBinWidth(0);
  printf("Data vs fit: %12.1f %12.1f\n", NData, NBackground);


  // Define background and parameters
  ws.factory("p1[0]");
  ws.var("p1")->setRange(FitFunction->GetParameter(1) - 5 * FitFunction->GetParError(1), FitFunction->GetParameter(1) + 5 * FitFunction->GetParError(1));
  ws.var("p1")->setVal(FitFunction->GetParameter(1));
  ws.var("p1")->setConstant(true);
  ws.var("p1")->Print();
  ws.factory("p2[0]");
  ws.var("p2")->setRange(FitFunction->GetParameter(2) - 5 * FitFunction->GetParError(2), FitFunction->GetParameter(2) + 5 * FitFunction->GetParError(2));
  ws.var("p2")->setVal(FitFunction->GetParameter(2));
  ws.var("p2")->setConstant(true);
  ws.var("p2")->Print();
  ws.factory("p3[0]");
  ws.var("p3")->setRange(FitFunction->GetParameter(3) - 5 * FitFunction->GetParError(3), FitFunction->GetParameter(3) + 5 * FitFunction->GetParError(3));
  ws.var("p3")->setVal(FitFunction->GetParameter(3));
  ws.var("p3")->setConstant(true);
  ws.var("p3")->Print();


  // define priors (constraints) for these params
  ws.factory("RooLognormal::p1_prior(p1, p1M0[0], p1S0[1])");
  ws.var("p1S0")->setVal(1.0 + FitFunction->GetParError(1) / FitFunction->GetParameter(1));
  ws.var("p1S0")->setConstant(true);
  ws.var("p1M0")->setVal(FitFunction->GetParameter(1));
  ws.var("p1M0")->setConstant(true);
  ws.factory("RooLognormal::p2_prior(-p2, p2M0[0], p2S0[1])");
  ws.var("p2S0")->setVal(1.0 - FitFunction->GetParError(2) / FitFunction->GetParameter(2));
  ws.var("p2S0")->setConstant(true);
  ws.var("p2M0")->setVal(-FitFunction->GetParameter(2));
  ws.var("p2M0")->setConstant(true);
  ws.factory("RooLognormal::p3_prior(-p3, p3M0[0], p3S0[1])");
  ws.var("p3S0")->setVal(1.0 - FitFunction->GetParError(3) / FitFunction->GetParameter(3));
  ws.var("p3S0")->setConstant(true);
  ws.var("p3M0")->setVal(-FitFunction->GetParameter(3));
  ws.var("p3M0")->setConstant(true);


  // Background function
  ws.factory("RooGenericPdf::background('(((1. - mjjj/7000.)^p1)/((mjjj/7000.)^(p2+ p3*log(mjjj/7000.))))', {mjjj, p1, p2, p3})");


  // Define lumi and lumi prior
  ws.factory("RooLognormal::lumi_prior(lumi[0], lumiM0[0], lumiS0[0])");
  ws.var("lumi")->setVal(LUMINOSITY);
  ws.var("lumi")->setRange(LUMINOSITY * (1. - 3. * LUMIERROR), LUMINOSITY * (1. + 3. * LUMIERROR));
  ws.var("lumiM0")->setVal(LUMINOSITY);
  ws.var("lumiS0")->setVal(1.0 + LUMIERROR);


  // cross section prior and set the allowed range for the cross section
  ws.factory("xs[0]");
  ws.var("xs")->setRange(MINXS, MAXXS);
  ws.var("xs")->setConstant(false);
  ws.factory("RooUniform::xs_prior(xs)");

  // define number of signal as xs*lumi*acceptance
  ws.factory("prod::nsig(xs, lumi, acceptance[0,1])");
  ws.factory("nsigX[0,500]");
  ws.defineSet("POI","xs");

  float const Acceptance = GetAcceptanceForMjjj(SignalMass);
  float const AcceptanceError = GetAcceptanceError(SignalMass);
  printf("Acceptance for mass %E is %E +/- %E\n", SignalMass, Acceptance, Acceptance * AcceptanceError);
  ws.factory("RooLognormal::acceptance_prior(acceptance, acceptanceM0[0], acceptanceS0[0])");
  ws.var("acceptanceS0")->setVal(1.0 + AcceptanceError);
  ws.var("acceptanceS0")->setConstant(true);
  ws.var("acceptance")->setRange(Acceptance * (1 - 3. * AcceptanceError), Acceptance * (1. + 3. * AcceptanceError));
  ws.var("acceptance")->setVal(Acceptance);
  ws.var("acceptance")->setConstant(true);
  ws.var("acceptanceM0")->setVal(Acceptance);
  ws.var("acceptanceM0")->setConstant(true);



  // Signal shape
  ws.factory("RooGaussian::signal(mjjj, sigMean[0], sigWidth[0,500])");
  ws.var("sigWidth")->setVal(SignalMass * 0.065);
  ws.var("sigWidth")->setRange( (SignalMass * 0.065) * (1.0 - 0.10), (SignalMass * 0.065) * (1.0 + 0.10));
  ws.var("sigWidth")->setConstant(true);
  ws.var("sigMean")->setVal(SignalMass);
  ws.var("sigMean")->setConstant(true);

  // prior on the signal width
  ws.factory("RooLognormal::sigWidth_prior(sigWidth, sigWidthM0[0], sigWidthS0[0])");
  ws.var("sigWidthS0")->setVal( 1 + 0.03 );  // 3% error
  ws.var("sigWidthS0")->setConstant(true);
  ws.var("sigWidthM0")->setVal(ws.var("sigWidth")->getVal());
  ws.var("sigWidthM0")->setConstant(true);

  //ws.factory("RooUniform::sigWidth_prior(sigWidth)");

  // nbkg prior
  ws.factory("RooLognormal::nbkg_prior(nbkg[0,60000], nbkgM0[0], nbkgS0[0])");
  ws.var("nbkgS0")->setVal(1. + 0.03);  // 3% error
  ws.var("nbkgS0")->setConstant(true);
  ws.var("nbkgM0")->setVal(NBackground);
  ws.var("nbkgM0")->setConstant(true);

  ws.var("nbkg")->setVal(NBackground);
  ws.var("nbkg")->setRange(NBackground * (1 - 3*0.03), NBackground * (1 + 3*0.03));
  ws.var("nbkg")->setConstant(true);

  // Define models (sig+bg, and bg only)
  ws.factory("SUM::sbmodel_noprior(nsig*signal, nbkg*background)");
  ws.factory("SUM::bgmodel_noprior(nbkg*background)");


  // Pick which constraints and nuisance params you want to use
  switch (7) {
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
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      break;
    case 4:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior,p3_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      break;
    case 5:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior,p3_prior,lumi_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3,lumi");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      ws.var("lumi")->setConstant(false);
      break;
    case 6:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior,p3_prior,lumi_prior,acceptance_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3,lumi,acceptance");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      ws.var("lumi")->setConstant(false);
      ws.var("acceptance")->setConstant(false);
      break;
    case 7:
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior,p3_prior,lumi_prior,acceptance_prior,sigWidth_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3,lumi,acceptance,sigWidth");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      ws.var("lumi")->setConstant(false);
      ws.var("acceptance")->setConstant(false);
      ws.var("sigWidth")->setConstant(false);
      break;
    default:
      std::cerr << "pick something I have please" << std::endl;
      throw;
  }


  // Build modelconfig for bgmodel and set current workspace
  RooStats::ModelConfig ModelConfigBG("ModelConfigBG");
  ModelConfigBG.SetWorkspace(ws);

  // Setup this model
  ModelConfigBG.SetPdf(*ws.pdf("bgmodel_noprior"));
  ModelConfigBG.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigBG.SetObservables(*ws.var("mjjj"));
  ModelConfigBG.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigBG.SetNuisanceParameters(*ws.set("nuisance"));

  // Set the observable to zero and add POI snapsnot, import this modelconfig to the workspace
  ws.var("xs")->setVal(0);
  ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Extended(kTRUE));
  TCanvas CanFitBG("FitBG", "FitBG");
  CanFitBG.cd();
  RooPlot* ThisDataFitBG = ws.var("mjjj")->frame();
  ws.data("Data")->plotOn(ThisDataFitBG);
  ws.pdf("bgmodel_noprior")->plotOn(ThisDataFitBG);
  ThisDataFitBG->SetTitle(TString::Format("Background Fit M_{jjj} = %i", (int) SignalMass));
  ThisDataFitBG->Draw();
  CanFitBG.SetLogy(1);
  OutRootFile.cd();
  CanFitBG.Write();
  CanFitBG.SaveAs(TString::Format("Fit_%i_BG.eps", (int) SignalMass));


  RooArgSet POIAndNuisBG("POIAndNuisBG");
  POIAndNuisBG.add(*ModelConfigBG.GetParametersOfInterest());
  ModelConfigBG.SetSnapshot(POIAndNuisBG);
  ModelConfigBG.SetGlobalObservables( RooArgSet() );
  ws.import(ModelConfigBG);

  // Build modelconfig for sbmodel and set current workspace
  RooStats::ModelConfig ModelConfigSB("ModelConfigSB");
  ModelConfigSB.SetWorkspace(ws);

  // Setup this model
  ModelConfigSB.SetPdf(*ws.pdf("sbmodel_noprior"));
  ModelConfigSB.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigSB.SetObservables(*ws.var("mjjj"));
  ModelConfigSB.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigSB.SetNuisanceParameters(*ws.set("nuisance"));

  // Fit model and add POI snapsnot, import this modelconfig to the workspace
  ws.pdf("sbmodel_noprior")->fitTo(*ws.data("Data"),  RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit", "minimize"));
  TCanvas CanFitSB("FitSB", "FitSB");
  CanFitSB.cd();
  RooPlot* ThisDataFit = ws.var("mjjj")->frame();
  ws.data("Data")->plotOn(ThisDataFit);
  ws.pdf("sbmodel_noprior")->plotOn(ThisDataFit);
  ThisDataFit->SetTitle(TString::Format("Signal+Background Fit M_{jjj} = %i", (int) SignalMass));
  ThisDataFit->Draw();
  CanFitSB.SetLogy(1);
  OutRootFile.cd();
  CanFitSB.Write();
  CanFitSB.SaveAs(TString::Format("Fit_%i_SB.eps", (int) SignalMass));

  RooArgSet POIAndNuisSB("POIAndNuisSB");
  POIAndNuisSB.add(*ModelConfigSB.GetParametersOfInterest());
  ModelConfigSB.SetSnapshot(POIAndNuisSB);
  ModelConfigSB.SetGlobalObservables( RooArgSet() );
  ws.import(ModelConfigSB);

  // For debugging, why not save the workspace
  ws.Print();
  OutRootFile.cd();
  ws.Write();

  ws.var("acceptance")->Print();
  ws.var("lumi")->Print();
  ws.var("xs")->Print();
  ws.var("nbkg")->Print();


  // Parameters of the CLs method we'll call
  int   const calculatorType    = 2;
  int   const testStatType      = 3;
  bool  const useCls            = true;
  int   const npoints           = 50;
  float const poimin            = MINPOI;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
  float const poimax            = MAXPOI; //1;//60 / (LUMINOSITY * GetAcceptanceForMjjj(SignalMass));
  int   const ntoys             = 500;
  bool  const useNumberCounting = false;
  const char* nuisPriorName     = "";

  // Run the actual CLs
  RooStats::HypoTestInvTool HTIT;
  HTIT.SetParameter("GenerateBinned", true);
  RooStats::HypoTestInverterResult* MyResult = HTIT.RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, useCls, npoints, poimin, poimax, ntoys, useNumberCounting, nuisPriorName);

  // Save the HypoTestInverterResult object to root file
  OutRootFile.cd();
  MyResult->Write();

  // Close root file
  OutRootFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section] [SeedOffset]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Section = atoi(argv[2]);
  int const SeedOffset = atoi(argv[3]);

  RunMultJetsCLsSplit(InFileName, Section, SeedOffset);

  return 0;
}
