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
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SequentialProposal.h"
#include "RooRandom.h"


// Min and max for observable
float const MJJJMIN =  230;
float const MJJJMAX = 1500;

// Luminosity
//float const LUMINOSITY = 2177.0;
float const LUMINOSITY = 1127.0;
float const LUMIERROR  =   0.045;


TH1F* GetPE4Param (int const N, float const p1, float const p2, float const p3, TString const label = "blah")
{
  TH1F* h = new TH1F("PE", "PE", 65, 230, 1520);

  TF1 f("myfuncPE", "[0]*(((1. - x/7000.)^p1)/((x/7000.)^(p2+ p3*log(x/7000.))))", 230, 1520);
  f.SetParameter(0, 1);
  f.SetParameter(1, p1);
  f.SetParameter(2, p2);
  f.SetParameter(2, p3);
  f.SetParameter(0, 1.0/f.Integral(230,1520));
  std::cout << "PE NORM" << "  " << f.Integral(230,1520) << std::endl;

  for (int i = 0; i < N; ++i) {
    h->Fill(f.GetRandom());
  }

  TF1 ff("mypseudo1", "[0]*(((1. - x/7000.)^p1)/((x/7000.)^(p2+ p3*log(x/7000.))))", 230, 1520);
  ff.SetParameter(0, 1);
  ff.SetParameter(1, p1);
  ff.SetParameter(2, p2);
  ff.SetParameter(2, p3);
  ff.SetParameter(0, 10*N/ff.Integral(230,1520));
  std::cout << N << "  " << ff.Integral(230,1520) << std::endl;

  if (false) {
    TCanvas c;
    c.cd();
    h->Draw();
    ff.Draw("same");
    c.SaveAs(TString("MyPE_")+label+".pdf");
  }

  return h;
}



TH1F* GetPE3Param (int const N, float const p1, float const p2, TString const label = "blah")
{
  TH1F* h = new TH1F("PE", "PE", 65, 230, 1520);

  TF1 f("myfuncPE", "[0]*(((1. - x/7000.)^[1])/((x/7000.)^[2]))", 230, 1520);
  f.SetParameter(0, 1);
  f.SetParameter(1, p1);
  f.SetParameter(2, p2);
  f.SetParameter(0, 1.0/f.Integral(230,1520));
  std::cout << "PE NORM" << "  " << f.Integral(230,1520) << std::endl;

  for (int i = 0; i < N; ++i) {
    h->Fill(f.GetRandom());
  }

  TF1 ff("mypseudo1", "[0]*(((1. - x/7000.)^[1])/((x/7000.)^[2]))", 230, 1520);
  ff.SetParameter(0, 1);
  ff.SetParameter(1, p1);
  ff.SetParameter(2, p2);
  ff.SetParameter(0, 10*N/ff.Integral(230,1520));
  std::cout << N << "  " << ff.Integral(230,1520) << std::endl;

  if (false) {
    TCanvas c;
    c.cd();
    h->Draw();
    ff.Draw("same");
    c.SaveAs(TString("MyPE_")+label+".pdf");
  }

  return h;
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

  return (0.140834 + 0.000392742*m - -7.49412e-07*m*m + 4.50564e-10*m*m*m);
}









float RunMultJetsMCMC(TString const InFileName, int const Section, float const ThisMass)
{
  // Get the mass for this section
  float const MassStepSize = 20;
  float const SignalMass = ThisMass;
  std::cout << "SignalMass = " << SignalMass << std::endl;

  float const MINXS      =      0;
  float const MAXXS      =    100;

  float const MINPOI     =      0;
  float const MAXPOI     =   
    SignalMass <  300 ?  35   :
    SignalMass <  400 ?  13   :
    SignalMass <  500 ?  10   :
    SignalMass <  600 ?   6   :
    SignalMass <  650 ?   3   :
    SignalMass <  700 ?   2   :
    SignalMass <  780 ?   1.8 :
    SignalMass <  840 ?   1.3 :
    SignalMass < 860 ?   1.0 :
    SignalMass < 950 ?   0.7 :
    SignalMass < 1030 ?   0.5 :
    SignalMass < 1110 ?   0.4 :
    SignalMass < 1270 ?   0.3 :
    SignalMass < 1370 ?   0.18 :
    SignalMass < 1500 ?   0.16 :
    0.16;


  // Setup output root file based on mass
  TFile OutRootFile(TString::Format("MultiJetsCLs_%i.root", (int) SignalMass), "recreate");
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
  TF1* FitFunction = (TF1*) DataHist->GetFunction("g4");
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
  //ws.var("p2")->setRange(FitFunction->GetParameter(2) - 5 * FitFunction->GetParError(2), FitFunction->GetParameter(2) + 5 * FitFunction->GetParError(2));
  ws.var("p2")->setRange(0, FitFunction->GetParameter(2) + 5 * FitFunction->GetParError(2));
  ws.var("p2")->setVal(FitFunction->GetParameter(2));
  //ws.var("p2")->setVal(-FitFunction->GetParameter(2));
  ws.var("p2")->setConstant(true);
  ws.var("p2")->Print();
  ws.factory("p3[0]");
  //ws.var("p3")->setRange(FitFunction->GetParameter(3) - 5 * FitFunction->GetParError(3), FitFunction->GetParameter(3) + 5 * FitFunction->GetParError(3));
  ws.var("p3")->setRange(0, FitFunction->GetParameter(3) + 5 * FitFunction->GetParError(3));
  ws.var("p3")->setVal(FitFunction->GetParameter(3));
  ws.var("p3")->setConstant(true);
  ws.var("p3")->Print();


  // define priors (constraints) for these params
  ws.factory("RooLognormal::p1_prior(p1, p1M0[0], p1S0[1])");
  ws.var("p1S0")->setVal(1.0 + FitFunction->GetParError(1) / FitFunction->GetParameter(1));
  ws.var("p1S0")->setConstant(true);
  ws.var("p1M0")->setVal(FitFunction->GetParameter(1));
  ws.var("p1M0")->setConstant(true);
  ws.factory("RooLognormal::p2_prior(p2, p2M0[0], p2S0[1])");
  ws.var("p2S0")->setVal(1.0 + FitFunction->GetParError(2) / FitFunction->GetParameter(2));
  ws.var("p2S0")->setConstant(true);
  ws.var("p2M0")->setVal(FitFunction->GetParameter(2));
  ws.var("p2M0")->setConstant(true);
  ws.factory("RooLognormal::p3_prior(p3, p3M0[0], p3S0[1])");
  ws.var("p3S0")->setVal(1.0 + FitFunction->GetParError(3) / FitFunction->GetParameter(3));
  ws.var("p3S0")->setConstant(true);
  ws.var("p3M0")->setVal(FitFunction->GetParameter(3));
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
  ws.var("acceptance")->setRange(0, Acceptance * (1. + 3. * AcceptanceError));
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
  ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJJMIN, MJJJMAX), RooFit::Extended(kTRUE));
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
  RooFitResult* FitResult = ws.pdf("sbmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJJMIN, MJJJMAX), RooFit::Extended(kTRUE));
  if (FitResult == 0x0) {
    std::cout << "FitResult = " << FitResult << std::endl;
    exit(0);
  }
  TCanvas CanFit("Fit", "Fit");
  CanFit.cd();
  RooPlot* ThisDataFit = ws.var("mjjj")->frame();
  ws.data("Data")->plotOn(ThisDataFit);
  ws.pdf("bgmodel_noprior")->plotOn(ThisDataFit);
  ThisDataFit->SetTitle(TString::Format("M_{jjj} = %i", (int) SignalMass));
  ThisDataFit->Draw();
  CanFit.SetLogy(1);
  OutRootFile.cd();
  CanFit.Write();
  CanFit.SaveAs(TString::Format("Fit_%i.eps", (int) SignalMass));

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


  float UpperLimit = -999;
  RooStats::MCMCCalculator* mc;
  RooStats::MCMCInterval* interval;
  if (Section <= 0) {
    // implement cov matrix for data here
    RooStats::ProposalHelper ph;
    ph.SetVariables((RooArgSet&) FitResult->floatParsFinal());
    //ph.SetCovMatrix(FitResult->covarianceMatrix());
    exit(0);

    
    mc = new RooStats::MCMCCalculator(*ws.data("Data"), ModelConfigSB);
    mc->SetNumBins(50);
    mc->SetConfidenceLevel(0.95);
    mc->SetLeftSideTailFraction(0.0);
    RooStats::SequentialProposal sp(0.1);
    mc->SetProposalFunction(sp);
    mc->SetNumIters(5000000);
    mc->SetNumBurnInSteps(5);





    interval = (RooStats::MCMCInterval*) mc->GetInterval();
  } else {

    float const rand_NBGThisPE = NBackground + (ws.var("nbkgS0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);

    // Grab a toy.  If you can't get a toy there is something wrong
    int const NBGThisPE = RooRandom::randomGenerator()->Poisson(rand_NBGThisPE);
    ws.var("nbkg")->setVal(NBGThisPE);

    //float rand_p0 = 999;
    //while (rand_p0 < ws.var("p0")->getMin() || rand_p0 > ws.var("p0")->getMax()) {
    //  rand_p0 = ws.var("p0")->getVal() + (ws.var("p0S0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);
    //  std::cout << "rand_p0:  " << rand_p0 << std::endl;
    //}

    float rand_p1 = 999;
    while (rand_p1 < ws.var("p1")->getMin() || rand_p1 > ws.var("p1")->getMax()) {
      rand_p1 = ws.var("p1")->getVal() + (ws.var("p1S0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);
      std::cout << "rand_p1:  " << rand_p1 << std::endl;
    }

    float rand_p2 = 999;
    while (rand_p2 < ws.var("p2")->getMin() || rand_p2 > ws.var("p2")->getMax()) {
      rand_p2 = ws.var("p2")->getVal() + (ws.var("p2S0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);
      std::cout << "rand_p2:  " << rand_p2 << std::endl;
    }

    float rand_p3 = 999;
    while (rand_p3 < ws.var("p3")->getMin() || rand_p3 > ws.var("p3")->getMax()) {
      rand_p3 = ws.var("p3")->getVal() + (ws.var("p3S0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);
      std::cout << "rand_p3:  " << rand_p3 << std::endl;
    }

    int NTest = RooRandom::randomGenerator()->Poisson(
      1.*(NBackground + (ws.var("nbkgS0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1))
        );
    TH1F* hPEF = GetPE4Param(NTest, rand_p1, rand_p2, rand_p3, "blah");
    RooDataHist hPE("DataToFit", "dataset with x", *ws.var("mjjj"), hPEF);
    ws.import(hPE);

    mc = new RooStats::MCMCCalculator(*ws.data("DataToFit"), ModelConfigSB);
    mc->SetNumBins(50);
    mc->SetConfidenceLevel(0.95);
    mc->SetLeftSideTailFraction(0.0);
    mc->SetNumIters(500);
    mc->SetNumBurnInSteps(5);
    interval = (RooStats::MCMCInterval*) mc->GetInterval();
  }
  UpperLimit = interval->UpperLimit(*ws.var("xs"));

  delete mc;
  delete interval;

  // Upper limit
  printf("UPPER LIMIT: %12.3f\n", UpperLimit);


  // Close root file
  OutRootFile.Close();

  return UpperLimit;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Section = atoi(argv[2]);

  int const NPerSection = Section <= 0 ? 1 : 5;

  RooRandom::randomGenerator()->SetSeed(339723*(Section+2));
  gRandom->SetSeed(276823*(Section+2));

  // Setup output file
  TString OutName;
  if (Section <= 0) {
    OutName = "Limits_Data.dat";
  } else {
    OutName = TString::Format("PELimits_%i.dat", Section);
  }
  FILE* Out = fopen(OutName.Data(), "w");

  for (float Mass = 280; Mass <= 1480; Mass += 100) {
    fprintf(Out, "%10.3f ", Mass);
  }
  fprintf(Out, "\n");

  for (int ipe = Section * NPerSection; ipe < (Section + 1) * NPerSection; ++ipe) {
    for (float Mass = 280; Mass <= 1480; Mass += 100) {
      fprintf(Out, "%10E ", RunMultJetsMCMC(InFileName, Section, Mass));
    }
    fprintf(Out, "\n");
  }
  fclose(Out);

  return 0;
}
