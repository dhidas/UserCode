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
#include "TMatrixDSym.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalFunction.h"
#include "RooStats/SequentialProposal.h"
#include "TMatrixTSym.h"
#include "RooFitResult.h"

// Min and max for observable
//command out by yik
//float const MJJMIN =  280;
//float const MJJMAX = 1200;

float  MJJMIN =  280;
float  MJJMAX = 1240;


// Luminosity
/*  //2011A
float const LUMINOSITY = 2164.0+35.5;
float const LUMIERROR  =   0.045;
*/

float const LUMINOSITY = 4974.0;
float const LUMIERROR  =   0.022;

float GetAcceptanceForMjj (float const Mjj)
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

  //return -0.01027 + 4.38331e-05*Mjj + 4.43791e-08*Mjj*Mjj - 3.69426e-11*Mjj*Mjj*Mjj;
  //return -0.0173967 + 8.54121e-05 * Mjj + -2.44194e-08 * Mjj * Mjj;
  return 1.0;

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

  //  return (0.140834 + 0.000392742*m - -7.49412e-07*m*m + 4.50564e-10*m*m*m);
  return 0.10;

}




int RunFourJetsCLs (TString const InFileName, int const Section, int const SeedOffset)
{
  // Get the mass for this section
//  float const SignalMass = 300 + (40 * Section);
  float const SignalMass = 300 + (40 * Section);

  std::cout << "SignalMass = " << SignalMass << std::endl;

 // Set the roostats random seed based on secton number..fine..
  RooRandom::randomGenerator()->SetSeed(7723*(Section+SeedOffset));
  std::cout << "random seed= = " << (7723*(Section+SeedOffset)) << std::endl;

//added by yik
MJJMIN=280;
if(SignalMass>1100) MJJMAX=1400;
//added end

  float const MINXS      =      0;
  float const MAXXS      =    5;

//  float const MINPOI     =      0;
  float MINPOI     =
    SignalMass <  250 ?  0.001   :
    SignalMass <  300 ?  0.001   :
    SignalMass <  330 ?  0.001   :
    SignalMass <  370 ?  0.001   :
    SignalMass <  410 ?  0.001   :
    SignalMass <  450 ?  0.001   :
    SignalMass <  490 ?  0.0001   :    
    SignalMass <  530 ?  0.0001   :    
    SignalMass <  570 ?   0.001  :
    SignalMass <  610 ?   0.001   :
    SignalMass <  650 ?   0.001   :
    SignalMass <  690 ?   0.001   :   
    SignalMass <  730 ?   0.001 :
    SignalMass <  770 ?   0.001 :
    SignalMass <  810 ?   0.001 :
    SignalMass < 850 ?   0.001 :
    SignalMass < 890 ?   0.001 :
    SignalMass < 930 ?   0.001 :
    SignalMass < 970 ?   0.001 :   
    SignalMass < 1010 ?   0.001 :
    SignalMass < 1050 ?   0.001 :
    SignalMass < 1090 ?   0.001 :
    SignalMass < 1130 ?   0.001 :
    SignalMass < 1170 ?   0.001 :
    SignalMass < 1210 ?   0.001 :
    SignalMass < 1500 ?   0.001 :
    0.001;
  MINPOI = 0;

  

  float const MAXPOI     =   
    SignalMass <  250 ?  0.08   :
    SignalMass <  300 ?  0.08   :
    SignalMass <  330 ?  0.08   :
    SignalMass <  370 ?  0.08   :
    SignalMass <  410 ?  0.08   :
    SignalMass <  450 ?  0.08   :
    SignalMass <  490 ?  0.06   :    //to get observed limit accurately
    SignalMass <  530 ?  0.06   :    //to get observed limit accurately
    SignalMass <  570 ?   0.05  :
    SignalMass <  610 ?   0.05   : 
    SignalMass <  650 ?   0.04   :
    SignalMass <  690 ?   0.04   :
    SignalMass <  730 ?   0.03 :
    SignalMass <  770 ?   0.03 :
    SignalMass <  810 ?   0.025 :
    SignalMass < 850 ?   0.02 :
    SignalMass < 890 ?   0.02 :
    SignalMass < 930 ?   0.02 :
    SignalMass < 970 ?   0.01 :   //to get observed limit accurately
    SignalMass < 1010 ?   0.015 :
    SignalMass < 1050 ?   0.012 :
    SignalMass < 1090 ?   0.01 :
    SignalMass < 1130 ?   0.01 :
    SignalMass < 1170 ?   0.01 :
    SignalMass < 1210 ?   0.01 :
    SignalMass < 1500 ?   0.008 :
    SignalMass < 1550 ?   0.008 :
    SignalMass < 1600 ?   0.008 :
    SignalMass < 1650 ?   0.006 :
    SignalMass < 1700 ?   0.006 :
    SignalMass < 1750 ?   0.006 :
    SignalMass < 1800 ?   0.006 :
    SignalMass < 1850 ?   0.006 :
    SignalMass < 1900 ?   0.006 :
    SignalMass < 1950 ?   0.006 :
    0.006;

  
  // Open output file for limits
  std::ofstream OutFile(TString::Format("Limits_CLs_%i.dat", (int) SignalMass).Data());
  if (!OutFile.is_open()) {
    std::cerr << "ERROR: cannot open output file Limits_CLs_*.dat for mass: " << SignalMass << std::endl;
    throw;
  }

  // Setup output root file based on mass
  TFile OutRootFile(TString::Format("MultiJetsCLs_%i_%i.root", (int) SeedOffset, (int) SignalMass), "recreate");
  if (!OutRootFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output root file for mass: " << SignalMass << std::endl;
    throw;
  }
  OutRootFile.cd();



  // Start a workspace
  RooWorkspace ws("ws");

  // Obseravable
  ws.factory("mjj[0]");
  ws.var("mjj")->setRange(MJJMIN, MJJMAX);

  // Get the data hist and import it to the workspace
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    exit(1);
  }
  TH1* DataHist = (TH1*) InFile.Get("mypfbestPairdijetmassAvg");
  if (DataHist == 0x0) {
    std::cerr << "ERROR: cannot get data histogram" << std::endl;
    throw;
  }
  std::cout << "Got data hist with number of entries: " << DataHist->GetEntries() << std::endl;
  RooDataHist Data("Data", "dataset with x", *ws.var("mjj"), DataHist);
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
  float const NData = DataHist->Integral( DataHist->FindBin(MJJMIN), DataHist->FindBin(MJJMAX) );
  float const NBackground = FitFunction->Integral(MJJMIN, MJJMAX) / DataHist->GetBinWidth(0);
  printf("Data vs fit: %12.1f %12.1f\n", NData, NBackground);

  // Define background and parameters
  ws.factory("p1[0]");
  ws.var("p1")->setRange(FitFunction->GetParameter(1) - 5 * FitFunction->GetParError(1), FitFunction->GetParameter(1) + 5 * FitFunction->GetParError(1));
  ws.var("p1")->setVal(FitFunction->GetParameter(1));
  ws.var("p1")->setConstant(true);
  ws.var("p1")->Print();
  ws.factory("p2[0]");
  ws.var("p2")->setRange(FitFunction->GetParameter(2) - 5 * FitFunction->GetParError(2), FitFunction->GetParameter(2) + 5 * FitFunction->GetParError(2));
  ws.var("p2")->setVal(-FitFunction->GetParameter(2));
  ws.var("p2")->setConstant(true);
  ws.var("p2")->Print();
  ws.factory("p3[0]");
  ws.var("p3")->setRange(FitFunction->GetParameter(3) - 5 * FitFunction->GetParError(3), FitFunction->GetParameter(3) + 5 * FitFunction->GetParError(3));
  ws.var("p3")->setVal(-FitFunction->GetParameter(3));
  ws.var("p3")->setConstant(true);
  ws.var("p3")->Print();


  // Background function
  ws.factory("RooGenericPdf::background('(((1. - mjj/7000.)^p1)/((mjj/7000.)^(p2+ p3*log(mjj/7000.))))', {mjj, p1, p2, p3})");


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
  //  ws.factory("prod::nsig(xs, lumi, acceptance[0,1])");   //?? by yik, set a new low limit? changed by yik
  ws.factory("prod::nsig(xs, lumi, acceptance[0.9,1.1])"); 
  ws.factory("nsigX[0,500]");
  ws.defineSet("POI","xs");

  float const Acceptance = GetAcceptanceForMjj(SignalMass);
  float const AcceptanceError = GetAcceptanceError(SignalMass);
  printf("Acceptance for mass %E is %E +/- %E\n", SignalMass, Acceptance, Acceptance * AcceptanceError);
  ws.factory("RooLognormal::acceptance_prior(acceptance, acceptanceM0[0], acceptanceS0[0])");
  ws.var("acceptanceS0")->setVal(1.0 + AcceptanceError);
  ws.var("acceptanceS0")->setConstant(true);
  //  ws.var("acceptance")->setRange(0, Acceptance * (1. + 3. * AcceptanceError));  //changed by yik
  ws.var("acceptance")->setRange(Acceptance * (1. - 3. * AcceptanceError), Acceptance * (1. + 3. * AcceptanceError)); 
  ws.var("acceptance")->setVal(Acceptance);
  ws.var("acceptance")->setConstant(true);
  ws.var("acceptanceM0")->setVal(Acceptance);
  ws.var("acceptanceM0")->setConstant(true);

  // Signal shape 
  //changed by yik to full shape
  /*
  ws.factory("RooGaussian::signal(mjj, sigMean[0], sigWidth[0,500])");
  ws.var("sigWidth")->setVal(SignalMass * 0.065);
  ws.var("sigWidth")->setRange( (SignalMass * 0.065) * (1.0 - 0.10), (SignalMass * 0.065) * (1.0 + 0.10));
  ws.var("sigWidth")->setConstant(true);
  ws.var("sigMean")->setVal(SignalMass);
  ws.var("sigMean")->setConstant(true);  
  */
  //get the numeric value for the pamaraterization, PFAK5
  float mygaus1mean,mygaus1sigma,mygaus2mean,mygaus2sigma, mygaus1frac;
  //gaus1
  mygaus1mean=SignalMass;
  mygaus1sigma=4.55992+0.0393607*SignalMass-5.71813e-06*SignalMass*SignalMass;
  //gaus2
 if(SignalMass<600.0)
   {mygaus2mean=2792.46-16.8088*SignalMass+0.0401528*SignalMass*SignalMass-3.93535e-05*SignalMass*SignalMass*SignalMass+1.39073e-08*SignalMass*SignalMass*SignalMass*SignalMass;}
 else
   {mygaus2mean=-43.8974+1.02325*SignalMass-0.000294366*SignalMass*SignalMass;}
 // cout<<"SignalMass="<<SignalMass<<", mygaus2mean="<<mygaus2mean<<endl;
  if(SignalMass<600.0)
    { mygaus2sigma=-3835.75+32.4216*SignalMass-0.0942709*SignalMass*SignalMass+0.00011812*SignalMass*SignalMass*SignalMass-5.41591e-08*SignalMass*SignalMass*SignalMass*SignalMass; }
  else
    { mygaus2sigma=158.513-0.0153738*SignalMass+6.72179e-05*SignalMass*SignalMass; }
  //  cout<<"SignalMass="<<SignalMass<<", mygaus1sigma="<<mygaus1sigma<<", mygaus1sigma="<<mygaus1sigma<<", mygaus2sigma="<<mygaus2sigma<<<<endl;
  //gaus1fraction
  if(SignalMass<600.0) 
    { mygaus1frac=1.12214-0.000786882*SignalMass-1.80846e-07*SignalMass*SignalMass; }
  else 
    { mygaus1frac=0.988864-0.00083739*SignalMass+2.80641e-07*SignalMass*SignalMass; }
  //  cout<<"SignalMass="<<SignalMass<<", mygaus1frac="<<mygaus1frac<<endl;
  std::cout<<"SignalMass="<<SignalMass<<", mygaus1mean="<<mygaus1mean<<", mygaus1sigma="<<mygaus1sigma<<", mygaus2mean="<<mygaus2mean<<", mygaus2sigma="<<mygaus2sigma<<", mygaus1frac="<<mygaus1frac<<std::endl;
  ws.factory("Gaussian::myg1(mjj, sigMean1[0], sigWidth1[0,60])");
  ws.factory("Gaussian::myg2(mjj, sigMean2[0], sigWidth2[70,1500])");
  //  ws.Print();
  ws.factory("SUM::signal(g1frac[0.5,0,1]*myg1,myg2)") ;
  ws.var("sigMean1")->setVal(mygaus1mean);
  ws.var("sigMean1")->setConstant(true);
  ws.var("sigWidth1")->setVal(mygaus1sigma);
  ws.var("sigWidth1")->setConstant(true);
  ws.var("sigMean2")->setVal( mygaus2mean );
  ws.var("sigMean2")->setConstant(true);
  ws.var("sigWidth2")->setVal( mygaus2sigma );
  ws.var("sigWidth2")->setConstant(true);
  ws.var("g1frac")->setVal( mygaus1frac );
  ws.var("g1frac")->setConstant(true);
  
  // prior on the signal width
  ws.factory("RooUniform::sigWidth1_prior(sigWidth1)");

  std::cout << 1.0 + FitFunction->GetParError(1) / FitFunction->GetParameter(1) << std::endl;
  std::cout << 1.0 + FitFunction->GetParError(2) / FitFunction->GetParameter(2) << std::endl;
  std::cout << 1.0 + FitFunction->GetParError(3) / FitFunction->GetParameter(3) << std::endl;

  // nbkg prior
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

  /*
  //added by yik, try to use multigaussian function
  //TMatrixDSym covFit=ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Extended(kTRUE),RooFit::Save())->covarianceMatrix();
  //myresult->covarianceMatrix().Print();
  //TMatrixDSym mycovFit=myresult->covarianceMatrix();
  //mycovFit=ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Extended(kTRUE),RooFit::Save() )->covarianceMatrix();
  //  ws.factory("RooMultiVarGaussian::mybkg_prior(RooArgList(nbkg,p1,p2,p3), myresult)");
  //  ws.factory("RooMultiVarGaussian::mybkg_prior( {nbkg,p1,p2,p3}, myresult)");

  //  RooFitResult * myresult=ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Extended(kTRUE),RooFit::Save(true) );
  RooFitResult * myresult=ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Save(true) );
  RooStats::SequentialProposal sp(0.1);
  RooStats::ProposalHelper ph;
  ph.SetVariables( myresult->floatParsFinal() );
  TMatrixDSym covFit=myresult->covarianceMatrix();
  //ph.SetCovMatrix(  myresult->covarianceMatrix() );
  ph.SetUpdateProposalParameters(kTRUE);
  ph.SetCacheSize(100);
  RooStats::ProposalFunction* mybkg_prior=ph.GetProposalFunction();
  //added by yik end
  */

  
  // Pick which constraints and nuisance params you want to use
  switch (6) {    //was 6 for default,  7 to check gaus1 as a nuisance parameter
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
      ws.factory("PROD::constraints(nbkg_prior,p1_prior,p2_prior,p3_prior,lumi_prior,acceptance_prior,sigWidth1_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3,lumi,acceptance,sigWidth1");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      ws.var("lumi")->setConstant(false);
      ws.var("acceptance")->setConstant(false);
      ws.var("sigWidth1")->setRange( mygaus1sigma * (1.0 - 0.10), mygaus1sigma * (1.0 + 0.10));
      ws.var("sigWidth1")->setConstant(false);
      break;
    case 8:
      ws.factory("PROD::constraints(mybkg_prior,lumi_prior,acceptance_prior,sigWidth1_prior)");
      ws.defineSet("nuisance", "nbkg,p1,p2,p3,lumi,acceptance,sigWidth1");
      ws.var("nbkg")->setConstant(false);
      ws.var("p1")->setConstant(false);
      ws.var("p2")->setConstant(false);
      ws.var("p3")->setConstant(false);
      ws.var("lumi")->setConstant(false);
      ws.var("acceptance")->setConstant(false);
      ws.var("sigWidth1")->setRange( mygaus1sigma * (1.0 - 0.10), mygaus1sigma * (1.0 + 0.10));
      ws.var("sigWidth1")->setConstant(false);
      break;
    default:
      std::cerr << "pick something I have please" << std::endl;
      throw;
  }
  ws.var("p1")->Print();
  ws.var("p2")->Print();
  ws.var("p3")->Print();
  ws.var("nbkg")->Print();
  ws.var("lumi")->Print();
  ws.var("acceptance")->Print();

//the old wa to set CLs limit, then we will subtract the possible signal to background shape by Tomasso's opinion
/*
  // Build modelconfig for bgmodel and set current workspace
  RooStats::ModelConfig ModelConfigBG("ModelConfigBG");
  ModelConfigBG.SetWorkspace(ws);

  // Setup this model
  ModelConfigBG.SetPdf(*ws.pdf("bgmodel_noprior"));
  ModelConfigBG.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigBG.SetObservables(*ws.var("mjj"));
  ModelConfigBG.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigBG.SetNuisanceParameters(*ws.set("nuisance"));

  // Set the observable to zero and add POI snapsnot, import this modelconfig to the workspace
  ws.var("xs")->setVal(0);
  ws.pdf("bgmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Extended(kTRUE));
  RooArgSet POIAndNuisBG("POIAndNuisBG");
  POIAndNuisBG.add(*ModelConfigBG.GetParametersOfInterest());
  ModelConfigBG.SetSnapshot(POIAndNuisBG);
  ModelConfigBG.SetGlobalObservables( RooArgSet() );
  ws.import(ModelConfigBG);
*/

  // Build modelconfig for sbmodel and set current workspace
  RooStats::ModelConfig ModelConfigSB("ModelConfigSB");
  ModelConfigSB.SetWorkspace(ws);

  // Setup this model
  ModelConfigSB.SetPdf(*ws.pdf("sbmodel_noprior"));
  ModelConfigSB.SetParametersOfInterest(*ws.set("POI"));
  ModelConfigSB.SetObservables(*ws.var("mjj"));
  ModelConfigSB.SetPriorPdf(*ws.pdf("constraints"));
  ModelConfigSB.SetNuisanceParameters(*ws.set("nuisance"));

  // Fit model and add POI snapsnot, import this modelconfig to the workspace
  ws.pdf("sbmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJMIN, MJJMAX), RooFit::Extended(kTRUE));
  TCanvas CanFit("Fit", "Fit");
  CanFit.cd();
  RooPlot* ThisDataFit = ws.var("mjj")->frame();
  ws.data("Data")->plotOn(ThisDataFit);
  ws.pdf("bgmodel_noprior")->plotOn(ThisDataFit);
  ThisDataFit->SetTitle(TString::Format("M_{jj} = %i", (int) SignalMass));
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

//the new way to set CLs limit by Tommaso's method--subtract possible signal from data for background shape.
      RooStats::ModelConfig * ModelConfigBG = ( RooStats::ModelConfig*) ModelConfigSB.Clone();
      ModelConfigBG->SetName("ModelConfigBG");     
      RooRealVar * var = dynamic_cast<RooRealVar*>(ModelConfigBG->GetParametersOfInterest()->first());
      var->setVal(0);
      ModelConfigBG->SetSnapshot( RooArgSet(*var)  );


  // Parameters of the CLs method we'll call
  int   const calculatorType    = 2;
  int   const testStatType      = 3;
  bool  const useCls            = true;
  int   const npoints           = 200;   //20
  float const poimin            = MINPOI;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
  float const poimax            = MAXPOI; //1;//60 / (LUMINOSITY * GetAcceptanceForMjj(SignalMass));
  int   const ntoys             = 1000; //use 100, 500
  bool  const useNumberCounting = false;
  const char* nuisPriorName     = "";


  // Run the actual CLs
  RooStats::HypoTestInvTool HTIT;
  HTIT.SetParameter("GenerateBinned", true);
  RooStats::HypoTestInverterResult* MyResult = HTIT.RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, useCls, npoints, poimin, poimax, ntoys, useNumberCounting, nuisPriorName);
  /*
  cout<<"Hello 2"<<endl;
  // Print the limits
  printf(" expected limit (-2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-2));
  printf(" expected limit (-1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-1));
  printf(" expected limit (median) %12.3E\n", MyResult->GetExpectedUpperLimit(0) );
  printf(" expected limit (+1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(1) );
  printf(" expected limit (+2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(2) );
  printf(" observed limit          %12.3E +/- %12.3E\n", MyResult->UpperLimit(), MyResult->UpperLimitEstimatedError());


  // Write results and close file
  OutFile << SignalMass << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-2) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 0) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 2) << std::endl;
  OutFile << MyResult->UpperLimit() << std::endl;
  OutFile.close();
 
  // Number of entries in result
  const int NEntries = MyResult->ArraySize();

  // Just some names
  const char* TypeName = "Hybrid";
  const char* ResultName = MyResult->GetName();
  TString PlotTitle = TString::Format("%s CL Scan for workspace %s", TypeName, ResultName);

  // Grab the result plot
  RooStats::HypoTestInverterPlot *Plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", PlotTitle, MyResult);
  TCanvas CanCLb("CLb_2CL", "CLb_2CL");
  CanCLb.cd();
  Plot->SetTitle(TString::Format("Hybrid CL Scan for M_{jj} = %i", (int) SignalMass));
  Plot->Draw("CLb 2CL");  // plot all and Clb
  CanCLb.SaveAs(TString::Format("CLb2L_%i.eps", (int) SignalMass));
  


  OutRootFile.cd();
  CanCLb.Write();

  // Draw the sampling distributions
  TCanvas CanHTI("HTI_Result", "HTI_Result");
  CanHTI.Divide(3, (int) TMath::Ceil(NEntries/3));
  for (int i = 0; i < NEntries; ++i) {
    CanHTI.cd(i + 1);
    RooStats::SamplingDistPlot * SamplingPlot = Plot->MakeTestStatPlot(i);
    //SamplingPlot->SetLogYaxis(true);
    delete SamplingPlot;
  }
  CanHTI.SaveAs(TString::Format("HTI_Result_%i.eps", (int) SignalMass));
  OutRootFile.cd();
  CanHTI.Write();
  */



//comment out by yik
/*
  // Print the limits
  printf(" expected limit (-2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-2));
  printf(" expected limit (-1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-1));
  printf(" expected limit (median) %12.3E\n", MyResult->GetExpectedUpperLimit(0) );
  printf(" expected limit (+1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(1) );
  printf(" expected limit (+2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(2) );
  printf(" observed limit          %12.3E +/- %12.3E\n", MyResult->UpperLimit(), MyResult->UpperLimitEstimatedError()); 


  // Write results and close file
  OutFile << SignalMass << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-2) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 0) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 2) << std::endl;
  OutFile << MyResult->UpperLimit() << std::endl;
  OutFile.close();
*/


  // Close root file

  OutRootFile.cd();
  MyResult->Write();
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

  RunFourJetsCLs(InFileName, Section, SeedOffset);

  return 0;
}




/*
int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Section = atoi(argv[2]);

  RunFourJetsCLs(InFileName, Section, SeedOffset);

  return 0;
}
*/
