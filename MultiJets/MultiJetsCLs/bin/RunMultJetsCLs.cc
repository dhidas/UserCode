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


// Min and max for observable
float const MJJJMIN =  230;
float const MJJJMAX = 1500;

// Luminosity
float const LUMINOSITY = 2177.0;
float const LUMIERROR  =   0.045;




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









int RunMultJetsCLs (TString const InFileName, int const Section)
{
  // Get the mass for this section
  float const SignalMass = 250 + (20 * Section);
  std::cout << "SignalMass = " << SignalMass << std::endl;

  float const MINXS      =      0;
  float const MAXXS      =    100;

  float const MINPOI     =      0;
  float const MAXPOI     =   
    SignalMass <  300 ?  30   :
    SignalMass <  400 ?  10   :
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


  // Open output file for limits
  std::ofstream OutFile(TString::Format("Limits_CLs_%i.dat", (int) SignalMass).Data());
  if (!OutFile.is_open()) {
    std::cerr << "ERROR: cannot open output file Limits_CLs_*.dat for mass: " << SignalMass << std::endl;
    throw;
  }

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
  ws.var("p2")->setRange(FitFunction->GetParameter(2) - 5 * FitFunction->GetParError(2), FitFunction->GetParameter(2) + 5 * FitFunction->GetParError(2));
  ws.var("p2")->setVal(-FitFunction->GetParameter(2));
  ws.var("p2")->setConstant(true);
  ws.var("p2")->Print();
  ws.factory("p3[0]");
  ws.var("p3")->setRange(FitFunction->GetParameter(3) - 5 * FitFunction->GetParError(3), FitFunction->GetParameter(3) + 5 * FitFunction->GetParError(3));
  ws.var("p3")->setVal(-FitFunction->GetParameter(3));
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
  switch (4) {
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
  ws.pdf("sbmodel_noprior")->fitTo(*ws.data("Data"), RooFit::Range(MJJJMIN, MJJJMAX), RooFit::Extended(kTRUE));
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

  // Parameters of the CLs method we'll call
  int   const calculatorType    = 1;
  int   const testStatType      = 3;
  bool  const useCls            = true;
  int   const npoints           = 20;
  float const poimin            = MINPOI;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
  float const poimax            = MAXPOI; //1;//60 / (LUMINOSITY * GetAcceptanceForMjjj(SignalMass));
  int   const ntoys             = 600;
  bool  const useNumberCounting = false;
  const char* nuisPriorName     = "";

  // Run the actual CLs
  RooStats::HypoTestInverterResult* MyResult = RunInverter(&ws, "ModelConfigSB", "ModelConfigBG", "Data", calculatorType, testStatType, npoints, poimin, poimax, ntoys, useCls, useNumberCounting, nuisPriorName);

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
  Plot->SetTitle(TString::Format("Hybrid CL Scan for M_{jjj} = %i", (int) SignalMass));
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
    SamplingPlot->SetLogYaxis(true);
    delete SamplingPlot;
  }
  CanHTI.SaveAs(TString::Format("HTI_Result_%i.eps", (int) SignalMass));
  OutRootFile.cd();
  CanHTI.Write();

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

  // Close root file
  OutRootFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section]" << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];
  int const Section = atoi(argv[2]);

  RunMultJetsCLs(InFileName, Section);

  return 0;
}
