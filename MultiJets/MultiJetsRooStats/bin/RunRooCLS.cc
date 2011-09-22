////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Jan 19 12:26:55 CET 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "StandardHypoTestInvDemo.h"

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
#include "RooStats/HypoTestResult.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HypoTestInverterResult.h"


#include "TStopwatch.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"





// Utility functions


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










TH1D* GetHistForMjjj (float const Mjjj, TFile* File, TString const Suffix)
{
  // This function will return a histogram
  // Before, I had put all PEs in a file and ran over them.. now generate
  // on the go.. so this may or may not be used for more than data..

  char BUFF[200];
  sprintf(BUFF,  "Mjjj_70_20_160_" + Suffix);
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



float GetACCERROR (float const m)
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

  return 0.140834 + 0.000392742*m - -7.49412e-07*m*m + 4.50564e-10*m*m*m;
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



std::pair<float, float> GetGausWidthRange (float const Mjjj)
{
  // Get the range for the gaussian width you want to use for a given Mjjj
  // This will eventually be a parametrization

  // Talk with dan, we decided to go with 0.7 * mass: 0.065 - 0.075

  //return std::make_pair<float, float>(Mjjj * 0.07, Mjjj * 0.07);
  return std::make_pair<float, float>(Mjjj * 0.065, Mjjj * 0.075);

  //if (Mjjj < 250) return std::make_pair<float, float>(10, 15);
  //if (Mjjj < 350) return std::make_pair<float, float>(10 + 10 * (Mjjj - 250.) / 100., 15 + 10 * (Mjjj - 250.) / 100.);
  //return std::make_pair<float, float>(20, 25);
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







std::vector<float> DoFit (RooWorkspace& ws, RooStats::ModelConfig& modelConfig, RooStats::ModelConfig& modelConfigBG, TString const label, int const method, int const Section, float const MAXXS)
{
  // Used for the limits
  std::vector<float> Limits;

  // The different methods
  if(method==0) {
  } else if(method==1) {
  } else if(method==2) {


    int   const calculatorType = 1;
    int   const testStatType = 3;
    bool  const useCls = true;
    int   const npoints = 4;
    float const poimin = 0;   // Set to bigger than max and npoints to zero for search (observed makes sense, expected do on own )
    float const poimax = MAXXS;
    int   const ntoys = 5;
    bool  const useNumberCounting = false;
    const char* nuisPriorName = "prior5bWthXS";

    // For debugging
    ws.SaveAs(TString("DeanWS_")+label+".root");

    RooStats::HypoTestInverterResult* r = RunInverter(&ws, "modelConfig", "modelConfigBG", "data", calculatorType, testStatType, npoints, poimin, poimax, ntoys, useCls, useNumberCounting, nuisPriorName);

    double upper = r->UpperLimit();
    double ulError = r->UpperLimitEstimatedError();
    const int nEntries = r->ArraySize();



    std::cout << "The computed upper limit is: " << upper << " +/- " << ulError << std::endl;
    const char *  typeName = "Hybrid";
    const char * resultName = r->GetName();
    TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName,resultName);
    RooStats::HypoTestInverterPlot *plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
    plot->Draw("CLb 2CL");  // plot all and Clb

    if (Section == -1) {
      TCanvas c2;
      c2.Divide( 2, TMath::Ceil(nEntries/2));
      for (int i=0; i<nEntries; i++) {
        c2.cd(i+1);
        RooStats::SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
      }
      c2.SaveAs(TString("HTI_Result_")+label+".pdf");
    }

    Limits.push_back(r->GetExpectedUpperLimit(-2));
    Limits.push_back(r->GetExpectedUpperLimit(-1));
    Limits.push_back(r->GetExpectedUpperLimit(0));
    Limits.push_back(r->GetExpectedUpperLimit(1));
    Limits.push_back(r->GetExpectedUpperLimit(2));
    Limits.push_back(upper);

    printf(" expected limit (-2 sig) %12.3E\n", Limits[0]);
    printf(" expected limit (-1 sig) %12.3E\n", Limits[1]);
    printf(" expected limit (median) %12.3E\n", Limits[2]);
    printf(" expected limit (+1 sig) %12.3E\n", Limits[3]);
    printf(" expected limit (+2 sig) %12.3E\n", Limits[4]);
    printf(" observed limit          %12.3E\n", Limits[5]);

  } else {
    std::cerr << "ERROR: what method do you want anyway!??" << std::endl;
    throw;
  }

  return Limits;
}



TH1F* GetMeAHist (TH1F* h)
{
  float const xmin =  230;
  float const xmax = 1520;
  float const binsize = h->GetBinWidth(0);
  int   const nbins = (xmax - xmin)/binsize;
  TH1F* g = new TH1F("MyNewHist", "MyNewHist", nbins, xmin, xmax);

  int const firstbin = h->FindBin(xmin);
  int const lastbin  = h->FindBin(xmax);
  printf("nbins %i  min %f  max %f  first %i  last %i\n", nbins, xmin, xmax, firstbin, lastbin);

  for (int i = firstbin; i != lastbin; ++i) {
    for (int j = 0; j != h->GetBinContent(i); ++j) {
      g->Fill(binsize*i - 0.5);
    }
  }

  if (false) {
    TCanvas* c = new TCanvas();
    c->cd();
    h->Draw();
    g->SetLineColor(6);
    g->Draw("same");
    c->SaveAs("NewHist.pdf");
    delete c;
  }
  return g;
}






std::vector<float> RunMultiJetsRooStats (TString const InFileName, float const SignalMass, int const method, int const statLevel, int const Section)
{
  // Declare some constants
  float const MININVMASS =    230;
  float const MAXINVMASS =   1520;

  float const LUMINOSITY = 2177.0;
  float const LUMIERROR  =   0.045;

  //float const ACCERROR   =  0.13; // Old
  float const ACCERROR   =   GetACCERROR(SignalMass); // Include MC stat pileup and JES

  float const MINXS      =      0;
  float const MAXXS      =   SignalMass <  350 ? 100   :
                             SignalMass <  500 ?  30   : 
                             SignalMass <  800 ?  20   :
                             SignalMass <  900 ?  15   :
                             SignalMass < 1200 ?  10   :
                             SignalMass < 1400 ?  10   :
                             10;

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
  RooDataHist binnedData("data", "dataset with x", *ws.var("mjjj"), DataTH1);
  RooDataHist fitData("DataToFit", "dataset with x", *ws.var("mjjj"), DataTH1);
  ws.import(binnedData);

  // Get the fitted function and print values
  TF1* Func = DataTH1F->GetFunction("g3");
  printf("Fit parameters: %12.3f %12.3f %12.3f\n",
      Func->GetParameter(0),
      Func->GetParameter(1),
      Func->GetParameter(2)
      );

  // Define range and bins
  float const XMIN = 230;//Func->GetXmin();
  float const XMAX = 1520;//Func->GetXmax();
  float const BINSIZE = DataTH1F->GetBinWidth(0);
  int   const BINSMJJJ   = (int) (XMAX - XMIN) / BINSIZE;
  ws.var("mjjj")->setBins(BINSMJJJ);
  printf("Fit Min/Max: %12.2f %12.2f\n", XMIN, XMAX);

  float const NBGFromFit = DataTH1->Integral( DataTH1->FindBin(XMIN), DataTH1->FindBin(XMAX)   );
  printf("Data vs fit: %12.1f %12.1f\n", NBGFromFit, Func->Integral(XMIN, XMAX) / DataTH1->GetBinWidth(0));

  mjjj->setRange(XMIN, XMAX);

  // Other fit, Used for syst numbers
  TH1* SystTH1 = (TH1*) GetHistForMjjj(SignalMass, &InFile, "6jet");
  TF1* SystFunc = SystTH1->GetFunction("g3");
  printf("Syst Numbers: %12.3f %12.3f %12.3f\n",
      SystFunc->GetParameter(0) / Func->GetParameter(0),
      SystFunc->GetParameter(1) / Func->GetParameter(1),
      SystFunc->GetParameter(2) / Func->GetParameter(2));

  if (Section == -1) {
    TCanvas Can1;
    Can1.cd();
    DataTH1->Draw("ep");
    Func->Draw("samel");
    //SystFunc->Draw("lsame");
    Can1.SaveAs(TString("DefaultFit_")+label+".pdf");
  }

  // Define variables for the background function and the function itself
  //ws.factory("p0[0,2]");
  ws.factory("p1[40,60]");
  ws.factory("p2[-2,0]");
  ws.factory("RooGenericPdf::background('( ((1.0 - mjjj/7000.0)^p1) / (mjjj/7000.)^p2)', {mjjj, p1, p2})");




  // nbkg prior
  ws.factory("RooLognormal::nbkg_prior(nbkg[0,60000], nbkgM0[0], nbkgS0[0])");
  ws.var("nbkgS0")->setVal(1. + 0.03);
  ws.var("nbkg")->setVal(NBGFromFit);
  ws.var("nbkg")->setRange(NBGFromFit * (1 - 3*0.03), NBGFromFit * (1 + 3*0.03));
  ws.var("nbkgM0")->setVal(NBGFromFit);


  // prior on fit param 0
  //ws.factory("RooLognormal::p0_prior(p0, p0M0[0], p0S0[1])");
  //ws.var("p0S0")->setVal(1.05);
  //ws.var("p0")->setRange(Func->GetParameter(0)*(1-3*0.05), Func->GetParameter(0)*(1+3*0.05));
  //ws.var("p0")->setVal(Func->GetParameter(0));
  //ws.var("p0M0")->setVal(Func->GetParameter(0));

  // prior on fit param 1
  ws.factory("RooLognormal::p1_prior(p1, p1M0[0], p1S0[1])");
  ws.var("p1S0")->setVal(1.05);
  ws.var("p1")->setRange(Func->GetParameter(1)*(1-3*0.05), Func->GetParameter(1)*(1+3*0.05));
  ws.var("p1")->setVal(Func->GetParameter(1));
  ws.var("p1M0")->setVal(Func->GetParameter(1));

  // prior on fit param 0
  ws.factory("RooLognormal::p2_prior(-1.0*p2, -1.0*p2M0[0], p2S0[1])");
  ws.var("p2S0")->setVal(1.05);
  ws.var("p2")->setRange(Func->GetParameter(2)*(1+3*0.05), Func->GetParameter(2)*(1-3*0.05));
  ws.var("p2")->setVal(Func->GetParameter(2));
  ws.var("p2M0")->setVal(Func->GetParameter(2));


  // Lumi prior
  ws.factory("RooLognormal::lumi_prior(lumi[0], lumiM0[0], lumiS0[0])");
  ws.var("lumi")->setVal(LUMINOSITY);
  ws.var("lumi")->setRange(LUMINOSITY * (1. - 3. * LUMIERROR), LUMINOSITY * (1. + 3. * LUMIERROR));
  ws.var("lumiM0")->setVal(LUMINOSITY);
  ws.var("lumiS0")->setVal(1.0 + LUMIERROR);

  // cross section prior and set the allowed range for the cross section
  ws.factory("xs[0]");
  ws.var("xs")->setRange(MINXS, MAXXS);
  ws.factory("RooUniform::xs_prior(xs)");

  // define number of signal as xs*lumi*acceptance
  ws.factory("prod::nsig(xs, lumi, acceptance[0,1])");


  // Acceptance prior
  float const Acc = GetAcceptanceForMjjj(SignalMass);
  printf("Acceptance for mass %.1f is %E\n", SignalMass, Acc);
  ws.factory("RooLognormal::acceptance_prior(acceptance, acceptanceM0[0], acceptanceS0[0])");
  ws.var("acceptanceS0")->setVal(1.0 + ACCERROR);
  ws.var("acceptance")->setVal(Acc);
  ws.var("acceptance")->setRange(0, Acc * (1. + 3. * ACCERROR));
  ws.var("acceptanceM0")->setVal(Acc);


  // set the signal mass
  ws.factory("sigMass[0]");
  ws.var("sigMass")->setRange(MININVMASS, MAXINVMASS);
  ws.var("sigMass")->setVal(SignalMass);

  // define signal gaussian
  ws.factory("RooGaussian::signal(mjjj, sigMass, sigWidth[0,50])");
  ws.var("sigWidth")->setVal( GetGausWidthRange(SignalMass).first/2. + GetGausWidthRange(SignalMass).second/2.);
  ws.var("sigWidth")->setRange( GetGausWidthRange(SignalMass).first, GetGausWidthRange(SignalMass).second);

  ws.factory("RooUniform::sigWidth_prior(sigWidth)");


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
  ws.factory("PROD::prior5b(nbkg_prior,p1_prior,p2_prior,lumi_prior,acceptance_prior,sigWidth_prior)");
  ws.factory("PROD::prior5bWthXS(xs_prior,nbkg_prior,p1_prior,p2_prior,lumi_prior,acceptance_prior,sigWidth_prior)");
  ws.defineSet("nuisSet5b","nbkg,p1,p2,lumi,acceptance,sigWidth");
  ws.defineSet("nuisSet5bBG","nbkg,p1,p2");
  

  // let's not fix the cross section..just to make sure
  ws.var("xs")->setConstant(false);

  //  Fix everything and turn off individually based on what we need later..
  ws.var("lumi")->setConstant(true);
  ws.var("lumiM0")->setConstant(true);
  ws.var("lumiS0")->setConstant(true);
  ws.var("acceptance")->setConstant(true);
  ws.var("acceptanceM0")->setConstant(true);
  ws.var("acceptanceS0")->setConstant(true);
  //ws.var("p0")->setConstant(true);
  //ws.var("p0M0")->setConstant(true);
  //ws.var("p0S0")->setConstant(true);
  ws.var("p1")->setConstant(true);
  ws.var("p1M0")->setConstant(true);
  ws.var("p1S0")->setConstant(true);
  ws.var("p2")->setConstant(true);
  ws.var("p2M0")->setConstant(true);
  ws.var("p2S0")->setConstant(true);
  ws.var("nbkg")->setConstant(true);
  ws.var("nbkgM0")->setConstant(true);
  ws.var("nbkgS0")->setConstant(true);
  ws.var("sigMass")->setConstant(true);

  // setup model config
  RooStats::ModelConfig modelConfig("modelConfig");
  RooStats::ModelConfig modelConfigBG("modelConfigBG");
  modelConfig.SetWorkspace(ws);
  modelConfigBG.SetWorkspace(ws);

  // Which prior and nuis are we using
  if (statLevel == 0) {
    ws.factory("PROD::model(model_noprior,prior0)");
    ws.factory("PROD::background_model(background_noprior,prior0)");
    ws.var("sigWidth")->setConstant(false);
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet0"));
    modelConfigBG.SetNuisanceParameters(*ws.set("nuisSet0"));
  } else if (statLevel == 6) {
    ws.factory("PROD::model(model_noprior,prior5b)");
    ws.factory("PROD::background_model(background_noprior,prior5b)");
    ws.var("sigWidth")->setConstant(false);
    ws.var("nbkg")->setConstant(false);
    //ws.var("p0")->setConstant(false);
    ws.var("p1")->setConstant(false);
    ws.var("p2")->setConstant(false);
    ws.var("lumi")->setConstant(false);
    ws.var("acceptance")->setConstant(false);
    modelConfig.SetNuisanceParameters(*ws.set("nuisSet5b"));
    modelConfigBG.SetNuisanceParameters(*ws.set("nuisSet5bBG"));
  } else {
    std::cerr << "HELP!" << std::endl;
    exit(1);
  }

  modelConfig.SetPdf(*ws.pdf("model"));
  modelConfig.SetParametersOfInterest(*ws.set("POI"));
  modelConfig.SetObservables(*ws.var("mjjj"));
  ws.var("xs")->setVal(0.0);
  //ws.pdf("model")->fitTo(*ws.data("DataToFit"));
  ws.pdf("model")->fitTo(*ws.data("data"), RooFit::Range(XMIN,XMAX), RooFit::Extended(kTRUE));
  //ws.factory("poiAndNuis::set(POI,nuisSet5b)");
  RooArgSet* poiAndNuis = new RooArgSet("poiAndNuis");
  //poiAndNuis->add(*modelConfig.GetNuisanceParameters());
  poiAndNuis->add(*modelConfig.GetParametersOfInterest());
  modelConfig.SetPriorPdf(*ws.pdf("prior5bWthXS"));
  modelConfig.SetSnapshot(*poiAndNuis);
  //modelConfig.SetGlobalObservables( *poiAndNuis );
  modelConfig.SetGlobalObservables( RooArgSet() );
  delete poiAndNuis;


  modelConfigBG.SetPdf(*ws.pdf("background_model"));
  modelConfigBG.SetParametersOfInterest(*ws.set("POI"));
  modelConfigBG.SetObservables(*ws.var("mjjj"));
  modelConfigBG.SetPriorPdf(*ws.pdf("prior5b"));
  ws.pdf("background_model")->fitTo(*ws.data("data"));
  ws.var("xs")->setVal(0.0);
  poiAndNuis = new RooArgSet("poiAndNuisBG");
  //poiAndNuis->add(*modelConfigBG.GetNuisanceParameters());
  poiAndNuis->add(*modelConfigBG.GetParametersOfInterest());
  ws.pdf("background_model")->fitTo(*ws.data("data"), RooFit::Range(XMIN,XMAX), RooFit::Extended(kTRUE));
  //modelConfigBG.SetGlobalObservables( *poiAndNuis );
  modelConfigBG.SetGlobalObservables( RooArgSet() );
  modelConfigBG.SetSnapshot(*poiAndNuis);
  delete poiAndNuis;

  // import the modelconfig
  ws.import(modelConfig);
  ws.import(modelConfigBG);



  // Just some print statements
  ws.Print();
  modelConfig.Print();
  modelConfigBG.Print();
  ws.var("nbkg")->Print();
  ws.var("nbkgM0")->Print();
  ws.var("nbkgS0")->Print();
  //ws.var("p0")->Print();
  //ws.var("p0M0")->Print();
  //ws.var("p0S0")->Print();
  ws.var("p1")->Print();
  ws.var("p1M0")->Print();
  ws.var("p1S0")->Print();
  ws.var("p2")->Print();
  ws.var("p2M0")->Print();
  ws.var("p2S0")->Print();
  ws.var("lumi")->Print();
  ws.var("lumiM0")->Print();
  ws.var("lumiS0")->Print();
  ws.var("acceptance")->Print();
  ws.var("xs")->Print();
  ws.var("sigMass")->Print();
  ws.var("sigWidth")->Print();





  std::vector<float> Limits;
  if (Section == -2) {
    // Do nothing right now
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
    Can.SetLogy(1);
    Can.SaveAs(TString("Fit_")+label+".pdf");
    ws.var("nbkg")->Print();

    Limits = DoFit(ws, modelConfig, modelConfigBG, label, method, Section, MAXXS);
  } else {
    // Setup some stuff for toy MC experiments TOY, everyone likes TOYS, but not to be toyed with.
    // I can get you a toy, believe me.  There are ways, Dude.


    float const rand_NBGThisPE = NBGFromFit + (ws.var("nbkgS0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1);

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

    int NTest = RooRandom::randomGenerator()->Poisson(
      1.*(NBGFromFit + (ws.var("nbkgS0")->getVal() - 1.0) * RooRandom::randomGenerator()->Gaus(0,1))
        );
    TH1F* hPEF = GetPE3Param(NTest, rand_p1, rand_p2, label);
    RooDataHist hPE("DataToFit", "dataset with x", *ws.var("mjjj"), hPEF);

    //RooDataHist* hPE = PE->binnedClone("DataToFit");
    printf("Data vs PE: %12.1f %12.1f\n", NBGFromFit, hPE.sum(false));

    // Just for debug
    if (false) {
      TCanvas CanPE;
      CanPE.cd();
      RooPlot* PEplot = ws.var("mjjj")->frame();
      hPE.plotOn(PEplot);
      PEplot->Draw();
      CanPE.SaveAs(TString("PE_")+label+".pdf");
    }

    ws.import(hPE);
    if (false) {
      ws.var("nbkg")->Print();
      TCanvas Can;
      Can.cd();
      RooPlot* datafit = ws.var("mjjj")->frame();
      ws.data("DataToFit")->plotOn(datafit);
      //ws.pdf("model")->fitTo(*ws.data("DataToFit"));
      ws.pdf("model")->plotOn(datafit);
      datafit->Draw();
      Can.SaveAs(TString("Fit_")+label+".pdf");
      ws.var("nbkg")->Print();
    }
    Limits = DoFit(ws, modelConfig, modelConfigBG, label, method, Section, MAXXS);

    delete hPEF;

  }



  return Limits;
}


int main (int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [Section] (if -1, mass id)" << std::endl;
    return 1;
  }

  float const StepSize  =  10;

  int const Section = atoi(argv[2]);
  float const BeginMass = argc == 4 ?  250 + atof(argv[3])*StepSize : 250;
  std::cout << "BeginMass: " << BeginMass << std::endl;

  TString const InFileName = argv[1];

  int const Method      =  2;
  int const Systematics =  6;

  // Set the roostats random seed based on secton number..fine..
  RooRandom::randomGenerator()->SetSeed(121723*(Section+2));
  gRandom->SetSeed(711423*(Section+2));


  // Setup output file
  TString OutName;
  if (Section == -1) {
    OutName = TString::Format("Limits_Data_%i.dat", (int) BeginMass);
  } else {
    OutName = TString::Format("PELimits_%i.dat", Section);
  }
  FILE* Out = fopen(OutName.Data(), "w");
  if (Out == NULL) {
    std::cerr << "ERROR: cannot open output file: " << OutName << std::endl;
    exit(1);
  }



  // Print header to outfile


  if (Section == -1) {
    std::vector<float> Limits = RunMultiJetsRooStats(InFileName, BeginMass, Method, Systematics, Section);
    fprintf(Out, "%10E\n", BeginMass);
    fprintf(Out, "%10E\n", Limits[0]);
    fprintf(Out, "%10E\n", Limits[1]);
    fprintf(Out, "%10E\n", Limits[2]);
    fprintf(Out, "%10E\n", Limits[3]);
    fprintf(Out, "%10E\n", Limits[4]);
    fprintf(Out, "%10E\n", Limits[5]);
  } else {
    std::cout << "Oh no, what's going on" << std::endl;
    throw;
  }

  return 0;
}
