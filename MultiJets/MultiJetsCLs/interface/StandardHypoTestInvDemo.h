#ifndef GUARD_StandardHypoTestInvDemo_h
#define GUARD_StandardHypoTestInvDemo_h


#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"

#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

// internal routine to run the inverter
RooStats::HypoTestInverterResult* RunInverter(RooWorkspace * w, const char * modelSBName, const char * modelBName, const char * dataName,
                                     int type,  int testStatType, int npoints, double poimin, double poimax, int ntoys, 
                                     bool useCls, bool useNumberCounting, const char * nuisPriorName);

void StandardHypoTestInvDemo(const char * fileName =0,
                             const char * wsName = "combined",
                             const char * modelSBName = "ModelConfig",
                             const char * modelBName = "",
                             const char * dataName = "obsData",                 
                             int calculatorType = 0,
                             int testStatType = 3, 
                             bool useCls = true ,  
                             int npoints = 5,   
                             double poimin = 0,  
                             double poimax = 5, 
                             int ntoys=1000,
                             bool useNumberCounting = false,
                             const char * nuisPriorName = 0);
void ReadResult(const char * fileName, const char * resultName, bool useCLs);


#endif
