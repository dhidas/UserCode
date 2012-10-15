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

namespace RooStats { 

   class HypoTestInvTool{

   public:
      HypoTestInvTool();
      ~HypoTestInvTool(){};

      HypoTestInverterResult * 
      RunInverter(RooWorkspace * w, 
                  const char * modelSBName, const char * modelBName, 
                  const char * dataName,
                  int type,  int testStatType, 
                  bool useCLs, 
                  int npoints, double poimin, double poimax, int ntoys, 
                  bool useNumberCounting = false, 
                  const char * nuisPriorName = 0);



      void
      AnalyzeResult( HypoTestInverterResult * r,
                     int calculatorType,
                     int testStatType, 
                     bool useCLs,  
                     int npoints,
                     const char * fileNameBase = 0 );

      void SetParameter(const char * name, const char * value);
      void SetParameter(const char * name, bool value);
      void SetParameter(const char * name, int value);
      void SetParameter(const char * name, double value);

   private:

      bool mPlotHypoTestResult;
      bool mWriteResult;
      bool mOptimize;
      bool mUseVectorStore;
      bool mGenerateBinned;
      bool mUseProof;
      bool mRebuild;
      int     mNWorkers;
      int     mNToyToRebuild;
      int     mPrintLevel;
      double  mNToysRatio;
      double  mMaxPoi;
      std::string mMassValue;
      std::string mMinimizerType;                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()

   };

} // end namespace RooStats


#endif
