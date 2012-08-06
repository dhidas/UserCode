#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <stdlib.h>

#include "LangauFitter.h"

using namespace std;

// Define the pixel structure for raw data readout.
struct pixel{
     int col;
     int row;
     int ana;
     float Vcal;
};

struct cluster{
     vector<pixel> vpix;
     int size, masked;
     float charge;
     float col,row;
};
   
class pulseHeightSpectrum{

  private:
  
   // objects
   
   TTree *t;
   
   char path[100];
   char fileName[200];
   char histoFileName[200];
   char maskFileName[200];
   int runNumber;
   
   // Quantities returned for each event
   unsigned int eventNumber;
   double realTimeStamp;
   unsigned int PixN, PixCol[1000], PixRow[1000];
   int PixAna[1000];
   float PixVcal[1000];
   
   int smallSigCut;
   

   vector<cluster> clu;
   
    
   
   // histograms
   //
   // containing info of all hit pixels
   TH2D *pOccupancy,                 // number if hits in each pixel

        *pHighChargeMap,              // map of pixels with charge in "extra" peaks

        *pOccupancyVcalNan,          // sometimes the Cal-function returns a nan-Value
	                             // therefore the occupancy for Vacl values is lower 
				     // in some pixels 
        *pAverageSignalADC,          // average signal ADC counts
        *pAverageSignalVcal;         //                Vcal
   TH1D *pNumberOfHitPixels;         // number of hit pixels per trigger
   TH1D *pPulsHeightDistributionADC, // puls height distribution in ADC counts
   	*pPulsHeightDistributionVcal,//                             Vcal
        *pAverageSignalDistributionADC, // average signal per pixel ADC counts
        *pAverageSignalDistributionVcal;//                          Vcal				     
	
   double occ[4160], aver[4160];	
   TGraph *OccupancyVsAverageSignal;	// scatter plot to separate hits from noise
	
				    
   TH2I *pMask;			     // mask of all pixels =0: pix ok
                                     //                    >0: do not use	
				    
   // containing info on clusters
   TH1D *cNumberOfClusters;
   TH1D *cSize;
   TH1D *cQ[8]; // clusterQ[0] --> charge distribution all clusters
                      // clusterQ[1] --> charge distribution of clusters size 1
		      // clusterQ[5] --> charge distribution of clusters size 5
		      // clusterQ[6] --> charge distribution of clusters size >5
		  
   TH2D *cCentre[8];
   TH2D *cAverageQ[8];
   
   TH2D *cCentreMasked[8];// same for masked clusters
   TH1D *cQMasked[8];  
   
   TH2D *cPosLowSignal, // position of 1-hit clusters < 350 Vcal
        *cPosHighSignal;
   
   // Fitfunctions
   TF1 *lanGauFit[8];
   //Fitparameters
   double fitParameters[8][4];
   double fitErrors[8][4];
   TH2D *hFitPar, *hFitErr;
   
   
   // methods
   public:
   
   pulseHeightSpectrum(int runNumber);
   int loadFile();
   void createMask();  // 1st event loop
   void clusterLoop(); // 2nd event loop
   void fitLandau();
   void writeHistograms();
   
   private:
      
   void initHistograms();

   void applyMaskFile();
   void maskPix(int col, int row, int val){pMask->Fill(col+1, row+1, val);};
   void maskRow(int row, int val){for(int colCounter=0;colCounter<52;colCounter++) pMask->Fill(colCounter, row, val);};
   void maskCol(int col, int val){for(int rowCounter=0;rowCounter<80;rowCounter++) pMask->Fill(col, rowCounter, val);};
   vector<cluster> findCluster();
   
    
     










};
