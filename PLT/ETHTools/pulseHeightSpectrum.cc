#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TTimer.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

#include <fstream>>
#include <iostream>
#include <iomanip.h>
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include<cmath>
#include<cstdlib>

#include "pulseHeightSpectrum.h"

// construnctor makes nothing but initialising the Histos
pulseHeightSpectrum::pulseHeightSpectrum(int rN){

  runNumber = rN;
  sprintf(path,"/home/l_tester/log/bt05r%06d/",runNumber);
  //sprintf(path,"/data/disk1/rohe/source_test_2008/rawdata/bt05r%06d/",runNumber);  
  sprintf(fileName,"%srun_%06d.root",path, runNumber);
  cout << "Analysis of file "<< fileName << endl;
  
  sprintf(histoFileName,"%sspectra_%06d.root",path,runNumber);
  cout << "Histogram file name: "<< histoFileName << endl;  
  
  sprintf(maskFileName,"%spixelMask.dat",path);
  cout << "Mask file name: "<< maskFileName << endl;    
  
  smallSigCut = 100;
  
  initHistograms();

}

// ====================================================================================
void pulseHeightSpectrum::initHistograms(){

   cout << "init histograms\n";

   
   // containing info of all hit pixels
   
   char hisName[50], hisTitle[100];
   sprintf(hisName,"pOccupancy_%d",runNumber);
   pOccupancy = new TH2D(hisName,"Number of Hit in Pixel",52, 0, 52, 80, 0, 80); 

   sprintf(hisName,"pHighChargeMap_%d",runNumber);
   pHighChargeMap = new TH2D(hisName,"Map of Pixels with Charge from 210 to 250 (Vcal)", 52, 0, 52, 80, 0, 80);

   sprintf(hisName,"pOccupancyVcalNan_%d",runNumber);   
   pOccupancyVcalNan = new TH2D(hisName,"Number of NaN Values returned from Cal-Fct",52, 0, 52, 80, 0, 80);  
   sprintf(hisName,"pAverageSignalADC_%d",runNumber);         
   pAverageSignalADC = new TH2D(hisName,"Average Signal in ADC counts",52, 0, 52, 80, 0, 80);
   sprintf(hisName,"pAverageSignalVcal_%d",runNumber);            
   pAverageSignalVcal = new TH2D(hisName,"Average Signal in Vcal units",52, 0, 52, 80, 0, 80);
   sprintf(hisName,"pNumberOfHitPixels_%d",runNumber);               
   pNumberOfHitPixels = new TH1D(hisName,"Number of hit pixels per Trigger", 50, 0, 50);
   sprintf(hisName,"pPulsHeightDistributionADC_%d",runNumber);               
   pPulsHeightDistributionADC = new TH1D(hisName,"Pulse height Distribution ADC counts", 401, -500, 1500);
   sprintf(hisName,"pPulsHeightDistributionVcal_%d",runNumber);               
   pPulsHeightDistributionVcal= new TH1D(hisName,"Pulse height Distribution in Vcal units", 301, -100, 1400);
   sprintf(hisName,"pAverageSignalDistributionADC_%d",runNumber);               
   pAverageSignalDistributionADC = new TH1D(hisName,"Distribution of pixel's average signal in ADC counts", 401, -500, 1500);   
    sprintf(hisName,"pAverageSignalDistributionVcal_%d",runNumber);               
   pAverageSignalDistributionVcal = new TH1D(hisName,"Distribution of pixel's average signal in Vcal units", 751, -100, 1400);  

   sprintf(hisName,"pMask_%d",runNumber);               
   pMask = new TH2I(hisName,"Pixel Mask",52, 0, 52, 80, 0, 80);     
	                          
   // containing info on clusters
   sprintf(hisName,"cNumberOfClusters_%d",runNumber);                  
   cNumberOfClusters  = new TH1D(hisName,"Number of Clusters per Trigger",50,0,50);
   sprintf(hisName,"cSize_%d",runNumber);                     
   cSize = new TH1D(hisName,"Cluster Size distribution", 50, 0, 50);

   for (int i=0; i<8; i++){
   
     char hisNameQ[20], hisTitleQ[50];
     char hisNameC[20], hisTitleC[50];
     char hisNameA[20], hisTitleA[50]; 
     char hisNameM[20], hisTitleM[50];      
     char hisNameQM[20], hisTitleQM[50];    
             
     switch (i){
     
       case 0: 
         sprintf(hisNameQ,"cQ%d_%d",i,runNumber);
         sprintf(hisTitleQ,"Charge Distribution for all Clusters Run %d", runNumber);
    
         sprintf(hisNameC,"cCentre%d_%d",i,runNumber);
         sprintf(hisTitleC,"Position of all Clusters Run %d", runNumber);
 
         sprintf(hisNameA,"cAverange%d_%d",i,runNumber);
         sprintf(hisTitleA,"Averange Charge of all Clusters Run %d",runNumber);  
	     
         sprintf(hisNameQM,"cQMasked%d_%d",i,runNumber);
         sprintf(hisTitleQM,"Charge Distribution for all masked Clusters Run %d", runNumber);
 	 
	 sprintf(hisNameM,"cCentreMasked%d_%d",i,runNumber);                            
	 sprintf(hisTitleM,"Position of all masked Clusters Run %d",runNumber);                            	 
         break;
       
       case 1:
       case 2:
       case 3:
       case 4:
       case 5:
         sprintf(hisNameQ,"cQ%d_%d",i,runNumber);
         sprintf(hisTitleQ,"Charge Distribution for Clusters of %d Pixels Run %d", i, runNumber);
    
         sprintf(hisNameC,"cCentre%d_%d",i,runNumber);
         sprintf(hisTitleC,"Position of Clusters with %d Pixels Run %d", i, runNumber);
 
         sprintf(hisNameA,"cAverange%d_%d",i,runNumber);
         sprintf(hisTitleA,"Averange Charge of Clusters with %d Pixels Run %d", i, runNumber);     
	 
         sprintf(hisNameQM,"cQMasked%d_%d",i,runNumber);
         sprintf(hisTitleQM,"Charge Distribution for masked Clusters of %d Pixels Run %d", i, runNumber);
 	 
	 sprintf(hisNameM,"cCentreMasked%d_%d",i,runNumber);                            
	 sprintf(hisTitleM,"Position of all masked Clusters with %d Pixels Run %d",i,runNumber);             
         break;
       
       case 6:  
   
         sprintf(hisNameQ,"cQ%d_%d",i, runNumber);
         sprintf(hisTitleQ,"Charge Distribution for Clusters of more than 5 Run %d", runNumber);
    
         sprintf(hisNameC,"cCentre%d_%d",i,runNumber);
         sprintf(hisTitleC,"Position of Clusters with more than 5 Pixels Run %d", runNumber);
 
         sprintf(hisNameA,"cAverange%d_%d",i,runNumber);
         sprintf(hisTitleA,"Averange Charge of Clusters with more than 5 Pixels Run %d", runNumber);  
   
         sprintf(hisNameQM,"cQMasked%d_%d",i, runNumber);
         sprintf(hisTitleQM,"Charge Distribution for masked Clusters of more than 5 Run %d", runNumber);
 	 
	 sprintf(hisNameM,"cCentreMasked%d_%d",i,runNumber);                            
	 sprintf(hisTitleM,"Position of all masked Clusters with more than 5 Pixels Run %d",runNumber);   	     	 
                  
         break;
	 
       case 7: 	 
       
         sprintf(hisNameQ,"cQ%d_%d",i, runNumber);
         sprintf(hisTitleQ,"Charge Distribution for Clusters with 1 or 2 Pixels  Run %d", runNumber);
    
         sprintf(hisNameC,"cCentre%d_%d",i,runNumber);
         sprintf(hisTitleC,"Position of Clusters with 1 or 2 Pixels Run %d", runNumber);
 
         sprintf(hisNameA,"cAverange%d_%d",i,runNumber);
         sprintf(hisTitleA,"Averange Charge of Clusters with 1 or 2 Pixels Run %d", runNumber);  
       
         sprintf(hisNameQM,"cQMasked%d_%d",i, runNumber);
         sprintf(hisTitleQM,"Charge Distribution for masked Clusters with 1 or 2 Pixels  Run %d", runNumber);

	 sprintf(hisNameM,"cCentreMasked%d_%d",i,runNumber);                            
	 sprintf(hisTitleM,"Position of all masked Clusters with 1 or 2 Pixels Run %d",runNumber);   	     	   
	 
	 break;
	 
       default:;      
     } 
     cQ[i]=new TH1D(hisNameQ,hisTitleQ, 501, -100, 2400);
     cCentre[i] = new TH2D(hisNameC,hisTitleC,52,0,52,80,0,80);
     cAverageQ[i] = new TH2D(hisNameA, hisTitleA, 52,0,52,80,0,80);
     
     cQMasked[i]=new TH1D(hisNameQM,hisTitleQM, 501, -100, 2400);     
     cCentreMasked[i] =  new TH2D(hisNameM,hisTitleM,52,0,52,80,0,80);     
     
   }
    
   sprintf(hisName,"cPosLowSignal_%d",runNumber);                     
   sprintf(hisTitle,"Pos of 1Hit Clusters < %d Vcal",smallSigCut);                        
   cPosLowSignal = new TH2D(hisName,hisTitle,52,0,52,80,0,80);
   
   sprintf(hisName,"cPosHighSignal_%d",runNumber);                     
   sprintf(hisTitle,"Pos of 1Hit Clusters > %d Vcal",smallSigCut);                        
   cPosHighSignal = new TH2D(hisName,hisTitle,52,0,52,80,0,80);   
   
   //sprintf(hisName,"cCentreMasked_%d",runNumber);                        
 
   sprintf(hisName,"hFitPar_%d",runNumber);                           
   hFitPar = new TH2D(hisName,"Fit Parameter",8,0,8,4,0,4);
   sprintf(hisName,"hFitErr_%d",runNumber);                              
   hFitErr = new TH2D(hisName,"Errors of Fit Parameter",8,0,8,4,0,4);
   
   cout<< ".... finished\n";
   return;
}

// ====================================================================================

int pulseHeightSpectrum::loadFile(){
    
    // open file
   
    cout << "open data file\n";
    TFile *f = new TFile(fileName);
    if(f->IsZombie()){
      cout<< fileName << " does not exist\n";
      delete(f);
      return 0;
    }   
    
    
    cout << "defining tree\n";
    t = (TTree*)f->Get("t1");

    t->SetBranchAddress("eventNumber",&eventNumber);
    t->SetBranchAddress("realTimeStamp",&realTimeStamp);
    t->SetBranchAddress("PixN",&PixN);
    t->SetBranchAddress("PixCol",&PixCol[0]);
    t->SetBranchAddress("PixRow",&PixRow[0]); 
    t->SetBranchAddress("PixAna",&PixAna[0]);
    t->SetBranchAddress("PixVcal",&PixVcal[0]);        
    //t->SetBranchAddress("pixelAnaArray",&pixelAnaArray[0][0]);
    //t->SetBranchAddress("pixelCalArray",&pixelCalArray[0][0]);
    
    //t->SetBranchAddress("CluN",&CluN);    
    //t->SetBranchAddress("CluSize",&CluSize[0]);        
    //t->SetBranchAddress("CluQ",&CluQ[0]);    
    //t->SetBranchAddress("CluCol",&CluCol[0]);        
    //t->SetBranchAddress("CluRow",&CluRow[0]);  
    
    cout<< ".... finished\n";
    return 1;      
        
}
// =======================================================================================
//
// event loop which fills the pulse heigh spectra and create a mask histogramm
//
void pulseHeightSpectrum::createMask(){
    
    cout << "Starting event loop for creating a pixel mask\n";
    int n = t->GetEntries();

    // event loop
    for (int k=0; k<n; k++){
    
      //cout<<"  open event "<< k <<endl;
    
      t->GetEntry(k);
    
      if ((k%10000)==0) cout << k << " Events processed\n";
    
      // fill the histos
      //
      // histos displaying with "hits"
    
      pNumberOfHitPixels->Fill(PixN);
    
      for (unsigned int hitCounter = 0; hitCounter < PixN; hitCounter++){
    
        pOccupancy->Fill(PixCol[hitCounter],PixRow[hitCounter]);  

	if (PixVcal[hitCounter]>210 && PixVcal[hitCounter]<250){
	  pHighChargeMap->Fill(PixCol[hitCounter],PixRow[hitCounter]);
	}

        pAverageSignalADC->Fill(PixCol[hitCounter],PixRow[hitCounter],PixAna[hitCounter]);    
	if (!isnan(PixVcal[hitCounter])){   
          pAverageSignalVcal->Fill(PixCol[hitCounter],PixRow[hitCounter],PixVcal[hitCounter]);
	  pPulsHeightDistributionVcal->Fill(PixVcal[hitCounter]);
        }else{
	  pOccupancyVcalNan->Fill(PixCol[hitCounter],PixRow[hitCounter]);  
	}    
        pPulsHeightDistributionADC->Fill(PixAna[hitCounter]);

      }
 
    }// end event loop
    cout << ".... event loop finished\n";
    
    cout << "Creating histograms and pixel mask\n";
    // Averaging
    
    
    for (int colCounter = 1; colCounter <= 52; colCounter++){
      for (int rowCounter = 1; rowCounter <= 80; rowCounter++){
        double tmpADC, tmpVcal;

	if (pOccupancy->GetBinContent(colCounter,rowCounter) > 0){
	  tmpADC = pAverageSignalADC->GetBinContent(colCounter,rowCounter) / pOccupancy->GetBinContent(colCounter,rowCounter);
	  
	  if (pOccupancy->GetBinContent(colCounter,rowCounter) > pOccupancyVcalNan->GetBinContent(colCounter,rowCounter)){
	    tmpVcal= pAverageSignalVcal->GetBinContent(colCounter,rowCounter) / (pOccupancy->GetBinContent(colCounter,rowCounter) - pOccupancyVcalNan->GetBinContent(colCounter,rowCounter));
	  }else{ 
	    tmpVcal = 0.;
	  }
	}else{ 
          tmpADC  = 0.; 
	  tmpVcal = 0.;
	}
	//cout << "Pix "<< colCounter-1 << " " << rowCounter-1 << " Occ: " << pOccupancy->GetBinContent(colCounter,rowCounter);
	//cout << " Sum ADC: "<< pAverageSignalADC->GetBinContent(colCounter,rowCounter);
	//cout << " Sum Vcal: "<<pAverageSignalVcal->GetBinContent(colCounter,rowCounter);
	//cout << " Av ADC " << tmpADC << " Av Vcal " << tmpVcal << endl;
	pAverageSignalADC->SetBinContent(colCounter,rowCounter,tmpADC);
	pAverageSignalVcal->SetBinContent(colCounter,rowCounter,tmpVcal);

	pAverageSignalDistributionADC->Fill(tmpADC);
	pAverageSignalDistributionVcal->Fill(tmpVcal);	
	
	// scatter plot occupancy vs average singal
	occ[(rowCounter-1)+((colCounter-1)*80)] = pOccupancy->GetBinContent(colCounter,rowCounter);
	aver[(rowCounter-1)+((colCounter-1)*80)] = tmpVcal;
	

      }  // end row
    }    // end col
    
    // create Scatter plot
    OccupancyVsAverageSignal = new TGraph(4160,occ,aver);

    // create a pixel mask

    for (int colCounter = 1; colCounter <= 52; colCounter++){
      for (int rowCounter = 1; rowCounter <= 80; rowCounter++){    
         // dead -> error code = 1
         if (pOccupancy->GetBinContent(colCounter,rowCounter) == 0){
           pMask->SetBinContent(colCounter,rowCounter,1); 
	   cout << "Mask pix "<< colCounter << " " << rowCounter << " --> dead\n";
	   continue;
         }
	 // calibration problematic --> error code 8 or 16

	 
	 //if ((pAverageSignalVcal->GetBinContent(colCounter,rowCounter) < (0.5*pPulsHeightDistributionVcal->GetMean())) ||
	 //    (pAverageSignalVcal->GetBinContent(colCounter,rowCounter) > (2.0*pPulsHeightDistributionVcal->GetMean()))){
	 if ((pAverageSignalVcal->GetBinContent(colCounter,rowCounter) < (0.1*pAverageSignalDistributionVcal->GetMean())) ||
	     (pAverageSignalVcal->GetBinContent(colCounter,rowCounter) > (10.0 *pAverageSignalDistributionVcal->GetMean()))){

	      pMask->Fill(colCounter-1,rowCounter-1,8); 
	      cout << "Mask pix "<< colCounter << " " << rowCounter << " --> Average Signal too different\n";	    
	 }	 
	 if (pOccupancyVcalNan->GetBinContent(colCounter,rowCounter) > (0.01 * pOccupancyVcalNan->GetEntries())){
	   pMask->Fill(colCounter-1,rowCounter-1,16); 
	   cout << "Mask pix "<< colCounter << " " << rowCounter << " --> Vcal too many NaN-Values\n";	    
	 } 
	 
	 // occupancy too high or low
	 int averageCounter = 0;
	 double averageSum = 0, averageOcc = 0;
	 // calculate the average Occupance in the region +/- 2 pix around the pixel
	 for (int avCol = colCounter - 2; avCol <= colCounter + 2; avCol++){
	   for (int avRow = rowCounter - 2; avRow <= rowCounter + 2; avRow++){
	   
	     if ((avRow < 1) || (avRow > 80) || (avCol < 1) || (avCol > 52)) continue;// outside the chip
	     if (pMask->GetBinContent(avCol,avRow) > 0) continue; // already known bad pixel
	     
	     averageCounter++;
	     averageSum += pOccupancy->GetBinContent(avCol,avRow);
	     // edge pixels count double, corner 4 times
	     if (avRow == 80) averageCounter++;
	     if ((avCol == 1) || (avCol == 52)) averageCounter++;
	     if (((avCol == 1) || (avCol == 52)) && (avRow == 80)) averageCounter += 2;
	   
	   }
	 }	 
	 
	 if (averageCounter >0){
	   averageOcc = averageSum / averageCounter; 
	   if (pOccupancy->GetBinContent(colCounter,rowCounter) < (0.5*averageOcc)){ 
	     pMask->Fill(colCounter-1,rowCounter-1,2); // too few hits
	     cout << "Mask pix "<< colCounter << " " << rowCounter << " Occupancy too low\n";	    
	   }
	   if (pOccupancy->GetBinContent(colCounter,rowCounter) > (2.0*averageOcc)){ 
	     pMask->Fill(colCounter-1,rowCounter-1,4); // hit too often 
	     cout << "Mask pix "<< colCounter << " " << rowCounter << " Occupancy too high\n";
	     //cout << "Number of hits: "<< pOccupancy->GetBinContent(colCounter,rowCounter);
	     //cout << " Average of neighbours = "<< averageSum <<" / "<< averageCounter << " = " << averageOcc << endl;
	   }
	       
	 }else{
	   pMask->Fill(colCounter-1,rowCounter-1,2); // all neighbours masked then mask also this guy. Should be impossible
	   //cout << "Mask pix "<< colCounter << " " << rowCounter << " all neigbours masked\n";
         }
	 
      }
    }
    // load Pixel Mask File to mask more pixels
    
    applyMaskFile();
    
    
    cout<< ".... finished\n";
    
    return;
}
// =======================================================================================
//
//  read in a file pixelMask.dat and add to the pMask Histogram 32=masked by file
//
//
void pulseHeightSpectrum::applyMaskFile(){
  int roc, col, row, val=32;	
  char keyWord[100], line[1000];

  ifstream maskFile;
  maskFile.open(maskFileName);
  
  if (maskFile.bad()) 
  {
      cout << "!!!!!!!!!  ----> Could not open file "<<maskFileName<<" to read pixel mask\n";
      return;
  }  
  
  cout << "Reading pixel mask from "<< maskFileName << endl;
  while(maskFile.good()){
       maskFile>>keyWord;
    if (strcmp(keyWord,"#")==0){ 
       maskFile.getline(line,60, '\n');
       cout << "# "<<line << endl;// ignore rows starting with "#" = comment
    }
    else if(strcmp(keyWord,"pix")==0){
       maskFile>>roc>>col>>row;
       cout << "Exclude "<<keyWord<<" "<<roc<<" "<<col<<" "<<row<<endl; 
       if ((col >= 0)&&(col < 52)&&(row >= 0)&&(row < 80)){
         maskPix(col,row, val);
       }else{
         cout << "!!!!!!!!!  ----> Pixel number out of range: "<<keyWord<<" "<<roc<<" "<<col<<" "<<row<<endl;
       }
    }else if(strcmp(keyWord,"col")==0){
       maskFile>>roc>>col;
       cout << "Exclude "<<keyWord<<" "<<roc<<" "<<col<<endl;
       if ((col >= 0)&&(col < 52)){
         maskCol(col, val);      
       }else{
         cout << "!!!!!!!!!  ----> Pixel number out of range: "<<keyWord<<" "<<roc<<" "<<col<<endl; 
       }
    }else if(strcmp(keyWord,"row")==0){
       maskFile>>roc>>row;
       cout << "Exclude "<<keyWord<<" "<<roc<<" "<<row<<endl;
       if ((row >= 0)&&(row < 80)){
         maskRow(row, val);   
       }else{
         cout << "!!!!!!!!!  ----> Pixel number out of range: "<<keyWord<<" "<<roc<<" "<<row<<endl;
       }
    }
    sprintf(keyWord," ");   
  
  }
  
  maskFile.close();
  
  return;
  

}


// =======================================================================================
//
// 2nd event loop which builts the clusters and fills all histogramms 
// containing cluster info
//
void pulseHeightSpectrum::clusterLoop(){

    cout << "Starting 2nd event loop\n";
    int n = t->GetEntries();

    // event loop
    for (int k=0; k<n; k++){
    
      //cout<<"  open event "<< k <<endl;
    
      t->GetEntry(k);
    
      if ((k%10000)==0) cout << k << " Events processed\n";
    
      // fill the histos
      //
      // histos with pixel hits already filled in 1. event loop (createMask) 
      
      //get Clusters
      clu = findCluster();
            
      // Clusters
      
      cNumberOfClusters->Fill(clu.size());      
      
      for (unsigned int cluCounter = 0; cluCounter < clu.size(); cluCounter++){
        if (clu.at(cluCounter).masked >0){
	  //cout<<"Masked not histogramed\n"; 
	  cCentreMasked[0]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	  cQMasked[0]->Fill(clu.at(cluCounter).charge);
	  if ( clu.at(cluCounter).size < 6){ 
	    cQMasked[clu.at(cluCounter).size]->Fill(clu.at(cluCounter).charge);
	    cCentreMasked[clu.at(cluCounter).size]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	    if ( clu.at(cluCounter).size < 3) {
	      cQMasked[7]->Fill(clu.at(cluCounter).charge);
	      cCentreMasked[7]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	    }else{
	      cQMasked[6]->Fill(clu.at(cluCounter).charge);
	      cCentreMasked[6]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	    }
	  }
	  continue;
	}
        cSize->Fill(clu.at(cluCounter).size);
	cQ[0]->Fill(clu.at(cluCounter).charge);
	cCentre[0]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	cAverageQ[0]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row,clu.at(cluCounter).charge);	
	if ( clu.at(cluCounter).size < 6){ 
	  cQ[clu.at(cluCounter).size]->Fill(clu.at(cluCounter).charge);
	  cCentre[clu.at(cluCounter).size]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	  cAverageQ[clu.at(cluCounter).size]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row, clu.at(cluCounter).charge);
	  
	  // fill his Nr 7 with 1+2 hit clusters
	  if ( clu.at(cluCounter).size < 3) {
	    cQ[7]->Fill(clu.at(cluCounter).charge);
	    cCentre[7]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	    cAverageQ[7]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row, clu.at(cluCounter).charge);	  
	  }
	  
	  
	}else{
	  cQ[6]->Fill(clu.at(cluCounter).charge);
	  cCentre[6]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	  cAverageQ[6]->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row, clu.at(cluCounter).charge);
	}
	// low signal histogram
	if ((clu.at(cluCounter).size == 1)&&(clu.at(cluCounter).charge < smallSigCut)){
	  cPosLowSignal->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	}
	// high signal histogram
	if ((clu.at(cluCounter).size == 1)&&(clu.at(cluCounter).charge >= smallSigCut)){
	  cPosHighSignal->Fill(clu.at(cluCounter).col,clu.at(cluCounter).row);
	}	
	
	
	
      }      
      
    }// end event loop
    
    cout <<".... 2nd event loop finished\n";
    
    // Averaging
    for (int colCounter = 1; colCounter <= 52; colCounter++){
      for (int rowCounter = 1; rowCounter <= 80; rowCounter++){
        double help1;
	
	for (int nClu = 0; nClu < 7; nClu++){
	  if (cCentre[nClu]->GetBinContent(colCounter,rowCounter) != 0){
	    help1 = cAverageQ[nClu]->GetBinContent(colCounter,rowCounter) / cCentre[nClu]->GetBinContent(colCounter,rowCounter);
	  }else{
	    help1 = 0;
	  }
	  cAverageQ[nClu]->SetBinContent(colCounter,rowCounter,help1);
	}
      }  // end row
    }    // end col
    
    return;
}

// =======================================================================================
vector<cluster> pulseHeightSpectrum::findCluster(){
  /* returns clusters with local coordinates and layer IDs	added 
	  the layermap is set at construction time
	*/



  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)
  // written by Wolfram 2005
  
  int fCluCut = 2;
  double fAnaMin = -500; //  no threshold
  
  pixel pb[1000];
  
  vector<cluster> v;
  if(PixN==0) return v;
  int* gone = new int[PixN];
  int* layer = new int[PixN];  
  for(unsigned int i=0; i<PixN; i++){
	 gone[i]=0;
	 // create an array of pixels pb
	 pb[i].col = PixCol[i];
	 pb[i].row = PixRow[i];
	 pb[i].ana = PixAna[i];
	 pb[i].Vcal= PixVcal[i];
  }
  unsigned int seed=0;
  while(seed<PixN){
    // start a new cluster
    cluster c;
    c.vpix.push_back(pb[seed]); gone[seed]=1;
    c.charge=0.; c.size=0; c.col=0; c.row=0;
	
    // let it grow as much as possible
    int growing;
    do{
      growing=0;
      for(unsigned int i=0; i<PixN; i++){
		  if(!gone[i]){
			 for(unsigned int p=0; p<c.vpix.size(); p++){
				int dr = c.vpix.at(p).row - pb[i].row;
				int dc = c.vpix.at(p).col - pb[i].col;
				if(    (dr>=-fCluCut) && (dr<=fCluCut) 
					 && (dc>=-fCluCut) && (dc<=fCluCut) )
				{
				  c.vpix.push_back(pb[i]); gone[i]=1;
				  growing=1;
				  break;//important!
				}
			 }
		  }
      }
    }while(growing);

    // added all I could. determine position and append it to the list of clusters
	 int nBig=0;
	 int clusterMasked = 0;
    for(vector<pixel>::iterator p=c.vpix.begin();  p!=c.vpix.end();  p++){
      double Qpix=p->Vcal;
      c.charge+=Qpix;
      c.col+=(*p).col*Qpix;
      c.row+=(*p).row*Qpix;
      if((*p).Vcal>fAnaMin){nBig++;}
      
      // if p is an edge pixel
      if (((*p).col <= 0)||((*p).col >= 51)||((*p).row <= 0)||((*p).row >= 79)) {
         //mask cluster
	 clusterMasked = 1;
	 //cout<<"Edge Pixel\n";
      }else{
        // check if vpix is = or adjecent to a masked pixel
        for (int iCol = ((*p).col - 1); iCol <= ((*p).col + 1); iCol++){
          for (int iRow = ((*p).row - 1); iRow <= ((*p).row + 1); iRow++){
	    
	    if (pMask->GetBinContent(iCol+1, iRow+1) > 0){
	      clusterMasked = 1; 
	      //cout<<"Next to masked pixel\n";
	    }
	  }
	}
      }
    } 
    c.size=c.vpix.size();
    c.masked=clusterMasked; 
    if(!(c.charge==0)){
	c.col=c.col/c.charge;
	c.row=c.row/c.charge;
    }else{
	c.col=(*c.vpix.begin()).col;
	c.row=(*c.vpix.begin()).row;
	cout << "pulseHeightSPectrum::findCluster>  cluster with zero charge" << endl;
    }
    //if (clusterMasked>0) cout << "Masked a cluster\n";
    if(nBig>0){
	v.push_back(c);
    }
    //look for a new seed
    while((++seed<PixN)&&(gone[seed]));
  }
  // nothing left,  return clusters
	 delete layer;
	 delete gone;
  return v;
}
// =======================================================================================
void pulseHeightSpectrum::fitLandau(){

  LangauFitter *lF;
  //fit Landau/Gauss curves to cQ[0,1,2,3]
  //larger clusters do not show nice Landaus
  lF = new LangauFitter();
  
  for (int i=0; i<8;i++){
    
    double lowLimit = 0, highLimit = 0;
  
    if (i==0) {lowLimit = 175.; highLimit = 400.;}
    if (((runNumber == 4463)||(runNumber == 4462)||(runNumber == 4461))
       &&(i==1)) 
    {
      lowLimit = 140.; 
      highLimit = 300.;
    }// determined by Hand
    
    lF->Fit(cQ[i], lowLimit, highLimit);
    
    lanGauFit[i]=lF->GetFitFct();
    cout<<"For Clusters of size "<<i<<endl;
    for (int ii=0; ii<4; ii++){  
      fitParameters[i][ii] = lF->GetParameter(ii);
      fitErrors[i][ii]= lF->GetError(ii);
      cout <<" Parameter Nr "<< ii<<" "<< fitParameters[i][ii] << " +/- " << fitErrors[i][ii] << endl;
    }
    cout << endl;
  }
  
  return;
}


// =======================================================================================

void pulseHeightSpectrum::writeHistograms(){

   TFile *fOut = new TFile(histoFileName,"RECREATE");
   
   pOccupancy->Write();

   pHighChargeMap->Write();

   pOccupancyVcalNan->Write();
   pAverageSignalADC->Write();          
   pAverageSignalVcal->Write();
   pAverageSignalDistributionADC->Write();
   pAverageSignalDistributionVcal->Write();	   
   pNumberOfHitPixels->Write();
   pPulsHeightDistributionADC->Write();
   pPulsHeightDistributionVcal->Write();
   
   char graphName[100];
   sprintf(graphName,"OccupancyVsAverageSignal_%d",runNumber);
   OccupancyVsAverageSignal->Write(graphName);
   
   pMask->Write();
   
   cNumberOfClusters->Write();
   cSize->Write();
   
   for (int nClu = 0; nClu < 8; nClu++){
     cQ[nClu]->Write(); 
     cCentre[nClu]->Write();
     cAverageQ[nClu]->Write();
     lanGauFit[nClu]->Write();
     
     cQMasked[nClu]->Write();   
     cCentreMasked[nClu]->Write();     
     
     for (int paraNr = 0; paraNr <4; paraNr++){
       hFitPar->SetBinContent(nClu+1,paraNr+1,fitParameters[nClu][paraNr]);
       hFitErr->SetBinContent(nClu+1,paraNr+1,fitErrors[nClu][paraNr]);
     }
     
   }
      
   // special histograms
   cPosLowSignal->Write();
   cPosHighSignal->Write();   
   
   hFitPar->Write();
   hFitErr->Write();
   
   cout << "Histograms written\n";
   
   fOut->Close();
   delete(fOut);

}



