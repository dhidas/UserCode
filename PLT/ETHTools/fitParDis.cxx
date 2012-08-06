#include <iostream>
#include <TSystem.h>
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"

#include "PHCalibration.h"
#include "PHCalibrationFit.h"

using namespace std;


 
int main(int argc, char **argv){
 
 
  int run = 0, trim = 60;
  char dir[200];
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-r")) {
      run=atoi(argv[++i]);
      cout << "run number: " << run << endl;
    }
    if(!strcmp(argv[i], "-trimVcal")) {
      trim = atoi(argv[++i]);
    }
    if(!strcmp(argv[i], "-dir")) {
      strcpy(dir,argv[++i]);
    }
  }  
 
  
  // reading in the hParOldameters from the "old" Fit
  // and Plot their distributions
  TH1D *hParOld[4];
  double fParOld[52][80][4];
  
  TH1D *hParNew[4];  
  double fParNew[52][80][4];  
  
  hParOld[0] = new TH1D("hParOld0","Distribution of hParOld 0", 500, 0, 1E-2);
  hParOld[1] = new TH1D("hParOld1","Distribution of hParOld 1", 500, 0, 5); 
  hParOld[2] = new TH1D("hParOld2","Distribution of hParOld 2", 500, 200., 800.);  
  hParOld[3] = new TH1D("hParOld3","Distribution of hParOld 3", 500, 500., 1000.);
   
  char path[200],fileOut[200];  
 
  //sprintf(path,"/data/disk1/rohe/source_test_2008/rawdata/bt05r%06d",run);
  sprintf(path,"/home/l_tester/singleROCTesting/output/%s",dir);   
  
  PHCalibration *fPHcal;
  fPHcal = new PHCalibration();
  fPHcal->LoadFitParameters(path,trim,0);//second argument is trim vcal value
  
  
  
  for(int colCounter = 0; colCounter < 52; colCounter++){
    for(int rowCounter = 0; rowCounter < 80; rowCounter++){
      for(int parCounter = 0; parCounter < 4; parCounter++){
      
         hParOld[parCounter]->Fill(fPHcal->GetFitPar(0,colCounter,rowCounter,parCounter));
	 fParOld[colCounter][rowCounter][parCounter] = fPHcal->GetFitPar(0,colCounter,rowCounter,parCounter);
	 
      }
    }
  }
  

  /* 
  // now get the "raw" data, redo the fits and compare
  PHCalibrationFit *fPhFit;
  fPhFit = new PHCalibrationFit(3);
  fPhFit->FitAllCurves(path,1);
  for(int colCounter = 0; colCounter < 52; colCounter++){
    for(int rowCounter = 0; rowCounter < 80; rowCounter++){
      for(int parCounter = 0; parCounter < 4; parCounter++){
         hParNew[parCounter]->Fill(fPHcal->GetFitPar(0,colCounter,rowCounter,parCounter));
	 cout << "first line" << endl;
	 fParNew[colCounter][rowCounter][parCounter] = fPHcal->GetFitPar(0,colCounter,rowCounter,parCounter);
	 //GetFitPar(0,colCounter,rowCounter,parCounter);
	 cout << "go to next loop" << endl;
      }
    }
  }
  */
  
  
  
  // write the file
  sprintf(fileOut,"%s/fitParameterDistributionRun_%d.root",path,trim);
  cout << "Output file: "<<fileOut<<endl;
  TFile f(fileOut,"RECREATE");
  TCanvas c1;
  c1.Divide(2,2);
  for(int parCounter = 0; parCounter < 4; parCounter++){
    c1.cd(parCounter+1);
    hParOld[parCounter]->Draw();
    hParOld[parCounter]->Write();
  } 
  f.Close();
  char gifOut[200];
  sprintf(gifOut, "%s/fitParameterDistribution_%d.gif", path, trim);
  c1.Print(gifOut);
  
  return 0;
  
}
