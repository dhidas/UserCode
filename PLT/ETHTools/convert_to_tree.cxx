#include "TApplication.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPaveLabel.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"


//#include "Viewer.h"
//#include "EventView.h"
#include "LangauFitter.h"
//#include "TGClient.h"
#include "PHCalibration.h"
#include "ConfigReader.h"
#include "BinaryFileReader.h"
#include<stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

// pixel hit and cluster struct
#include "pixelForReadout.h"


using namespace std;
//TCanvas* gCanvas;
 
/**********************************************************************
 *       program for offline analysis of sept 06 CERN  testbeam data  *
 *       udated and simplyfied Aug-?? 2008 for analising data taken   *
 *       summer 08 with the source setup                              *
 *                                                                    *
 *       author: Tilman using Wolfram's code                          *
 *                                                                    *
 *                                                                    *
 *  usage:                                                            *
 *        ./r -r <run-number>                                         *
 *  Options:                                                          *
 *     -v        verbose mode                                         *
 *     -l        bootstrap address levels (re-run without -l later)   *
 *     -s        <number of events> "sample" end after n events       *
 *     -t        <trim value> trim value of the chip used to          *
 *               create the correct file name for the ph-calibration  *
 *               data file. If no number is given 60 is assumed       *
 *                                                                    *
 *   address level decoding is based on the files                     *
 *         levels-module.dat                                          *
 *  and    levels-roc.dat                                             *
 *                                                                    *
 *                                                                    *
 *********************************************************************/

bool fexists(const char* filename)
{
  ifstream ifile(filename);
  bool qexists = ifile;
  ifile.close();
  return qexists;
} 
 
 int main(int argc, char **argv)
{

  //cout << "Size of long long int: "<< sizeof (long long int) << "byte" << endl;// output 8 byte

  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " [InFile] [RootOutFile] [TextOutFile]" << std::endl;
    return 1;
  }
  int run=0;
  int findlevels=0;
  int verbose=0;
  int sample=0;
  char* mtbfilename=argv[1];
  char* rootfileName=argv[2];
  char* textfileName=argv[3];
  char fConfigfile[101]="";
  int trimValueForCalibrationFile = 60;

  
  // -- command line arguments

  run = 1;
  
  //sprintf(mtbfilename,"mtb.bin");
  cout << "data file = "<<  mtbfilename << endl;
  
  ConfigReader *fConfig;
  sprintf(fConfigfile, "config.dat");  
  fConfig=new ConfigReader(fConfigfile,"config.dat");
  
  // create the output textfile
  if(fexists(textfileName)){
    cout<<"TextOutFile already exists. Choose a different name.\n";
    return 0;
  }
  FILE * textFile = fopen(textfileName,"w");
  cout << "rootfile = " << textfileName << endl;  
  

  // create the tree to fill
  //
  if(fexists(rootfileName)){
    cout<<"RootOutFile already exists. Choose a different name.\n";
    fclose(textFile);
    remove(textfileName);
    return 0;
  }
  TFile tf(rootfileName,"RECREATE");
  cout << "rootfile = " << rootfileName << endl;  

 
  // define all Variables to put into the tree
  vector<pixel> pixels;
  //const int nRow = 80, nCol = 52;    
  
  
  int eventNumber;
  //long long int timeStamp;  
  double realTimeStamp;
  //int   pixelAnaArray[nCol][nRow];// init with -500
  //float pixelCalArray[nCol][nRow];// init with -20
  
  unsigned int PixN, PixRow[1000], PixCol[1000];
  int PixAna[1000];
  float PixVcal[1000];
  
  //unsigned int CluN, CluSize[1000];
  //float CluQ[1000], CluCol[1000], CluRow[1000];
  
  
  //vector<cluster> clusters;  
      
  TTree t1("t1","Events taken");
  t1.Branch("eventNumber", &eventNumber, "eventNumber/i");
  t1.Branch("realTimeStamp", &realTimeStamp, "realTimeStamp/D");   
  t1.Branch("PixN", &PixN, "PixN/i"); 
  t1.Branch("PixRow", &PixRow, "PixRow[PixN]/i");
  t1.Branch("PixCol", &PixCol, "PixCol[PixN]/i");  
  t1.Branch("PixAna", &PixAna, "PixAna[PixN]/I");    
  t1.Branch("PixVcal", &PixVcal, "PixVcal[PixN]/F");  
    
  //t1.Branch("pixelAnaArray", &pixelAnaArray,"pixelAnaArray[52][80]/I" );
  //t1.Branch("pixelCalArray", &pixelCalArray,"pixelCalArray[52][80]/F" );

  //t1.Branch("CluN", &CluN, "CluN/i"); 
  //t1.Branch("CluSize", &CluSize, "CluSize[CluN]/i");//Cluster size  
  //t1.Branch("CluQ", &CluQ, "CluQ[CluN]/F");//charge
  //t1.Branch("CluCol", &CluCol, "CluCol[CluN]/F");//CoG Col
  //t1.Branch("CluRow", &CluRow, "CluRow[CluN]/F");//CoG Row       
  
  BinaryFileReader *fRoc, *fRocLv;
  
  // re-do level files if requested
  if(findlevels){
         cout << "Find the levels\n";
	 fRocLv=new BinaryFileReader("",mtbfilename,"0","RtbTemp", fConfig,"roc",1, trimValueForCalibrationFile);
	 //fRocLv->requireSync(0);//no synchronisation anyway
	 if((fRocLv->open()==0) ){
		while( (fRocLv!=NULL)&&(!fRocLv->eof()) ){fRocLv->readRecord();}
	 }
	 fRocLv->updateLevels();
	 fConfig->rewrite();
	 delete fRocLv;
  }
  

  
  fRoc=new BinaryFileReader("",mtbfilename,"0", Form("roc%05d",run),
			    fConfig,"roc",0, trimValueForCalibrationFile);

  if( !(fRoc->open()==0) ){
	 fRoc=NULL;
  }
  
  fRoc->setClusterCut(1);// no empty pixels within a cluster
  
  //==========Event loop ==============
  
  int rStat=0;  // 0 means no more data, 1=good >1=corrupt

  if( fRoc ) rStat=fRoc->readGoodDataEvent();
  bool more=  rStat ;
  
  int event_counter = (more?1:0);//counter in the event loop
  
  int readRoc=1;

  
  while (more){
        
    bool bRoc=fRoc && !fRoc->eof();
  
    // sort out what kind of data we have here
    if( bRoc ){
      readRoc=1;
    }else{
      // no data, what are we doing here?????
      cout << "EventReader::loop> nothing to be read, bailing out" << endl;
      exit(1);
    }
    
    if(verbose && fRoc && (0)){ // never show
      cout << "EventReader> ";
     
      if(readRoc){
	char cal=' ';
	if (fRoc->getTBMStatus()&0x02){ cal='c';}
	cout << Form("ROC:%10lld(x%2x)  TBM=%3x%C",
		     fRoc->getTrigBC(),fRoc->getTrigType(),fRoc->getTBMTrigger(),cal);
      }else{
	cout << "                       ";
      }
      cout << endl;
    }
    
    // ==========================================
    // get all the info written into the tree
    

    pixels = fRoc->getPixels();

    // reinitialze the arrays
    /*
    for (int ii=0; ii< nCol; ii++){
      for (int iii=0; iii< nRow; iii++){
        pixelAnaArray[ii][iii] = -500;
        pixelCalArray[ii][iii] = -20.;
      }
    }  
    */  
    // fill the variables for the tree
    eventNumber = event_counter;
    realTimeStamp = fRoc->getBC() * 2.5E-8;
    
    PixN = pixels.size();
    
    for (unsigned int countHits = 0; countHits < pixels.size(); countHits++){
    
      PixCol[countHits] = pixels.at(countHits).col;
      PixRow[countHits] = pixels.at(countHits).row;
      PixAna[countHits] = pixels.at(countHits).ana;
      PixVcal[countHits]= pixels.at(countHits).anaVcal;
      // coordinates returned from binaryFileReader are 0-52 and 0-79
      //pixelAnaArray[pixels.at(countHits).col][pixels.at(countHits).row] = pixels.at(countHits).ana;
      //pixelCalArray[pixels.at(countHits).col][pixels.at(countHits).row] = pixels.at(countHits).anaVcal;

      //Print to the text file. Choose to call this Channel 1 and ROC 0.
      fprintf(textFile, "1 0 %i %i %i %i\n",PixCol[countHits],PixRow[countHits],PixAna[countHits],eventNumber);

    }

    
    
    //clusters = fRoc->getHits();
    //CluN = clusters.size();
    /*for (unsigned int clusterCounter = 0; clusterCounter < CluN; clusterCounter++){
       CluSize[clusterCounter] = clusters.at(clusterCounter).size;       
       CluQ[clusterCounter]    = clusters.at(clusterCounter).charge;
       CluCol[clusterCounter]    = clusters.at(clusterCounter).col;       
       CluRow[clusterCounter]    = clusters.at(clusterCounter).row;       
    }*/
    

    t1.Fill();
    
    if (verbose && fRoc){ 
      int numberOfPixels = pixels.size();
      //int numberOfClusters = clusters.size();
      cout << numberOfPixels << " pixels hit " << endl;//numberOfClusters << " clusters found "<< endl;
    } 
            
    
    //==============================================
    
    if (event_counter%10000 == 0) cout << "Event: "<< event_counter << endl;
    
    if( fRoc && readRoc ) rStat=fRoc->readGoodDataEvent();
    more=(rStat);
    event_counter++;
    
    if ((sample > 0)&&(event_counter >= sample)) more = false; // interrupt if sample mode is called
    
  }
  //==========End event loop===========
  cout << "----------------------------------------------------------"<<endl;
  cout << "Events read: "<<event_counter<<endl;
  cout << "----------------------------------------------------------"<<endl;
  
  //fRoc->updateLevels();  
  //fConfig->rewrite();
 
  fRoc->printRunSummary();
  fRoc->printRunSummary2();
  
  fclose(textFile);
  tf.Write();
  tf.Close();
 
 }
