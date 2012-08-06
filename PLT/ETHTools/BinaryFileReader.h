#ifndef BINARYFILE_H
#define BINARYFILE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <deque>
#include <vector>
//#include "RocGeometry.h"
#include "ConfigReader.h"
#include "PHCalibration.h"

#include "pixelForReadout.h"

using namespace std;


#include "pixelForReadout.h"
class PHCalibration;
class TH1F;
class TH2F;

class BinaryFileReader {
  
 private:
  int fTBMHeader[8];
  int fTBMTrailer[8];
  int fNROC;
  int fLayerMap[16];
  int fNHit;
  int fHitROC[16];
  int fOffs[16];
  static const int fUnInitialized=99999;
  bool fIsModule;
  int fCluCut;
  double fAnaOffset[16];
  double fAnaSlope[16];
  int fAnaMin;
  //PHCalibration *rocPHCal;
  
 
  long long int fTime;
  long long int fTmin;
  long long int fTmax; 
  long long int fTLastCalInject;
  long long int fTLastReset;
  long long int fMostRecentTrigger;
  int fTBMTrigger;
  int fnTmaxError;
  int fTCT;
  int fSyncOk;
  int fResync;
  int fRequireSync;

  struct trigger_t{
	 long long int timeStamp;
	 int type;               // = header (internal or external), 0=bad time-stamp
	 int tbmTrigger;         // = expected TBM event number
	 int reset;              // = reset issued while trigger was on stack
  };
  deque<trigger_t> trigger;   // trigger stack
  long long int fTrigTime;   // most recent trigger
  int fTrigType;


  // levels
  int fUbTBM;
  int fTBM[4+1];
  int fUbROC[16];
  int fROC[16][6+1];
  double fDeconvolution[16];
  double fDeconvolution2[16];
  //RocGeometry fRocGeometry;


  public:
  int fDeconvTest;

  static const int kData           =0x01;
  static const int kExternalTrigger=0x02;
  static const int kInternalTrigger=0x04;
  static const int kReset          =0x08;
  static const int kCalInject      =0x10;
  static const int kTBMReset       =0x20;
  static const int kInfiniteRO     =0x40;
  static const int kOvflw          =0x80;

  // level histogram stuff
  static const int nBinLH=512;
  static const int LHMin=-1024;
  static const int LHMax=1024;
  TH1F* hUBTBM;
  TH1F* hLVLTBM;
  TH1F* hUBBROC[16];
  TH1F* h3rdClk[16];
  TH1F* hADROC[16];
  TH1F* hPHROC[16];// raw data without event filter
  TH1F* hPHVcalROC[16];// filled together with Maps

  TH1F* hLVLTBMUsed;
  TH1F* hADROCUsed[16];
  TH2F* hDeconv[16];

  // hit Maps
  TH2I *hModMapCal;
  TH2I *hRocMapCal[16];
  TH2I *hModMapExt;
  TH2I *hRocMapExt[16];
  TH2I *hModMapInt;
  TH2I *hRocMapInt[16];
  TH2I *hModMap;
  TH2I *hRocMap[16];
  TH1F *hNHit;
  TH1F *hNHitRoc[16];
  TH1F *hGino;
  
  // correlation histograms
  TH2I *hCorrelationRow[5][5];
  TH2I *hCorrelationCol[5][5];

 public:

  BinaryFileReader(const char* p, const char* f="mtb.bin",const char *layermap="0",
		   const char* tag="",
		   ConfigReader *cfg=0,  const char* cfgtag="",
		   const int levelMode=0,
		   const int trVal=60);
  ~BinaryFileReader();
  void init();
  //  void close();     // do once at the end , before program exit
  int  open();
  int  readRecord();              // read one record(=header+data) from file (Beat's format with TB header)
  int  readDataEvent();          // read next data event ( do not return triggers)
  int  readGoodDataEvent();       // same but skip bad events

  unsigned short readBinaryWord();     // read two chars and swap order
  void  nextBinaryHeader();            // starting from present position, find next header
  int  decodeBinaryData();             // called after header has been found => time and fData
  //int  getType();
  int  getNHit(){return fNHit;}
  int  getNROC(){return fNROC;}
  int* getTBMHeader(){return fTBMHeader;}
  void updateHistos();
  Double_t findLowestValue(TH1F* h, Float_t threshold=0);
  int  findLevels(TH1F* h, int n0, float* level, int algorithm=1);
  void Levels();
  void updateLevels();
  int  decode(int adc, int nLevel, int* level);
  void  decodePixels();
  void dump(int level=0);
  void printTrailer();
  long long getBC(){return fTime;}
  long long getTrigBC(){return fTrigTime;}
  int getTrigType(){return fTrigType;}
  int  getTBMTrigger();
  int  getTBMStatus();
  //void getHitsObsolete(int *buf);
  //float colToX(int col);
  //float rowToY(int row);
  //void toLocal(pixel& p);
  //vector<cluster> getHits();
  pixel pb[1001];
  //pixel* getPixels(){return pb;};//obsolete
  vector<pixel> getPixels(){vector<pixel> p; for(int i=0; i<fNHit; i++){p.push_back(pb[i]);};return p;};
  //void getPixels(pixel* pbuf){for(int i=0; i<fNHit; i++){pbuf[i]=pb[i];}};//obsolete
  void fillPixelMaps();
  void fillCorrelationPlots();
  bool convertDcolToCol(const int dcol,const int pix,
			int & colAdd, int & rowAdd) const;
  void printRunSummary();
  void printRunSummary2();
  
  void setClusterCut(int cut){fCluCut=cut;}
  void setAnaMin(int amin){fAnaMin=amin;}
  void setAnaOffset(int n,double* r){for(int i=0;i<n;i++){fAnaOffset[i]=r[i];}}
  void setAnaSlope(int n,double* r){for(int i=0;i<n;i++){fAnaSlope[i]=r[i];}}
  void setDeconvolution(double r){for(int i=0; i<16; i++){fDeconvolution[i]=r;}}
  void setDeconvolution(int n, double* r){
     for(int i=0; i< n; i++){fDeconvolution[i]=r[i]; fDeconvolution2[i]=r[i];}
  }
  void setDeconvolution2(int n, double* r){for(int i=0; i< n; i++){fDeconvolution2[i]=r[i];}}
  //void setPHCalibration(PHCalibration* c){fPHcal=c;}
  int  getTriggerStack(){return trigger.size();}
  void printHighEffPixels(int nevent, float thresh, TH2I** h);
  char* msgId();


 private:
  unsigned long fUpperTime, fLowerTime;
  int        fEvent;     // Event counter
  int        fMaxEvent;  // Event number limitation
  int        fHeader, fNextHeader;
  int        fEOF;
  int        fBufferSize;
  static const int NUM_DATA=10000;
  static const int MAX_PIXELS=1000;
  int        fBuffer[NUM_DATA];
  int        fData[NUM_DATA];
  ifstream   *fInputBinaryFile;
  char       fInputFileName[1000];
  char       fPath[1000];
  char       fTag[20];
  char       fLevelFileName[1000];
  char       msgBuf[1010];
  int        fBadTrailer;
  int        fWriteNewLevelFile;
  int        fLevelMode;
  ConfigReader* fcfg;
  char fcfgtag[100];
  // run statistic

  int fnRecord;
  int fnData;
  int fnDataWithHits;
  int fnTrig;
  int fnTrigExternal;
  int fnTrigInternal;
  int fnReset;
  int fDtReset;
  long long int fTlastReset;
  int fnOvflw;
  int fnInfiniteRO;
  int fnCalInject;
  int fnCalInjectHistogrammed;
  int fnBadTrailer;
  int fnCorrupt;
  int fNoTokenPass;
  int fnNoTokenPass;
  int fnInvalidAddress;
  int fNTS5;
  int fnResync;
  int fTruncated;
  int fnTruncated;

  int trimVal;
  PHCalibration *fPHcal;

 public:
  int  eof() {return fInputBinaryFile->eof();}
  int  getOverFlowCount(){return fnOvflw;};
  double getTriggerRate(){return float(fnTrig)/((fTmax-fTmin)*25e-9);}
  int getCalInjectCount(){return fnCalInject;};
  void requireSync(int k){fRequireSync=k;}
  void printPixel(int col, int row);
  void printVcalPeaks(const char* tag="ph");
};
#endif
