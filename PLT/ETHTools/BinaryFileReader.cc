#include <iostream>
#include <TSystem.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#include "TMath.h"


using namespace std;

#include "BinaryFileReader.h"
#include "PHCalibration.h"
#include "pixelForReadout.h"

//constructor used by fast_ana
BinaryFileReader::BinaryFileReader(const char* p, const char* f,// p=path, f=file
                                   const char *layermap, 
				   const char *tag, 
				   ConfigReader* cfg, const char *cfgtag,
				   int levelMode,
				   int trVal) {
  // new constructor, layer-info passed as an array
  // layer number per roc (readout-order, not I2C chipID)
  // array terminated by -1
  
  sprintf(fInputFileName,"%s/%s",p,f);
  sprintf(fPath,"%s",p); 
  
  fNROC=0;
  // is not needed reordering 1 chip makes no sense
  while( (fNROC<16) && !(layermap[fNROC]=='\0') )
  {
	 fLayerMap[fNROC]=int(layermap[fNROC]-'0');
	 fNROC++;
  }
  strncpy(fTag,tag,20);
  fcfg=cfg;
  strncpy(fcfgtag,cfgtag,20);
  fLevelMode=levelMode;  // 0 (default)=use levels from config, 1=re-determine levels
  
  trimVal = trVal;
  
  // PHCalibration
  //rocPHCal = new PHCalibration();
  //rocPHCal->LoadFitParameters(fPath,fNROC); 
  
  init();
}



void BinaryFileReader::init(){
  fPHcal = NULL;
  fMaxEvent = 9999999;
  fHeader = fNextHeader = -1;
  fEOF = 0;
  // init run statistics
  fnRecord            = 0;
  fnTrig              = 0;
  fnTrigExternal      = 0;
  fnTrigInternal      = 0;
  fnTrig              = 0;
  fnData              = 0;
  fnDataWithHits      = 0;
  fnReset             = 0;
  fDtReset            = 0;
  fTlastReset         = 0;
  fnOvflw             = 0;
  fnInfiniteRO        = 0;
  fnBadTrailer        = 0;
  fnCorrupt           = 0;
  fnNoTokenPass       = 0;
  fnInvalidAddress    = 0;
  fnCalInject         = 0;
  fnCalInjectHistogrammed=0;
  fnTruncated         =0;
  //fnResync    =0;
  fMostRecentTrigger=-1;
  fnTmaxError=0;
  fTime=0;
  fTmin=1000000000000LL;
  fTmax=0;
  fNTS5=0;
  fAnaMin=-500;
  fCluCut =2;
  fDeconvTest=-1;
  fTCT=-1;
  fTBMTrigger=-1;
  //fSyncOk=0;       // force resync
  //fResync=0;       // wait until after next reset
  //fRequireSync=1;

  // (un-)initialize levels
  fUbTBM=fUnInitialized;
  for(int roc=0; roc<fNROC; roc++){
    fUbROC[roc]=fUnInitialized;
         // values will be overwirtten with numbers
	 // from config.dat if available
	 fDeconvolution[roc]=0.02;
	 fDeconvolution2[roc]=0;
	 fAnaOffset[roc]=-200;
	 fAnaSlope[roc]=1;
  }

  // get some configs, if available
  if(fcfg){
    fcfg->get(Form("%s.clusterCut",fcfgtag),fCluCut,2);
    unsigned int nroc;
    nroc=fNROC;
    fcfg->geta(Form("%s.anaOffset",fcfgtag),nroc,fAnaOffset);	      
    fcfg->geta(Form("%s.anaSlope", fcfgtag),nroc,fAnaSlope);
    fcfg->geta(Form("%s.deconvolution",  fcfgtag),nroc,fDeconvolution);
    fcfg->geta(Form("%s.deconvolution2", fcfgtag),nroc,fDeconvolution2);
    // get levels from config file ultra black and levels in one row per TBM/ROC
    int buf[10]; // assemble UB and Levels here
    for (int j=0; j<10; j++) buf[j] = fUnInitialized;// set to very large value
    fcfg->geta(Form("%s.tbmLevels", fcfgtag),6,buf);
    fUbTBM=buf[0];     
    for(int i=0; i<5; i++) fTBM[i]=buf[i+1];
    
    //again set buf to fUnInitialized
    for (int j=0; j<10; j++) buf[j] = fUnInitialized;
    for(int roc=0; roc<fNROC; roc++){
      fcfg->geta(Form("%s.rocLevels", fcfgtag),roc,8,buf);
      fUbROC[roc]=buf[0];  for(int i=0; i<7; i++){ 
	fROC[roc][i]=buf[i+1];
      }
    }
    
  }
  
  
  // open PHCalibration
  fPHcal = new PHCalibration();
  fPHcal->LoadFitParameters(fPath,trimVal,0);
  
  // create control histograms

  char name[100], title[100]; char module[10]="";
 
  /*
  cout << strlen(fTag) << endl;
  char tag[21]="";
  if(strlen(fTag)>0){
	 snprintf(tag,20,"_%s",fTag);  // prepend an underscore if we have a tag
  }
  */
  
  sprintf(name,"%subtbm%s",module,fTag);  
  hUBTBM=new TH1F(name,"tbm ultrablack",nBinLH,LHMin,LHMax);

  sprintf(name,"%slvltbm%s",module,fTag);
  hLVLTBM=new TH1F(name,"tbm levles",nBinLH,LHMin,LHMax);

  sprintf(name,"%sNHits%s",module,fTag);
  hNHit=new TH1F(name,"number of hits",40,-0.5,39.5);

  sprintf(name,"%sGino%s",module,fTag);
  hGino=new TH1F(name,"unexpected hits",10,-0.5, 9.5);

  for(int roc=0; roc<fNROC; roc++){

    sprintf(name,"%sDeconv%d%s",module,roc,fTag);
    sprintf(title,"%sDeconv%d,",module,roc);
    hDeconv[roc]=new TH2F(name,title,350,-400,900, 350, -400,900); 
  
    sprintf(name,"%subb_%d%s",module,roc,fTag);
    sprintf(title,"ultrablack/black roc %d",roc);
    hUBBROC[roc]=new TH1F(name,title,nBinLH,LHMin,LHMax);
    
    sprintf(name,"%slvl_%d%s",module,roc,fTag);
    sprintf(title,"levels roc %d",roc);
    hADROC[roc]=new TH1F(name,title,nBinLH,LHMin,LHMax);

    sprintf(name,"%sph_%d%s",module,roc,fTag);
    sprintf(title,"pulse height %d",roc);
    hPHROC[roc]=new TH1F(name,title,nBinLH,LHMin,LHMax);
  
    sprintf(name,"%s3rdclk_%d%s",module,roc,fTag);
    sprintf(title,"3rd clk %d",roc);
    h3rdClk[roc]=new TH1F(name,title,nBinLH,LHMin,LHMax);

    sprintf(name,"%sphVal_%d%s",module,roc,fTag);
    sprintf(title,"pulse height in vcal DAC units %d",roc);
    hPHVcalROC[roc]=new TH1F(name,title,nBinLH,-100., 800.);

    // IMPORTANT
    sprintf(name,"%sMap_%d%s",module,roc,fTag);
    sprintf(title,"hitmap %d",roc);
    hRocMap[roc]=new TH2I(name,title,52,0,52,80,0,80);
     

    // IMPORTANT
    sprintf(name,"%sNHits_%d%s",module,roc,fTag);
    sprintf(title,"number of hits %d",roc);
    hNHitRoc[roc]=new TH1F(name,title,20,-0.5,19.5);

    // possible to fill?
    sprintf(name,"%sMapCal_%d%s",module,roc,fTag);
    sprintf(title,"Calinject Map %d",roc);
    hRocMapCal[roc]=new TH2I(name,title,52,0,52,80,0,80);

    // MPORTANT: Particles (??)
    sprintf(name,"%sMapExt_%d%s",module,roc,fTag);
    sprintf(title,"Ext Trigger Map %d",roc);
    hRocMapExt[roc]=new TH2I(name,title,52,0,52,80,0,80);

    // What is this ??
    sprintf(name,"%sMapInt_%d%s",module,roc,fTag);
    sprintf(title,"Int Trigger Map %d",roc);
    hRocMapInt[roc]=new TH2I(name,title,52,0,52,80,0,80);

  }
  
  
  //cout << " ..... done " << endl;

  // clone histograms to show used decoding levels
  double somelargenumber=100000.;   // dummy for filling histos
  hLVLTBMUsed= new TH1F(*hLVLTBM);
  hLVLTBMUsed->SetName(Form("tbmlv%s",fTag));
  for(int i=0; i<5; i++){ 
    hLVLTBMUsed->Fill(fTBM[i],somelargenumber);
  }
  for(int roc=0; roc<fNROC; roc++){
    hADROCUsed[roc]=new TH1F(*hADROC[roc]);
    hADROCUsed[roc]->SetName(Form("roclvl_%d%s",roc,fTag));
    for(int i=0; i<7; i++){
      hADROCUsed[roc]->Fill(fROC[roc][i],somelargenumber);
    }
  }

  // initialize geometry utitilty
  //fRocGeometry=RocGeometry(27*0.0150, 40.5*0.0100);
}

// ----------------------------------------------------------------------
BinaryFileReader::~BinaryFileReader(){
  // delete the biggest chunks
  for(int i=0; i<fNROC; i++){
	 delete hRocMap[i];
	 delete hRocMapInt[i];
	 delete hRocMapExt[i];
	 delete hRocMapCal[i];
  }
}

// ----------------------------------------------------------------------
int BinaryFileReader::open() {

  fInputBinaryFile = new ifstream(fInputFileName);

  if (fInputBinaryFile->is_open()) {

    cout << "--> reading from file " << fInputFileName << endl;

	 unsigned short word=readBinaryWord();
	 while( !((word&0xFF00)==0x8000) && (fEOF==0) ){
		word=readBinaryWord();
	 }

	 if( fEOF ) {
		fNextHeader=-1;
		cout << "--> file contains no headers:" << fInputFileName << endl;
		return 1;
	 }else if(word&0x8000){
		fNextHeader = word & 0x00FF;
		return 0;
	 }else{
		cout << "we should never be here" << endl;
		return 1;
	 }
  } else {
    cout << "--> ERROR: unable to open file " << fInputFileName << endl;
    return 1;
  }
};


// ----------------------------------------------------------------------
unsigned short BinaryFileReader::readBinaryWord() {
 
  if (fInputBinaryFile->eof()) { fEOF = 1; return 0; }
  unsigned char a = fInputBinaryFile->get();
  if (fInputBinaryFile->eof()) { fEOF = 1; return 0; }
  unsigned char b = fInputBinaryFile->get();
  if (fInputBinaryFile->eof()) { fEOF = 1; return 0; }
  unsigned short word  =  (b << 8) | a;

  //  cout << Form("readBinaryWord: %02x %02x word: %04x ", a, b, word) << endl;

  return word;
}


// ----------------------------------------------------------------------
void BinaryFileReader::nextBinaryHeader() {
  /* read data into fBuffer until the next header is encountered
	  afterwards :
	  fNextHeader = type of header encountered (=header of the next data block)
	  fHeader     = Header at the beginning of the block now in fBuffer
	  new triggers are put on the trigger stack 
  */
	  
  if(fEOF) return;
  
  //clear buffer
  for (int i = 0; i < NUM_DATA; ++i) {
	 fBuffer[i] = 0;
	 fData[i]=0;
  }

  // the header has already been read in by the previous call
  // it has been stored in fNextHeader
  fHeader=fNextHeader;
  // get at least three words (=time stamp)
  fBufferSize=0;
  for(fBufferSize=0; fBufferSize<3; fBufferSize++){
	 unsigned short word=0;
	 word=readBinaryWord();
	 if(fEOF) break;
	 fBuffer[fBufferSize]=word;
  }
  if(fEOF){
	 // no more data
	 return ;
  }

  // decode header time if there is one
  if((fHeader>=0)&&(fBufferSize==3)){

    unsigned short t0 = fBuffer[0];
    unsigned short t1 = fBuffer[1];
    unsigned short t2 = fBuffer[2];
    fUpperTime = t0;
    fLowerTime = (t1 << 16) | t2;
    long long int newTime
      =(((long long int) fUpperTime)<<32)+((long long int)fLowerTime);
    
    // sanity checks on the time-stamp, try to catch readout bit errors
    int tsOk=0;
    if (newTime<fTime){    //((newTime<fTime)||(newTime> (fTmax+40000000LL))){
      cout << msgId() << "bad time-stamp " 
	   << Form("%4x:  %4x  %4x  %4x",0x8000|fHeader, t0,t1,t2) << endl;
      cout << "newTime="<<newTime<<"    fTime="<<fTime<<"     fTmax="<<fTmax<<endl;
      newTime=0;
    }else{
      tsOk=1;
    }
    fTime=newTime;

	 
    // look at the stack in case of resets
    if( fHeader & kReset ){
      if (!trigger.empty()){
	//cout << msgId() << "Reset while trigger stack not empty" << endl;
	deque<trigger_t>::iterator it;
	for (it = trigger.begin(); it != trigger.end(); it++){
	  it->reset=1;
	}
      }
      fTLastReset=fTime;
      //  use only 1 roc --> no synchronisation
      //
      // resync on first event after reset if necessary
      // no sycronisation
      //if ((fSyncOk==0)&&(fRequireSync)) {
	//fResync=1;
	//cout << msgId() << "resync requested" << endl;
      //}
    }
    

    // keep track of triggers
    trigger_t thisTrigger;
    if( fHeader & (kExternalTrigger | kInternalTrigger) ){
      
      if(tsOk){
	thisTrigger.timeStamp=fTime;
	thisTrigger.type=fHeader & (kExternalTrigger | kInternalTrigger);
	thisTrigger.tbmTrigger=fTBMTrigger;
	// tbm marks the next n BCs as no token pass
	if(fTime-fTLastReset>8){  // according to TBM documentation
	  thisTrigger.reset=0;
	}else{
	  thisTrigger.reset=1;
	}
      }else{
	thisTrigger.timeStamp=0;
	thisTrigger.type=0;
	thisTrigger.tbmTrigger=-1;
	thisTrigger.reset=0;
      }
      
      // flag internal triggers with a matching cal-inject
      if ( fHeader & kInternalTrigger){
	int dtcal=int(fTime - fTLastCalInject);
	if( (dtcal==fTCT) || ((fTCT==-1)&&(dtcal < 256 )) ){
	  if(fTCT<0){ fTCT=dtcal; }// must be the first one, hope it's right
	  // matches tct, assume it's good, mark the trigger accordingly
	  thisTrigger.type |=kCalInject;
	}
      }
      
      if( (fTime < (fTLastCalInject+4)) ||  (fTime < fMostRecentTrigger+4) ){
	// trigger probably swallowd by FPGA
	cout << msgId() << "swallowed Trigger. Time between most recent trigger and `now': "<< fTime-fMostRecentTrigger<< endl;
      }else{
	trigger.push_back(thisTrigger);
	fMostRecentTrigger=fTime;
	fTBMTrigger++;		
	if(fTBMTrigger==256){fTBMTrigger=0;};
      }
    }
    
    // record min and max time for run duration
    //if((fTime<fTmax)||(fTime> (fTmax+40000000LL))){
    // delete 2nd condition as there are large time gaps between spills:
    if (fTime<fTmax){
      fnTmaxError++;
      if(fnTmaxError>5){
	cout << msgId() << "resetting Tmax from " << fTmax 
	     << " to " << fTime << endl;
	fTmax=fTime;
	fnTmaxError=0;
      }
    }else{
      if((fTime<fTmin)&&(tsOk)) fTmin=fTime;
      if((fTime>fTmax)&&(tsOk)) fTmax=fTime;
      fnTmaxError=0;
    }
    
    
    
    // keep track of the most recent calInject
    if( fHeader == kCalInject ){
      // no OR here, Calinject with simultaneous trigger/reset is not possible
      if( fTime < fMostRecentTrigger+4 ){
	// calinject probably swallowd by FPGA ?
	cout << msgId() << "swallowed calinject " << endl;
      }else{
	fTLastCalInject=fTime;
      }
    }else if (fHeader & kCalInject) {
      cout << msgId() << "  funny calinject header " << Form("%2x",fHeader) << endl;
    }
  }

  
  // get more data until a new header (or the end of file is encountered)    
  // headers are of the format 0x80XX
  // where the 8 bit of XX are
  // bit hex
  //  0     1 data
  //  1     2 ext trigger (CTB)   // new in november/testbeam
  //  2     4 int trigger
  //  3     8 reset
  //  4  0x10 cal
  //  5  0x20 TBM rest
  //  6  0x40 inf readout
  //  7  0x80 data overflow


  while (fEOF==0) {
	 unsigned short word=0;
	 word = readBinaryWord();
	 if (fEOF) break;
	 
	 if( (word&0x8000)==0 ){
		// not a header, keep adding
		if( fBufferSize < NUM_DATA) {
		  fBuffer[fBufferSize++] = word;
		}else{
		  // skip to avoid overrun and warn
		  cout <<  msgId() << "internal buffer overflow" << endl;
		  //fSyncOk=0;// no sychronisation any more
		}
	 }else{
		// header bit was set, was it a valid header?
		if( (word&0x7F00)==0) {
		  fNextHeader = word & 0x00FF;
		  break;
		}else{
		  cout << msgId() 
				 << "illegal header word ignored " << Form("%4x",word) 
				 <<endl;
		  if( fBufferSize < NUM_DATA) {
			 fBuffer[fBufferSize++] = word;
		  }else{
			 // skip to avoid overrun and warn
			 cout << msgId() << "internal buffer overflow" << endl;
		  }
		}
	 }
  }
}

// ----------------------------------------------------------------------
int BinaryFileReader::decodeBinaryData() { 

  int j(0);

  fNHit=0;
  fBadTrailer=0;
  fNoTokenPass=0;
  fTruncated=0;

  for(int roc=0; roc<fNROC; roc++){ fHitROC[roc]=0; }

  if (fHeader > 0) {  
    //cout << Form(" Event at time  %04x/%08x with Header %d", fUpperTime, fLowerTime, fHeader) << endl;
  } else {
    cout << "No valid header, skipping this event" << endl;
    return -1;
  }
  


  for (int i = 3; i < fBufferSize; i++) {
    int value = fBuffer[i] & 0x0fff;
    if (value & 0x0800) value -= 4096;
    fData[i-3] = value;
    ++j;
  }

  fBufferSize -=3;
  
  if (0) {
    for (int i = 0; i < fBufferSize; ++i) { 
      cout << Form(" %04x ", fData[i]);
    }
    cout << endl;
    
    for (int i = 0; i < fBufferSize; ++i) { 
      cout << Form(" %6i ", fData[i]);
    }
    cout << endl;
  }
  

  if(fHeader&kData){
    // level bootstrap mode
    if(fUbTBM==fUnInitialized){
      // guess TBM levels using the first header
      float a=-(fData[0]+fData[1]+fData[2])/3./4.;
      fUbTBM = int( -3*a);
      fTBM[0]= int( -2*a);
      fTBM[1]= int( -0.5*a);
      fTBM[2]= int(  0.5*a);
      //cout << fData[0] << " "  << fData[1] << " " << fData[2] << " " << fData[3] <<endl;
      cout << Form("UB/B(TBM)") << (fData[0]+fData[1]+fData[2])/3.  << " "  << fData[3];
      cout << "  set levels  UB<" << fUbTBM
	   << "   D>"  << fTBM[0]  <<endl;
    }
    

    // tbm header??
    if((fData[0]>fUbTBM)||(fData[1]>fUbTBM)||(fData[2]>fUbTBM)){
      cout << "bad TBM header ? Check levels !!" << endl;
      return 0;
    }else{
      // tbm header ok, copy 
      for(int i=0; i<8; i++){  fTBMHeader[i]=fData[i]; }
    }

    // does a tbm trailer follow right behind?
    int i=8;  // pointer, 8=fist element (UB) of 1. chip
    if((fData[8]<fUbTBM)&&(fData[9]<fUbTBM)&&
       (fData[10]>fTBM[1])&&(fData[10]<fTBM[2])){
      fNoTokenPass=1;
      fnNoTokenPass++;
    }else{
      // chop data along UBs
      for(int roc=0; roc<fNROC; roc++){
	// apply deconvolution for the black level and 3rd hit too
	fData[i+1]+= int (fDeconvolution[roc]*(fData[i+1]-fData[i  ]) + 0.49);
	fData[i+2]+= int (fDeconvolution[roc]*(fData[i+2]-fData[i+1]) + 0.49);
	// ROC header follows, are we in address level bootstrap mode?
	if(fUbROC[roc]==fUnInitialized){
	  cout << Form("UB/B(%3d)",roc) << fData[i] << " "  << fData[i+1];
	  fUbROC[roc] =int(fData[i]*0.8+fData[i+1]*0.2);
	  fROC[roc][0]=int(fData[i]*0.5+fData[i+1]*0.5);
	  cout << "  set levels  UB<" << fUbROC[roc] 
	       << "   D>"  << fROC[roc][0]  <<endl;
	}
	
	fOffs[roc]=i;   // keep pointers to chip data
	i=i+3;          // move on to data
	while((fData[i]>fROC[roc][0])&&(fNHit<MAX_PIXELS)&&(i<fBufferSize)){
	  // not an UB: hit data !!!
	  fNHit++;
	  fHitROC[roc]++;
	  //	cout << "hit roc=" << roc << "   :";
	  // change start
	  for(int k=0; k<6; k++){ 
	    if(k<5){ // address	      
	      fData[i]+= int (fDeconvolution[roc]*(fData[i]-fData[i-1]) +0.49); 
	      // a different way of devoncoltion: 
	      //fData[i]+= int ( fDeconvolution[roc]*fData[i]-fDeconvolution2[roc]*fData[i-1] +0.49); 
	    }else{// k==5: analogue pulse height
	      //deconvolution not necessarry, as already included in Ph-Calibration
	      //fData[i]+= int (fDeconvolution[roc]*(fData[i]-fData[i-1]) +0.49);
	    }
	    i++;
	  }
	  // change end
	}
      }// decode roc data
    }
    
    // TBM trailer
    for(int k=0; k<8; k++){  fTBMTrailer[k]=fData[i++]; }
    if(   (fTBMTrailer[0]<fUbTBM)
	  &&(fTBMTrailer[1]<fUbTBM)
	  &&(fTBMTrailer[2]>fUbTBM)&&(fTBMTrailer[2]<-fUbTBM) ){
      fBadTrailer=0;
      if(fNoTokenPass){ printTrailer();}
    }else{
      //so, we expected to find the trailer here but found something else
      // was this readout truncated by a reset ?
      //for(int i=-10; i<20; i++) cout << i << ")" << fData[fBufferSize-i]<< endl;
      if(    (fData[fBufferSize-8]<fUbTBM)
	     &&(fData[fBufferSize-7]<fUbTBM)
	     &&(fData[fBufferSize-6]>fUbTBM)&&(fData[fBufferSize-6]<-fUbTBM) 
	     ){
	// ok, there is a trailer. fill it and return (e.g. for status decoding), 
	// but mark the event as truncated
	for(int ii=0; ii<8; ii++){ fTBMTrailer[ii]=fData[fBufferSize-8+ii]; }
	fTruncated=1;
	fBadTrailer=0;	
	//change2 start
	if ( i < fBufferSize ){ //Trailer came "too late" (e.g. empty clock cycle?)
		cout << msgId() << "trailer too late" << Form("  next header = %4x",fNextHeader) << endl;	
	}else{ //Trailer came "too early": tuncated event	
		cout << msgId() << "truncated event" << Form("  next header = %4x",fNextHeader) << endl;
	}
	//change 2 end
	//printTrailer();
      }else if (fEOF){
	fBadTrailer=1; // last event of file, don't use, but don't complain
      }else{
	// something else went wrong
	fBadTrailer=1;
	cout << msgId() << "decodeBinaryData: bad trailer" << endl;
	dump(1);
	// force resync, just to be on the safe side
	// fSyncOk=0; // no synchronisation
      }
    }
  }
  return j;
}



// ----------------------------------------------------------------------
void BinaryFileReader::dump(int level){
  cout << endl;
  cout << "-------------------------------" << endl;
  cout << "raw("<< fBufferSize+3 <<"):"<< endl;
  cout << Form("[%4x] " , fHeader) << " ";
  for(int i=0; i<fBufferSize; i++){
    cout << Form("%2x ",fBuffer[i]) << " ";
	 if(i==3){ cout << Form("(%ld)",time);}
  }
  cout << endl;
  if(level>0){
    cout << "decoded:" << endl;
    for(int j=0; j<8; j++){
      cout << fData[j] << " ";
    }
    cout << endl;
    int i=8;
    for(int roc=0; roc<fNROC; roc++){
      i = fOffs[roc];
      cout << roc << ") ";
      cout << fData[i] << " " << fData[i+1] << " " << fData[i+2] << endl;
      i+=3;  // move pointer to hit data 
      for(int hit=0; hit<fHitROC[roc]; hit++){
		  cout << "      ";
		  for(int j=0; j<6; j++){
			 cout << fData[i++] << " ";
		  }
		  cout << endl;
      }
    }
    for(int j=0; j<8; j++){
      cout << fData[i++] << " ";
    }
    cout << endl;
  }
}

// ----------------------------------------------------------------------
void BinaryFileReader::printTrailer(){
  int stat=getTBMStatus();
  /*
  Trailer Status

    Bit 1) Stack Full LSB
    Bit 2) Pre-Calibrate Issued MSB
    Bit 3) Cleared Trigger Counter LSB
    Bit 4) Sync Trigger MSB
    Bit 5) Sync Trigger Error LSB
    Bit 6) Reset ROC MSB

    Bit 7) Reset TBMLSB
    Bit 8) No Token Pass MSB 
  */  cout << fInputFileName << "(" <<fTime << ") ";
  cout << "TBM status ";
  if(stat & 0x01){ cout <<"[Stack Full]  ";}
  if(stat & 0x02){ cout <<"[CalInject]   ";}
  if(stat & 0x04){ cout <<"[Cleared Trigger Counter] ";}
  if(stat & 0x08){ cout <<"[Sync] ";}
  if(stat & 0x10){ cout <<"[Sync Error] ";}
  if(stat & 0x20){ cout <<"[Reset ROC] ";}
  if(stat & 0x40){ cout <<"[Reset TBM] ";}
  if(stat & 0x80){ cout <<"[No Token Pass] ";}
  cout << endl;
}
// ----------------------------------------------------------------------

char *BinaryFileReader::msgId(){
  snprintf(msgBuf,100,"BinaryFileReader:%s (%10lld): ",fInputFileName,fTime);
  return msgBuf;
}


// ----------------------------------------------------------------------
/*
int BinaryFileReader::getType(){
  if((fHeader&kData)&&(fBadTrailer==1)) return 255;  // data with bad trailer 
  if(fNoTokenPass==1)                return   2;  // no token pass
  //  if((fHeader==80)&&(fNextHeader!=80)) return 1;  
  return fHeader;
};

*/

// ----------------------------------------------------------------------
void BinaryFileReader::updateHistos(){
  if (fTBMHeader[0]==0) dump();
  hUBTBM->Fill(fTBMHeader[0]);
  hUBTBM->Fill(fTBMHeader[1]);
  hUBTBM->Fill(fTBMHeader[2]);
  hLVLTBM->Fill(fTBMHeader[3]);
  hLVLTBM->Fill(fTBMHeader[4]);
  hLVLTBM->Fill(fTBMHeader[5]);
  hLVLTBM->Fill(fTBMHeader[6]);
  hLVLTBM->Fill(fTBMHeader[7]);
  hNHit->Fill(fNHit);
  for(int roc=0; roc<fNROC; roc++){
    hNHitRoc[roc]->Fill(fHitROC[roc]);
    int i = fOffs[roc];
    hUBBROC[roc]->Fill(fData[i]);
    hUBBROC[roc]->Fill(fData[i+1]);
    h3rdClk[roc]->Fill(fData[i+2]);
    i+=3;  // move pointer to hit data
    for(int hit=0; hit<fHitROC[roc]; hit++){
		hDeconv[roc]->Fill(fData[i-1],fData[i  ]);// i-1 is either 3rd clk or PH
		hDeconv[roc]->Fill(fData[i  ],fData[i+1]);
		hDeconv[roc]->Fill(fData[i+1],fData[i+2]);
		hDeconv[roc]->Fill(fData[i+2],fData[i+3]);
		hDeconv[roc]->Fill(fData[i+3],fData[i+4]);
		//hDeconv[roc]->Fill(fData[i+4],fData[i+5]);// pulse-height

      hADROC[roc]->Fill(fData[i  ]);
      hADROC[roc]->Fill(fData[i+1]);
      hADROC[roc]->Fill(fData[i+2]);
      hADROC[roc]->Fill(fData[i+3]);
      hADROC[roc]->Fill(fData[i+4]);
      hPHROC[roc]->Fill(fData[i+5]);

      i+=6;  //  next hit
    }
  }
}



// ----------------------------------------------------------------------
int BinaryFileReader::findLevels(TH1F* h, int n0, float* level, int algorithm){
  const int debug=0;
  int npeak=0;  // return number of peaks found
  if (debug) cout << "findLevels>" << endl;
  if(algorithm==1){ // only alorithm == 2 used
/*
    //Use TSpectrum to find the peak candidates
    TSpectrum *s = new TSpectrum(2*n0);
    //Int_t nfound = s->Search(h,1,"new");
    Int_t nfound = s->Search(h,1,"goff");
    Float_t *xpeaks = s->GetPositionX();
    for (int p=0;p<nfound;p++) {
      Float_t xp = xpeaks[p];
      Int_t bin = h->GetXaxis()->FindBin(xp);
      Float_t yp = h->GetBinContent(bin);
      if (yp/TMath::Sqrt(yp) <2) continue;
      level[npeak]=xpeaks[p];
      npeak++;
      //    cout << npeak << " " << xpeaks[p] << endl;
    }
*/
  }else{

    // clustering type peak finding
    const int nmiss=3;  // max number of missing/low  bins within a peak
    const float threshold=h->GetMaximum()/100.+1; // ignore bins below threshold
    int nlow=nmiss+1;  // to make sure the first entry starts a peak
    float s=0;
    float sx=0;
    if (debug) cout << " bin clustering threshold=" << threshold << endl;
    for(int b=0; b<h->GetNbinsX(); b++){
      if (debug) cout << Form("%5d  %10.1f   %10.1f   %3d %3d\n",b,h->GetBinCenter(b),h->GetBinContent(b),npeak,nlow);
      if(h->GetBinContent(b)>threshold){
		  if(nlow>nmiss){
			 // start a new peak
			 s=0;   // reset sums for mean 
			 sx=0;
			 npeak++;
		  }
		  // add to peak
		  s+=float(h->GetBinContent(b));
		  sx+=float(h->GetBinCenter(b))*float(h->GetBinContent(b));
		  level[npeak-1]=sx/s;
		  nlow=0;
      }else{
		  nlow++;
      }
    }
  }
  
  
  /* merge peaks if more than expected */
  while((n0>0)&&(npeak>n0)){
    int m=0;
    for(int i=1; i<(npeak-1); i++){
      if((level[i+1]-level[i])<level[m+1]-level[m]) m=i;
    }
    cout << "merged levels " << level[m] << " " << level[m+1] << endl;
    level[m]=(level[m]+level[m+1])/2.; // better if weighted?
    for(int i=m+1; i<(npeak-1); i++){
      level[i]=level[i+1];
    }
    npeak--;
  }
  return npeak;
}




// ----------------------------------------------------------------------
Double_t BinaryFileReader::findLowestValue(TH1F* h, Float_t threshold){
  const int debug=0;
  for(int b=0; b<h->GetNbinsX(); b++){
    if (debug) cout << Form("%5d  %10.1f   %10.1f\n",b,h->GetBinCenter(b),h->GetBinContent(b));
    if(h->GetBinContent(b)>threshold) return h->GetBinCenter(b);
  }
  return h->GetXaxis()->GetXmax();
}




// ----------------------------------------------------------------------
void BinaryFileReader::Levels(){
  float lubb[100];
  float l[100];
  int n;
   n=findLevels(hUBTBM, 1,l, 2 );
   if(n==1){
     fUbTBM=(int) (0.8*l[0]);
   }else{
     cout << "bad TBM ultrablack" << endl;
   }

   // TBM event counter levels
   n=findLevels(hLVLTBM,  4,l,2);
   if(n==4){
     //if(fUbTBM>l[0])  fUbTBM=l[0];
     fTBM[0]=(int) (0.5*(fUbTBM+l[0]));
     fTBM[1]=(int) (0.5*(l[0]+l[1]));
     fTBM[2]=(int) (0.5*(l[1]+l[2]));
     fTBM[3]=(int) (0.5*(l[2]+l[3]));
     fTBM[4]=(int) (0.5*(l[3]+LHMax));
   }else{
     cout << "problem with TBM levels" << endl;
     cout << " n= " << n << endl;
   }

   for (int roc=0; roc<fNROC; roc++){
     //cout <<" level finding roc " << roc << endl;
     n=findLevels(hUBBROC[roc], 2,lubb,2);
     if(n==2){
       fUbROC[roc]=(int) (0.8*lubb[0]+0.1*lubb[1]);
       Double_t phmin=findLowestValue(hPHROC[roc]);
       if(phmin <  fUbROC[roc]){
	 cout << "warning: lowest PH < UB cut    PHmin" << phmin
	      << "   Ub= "<< fUbROC[roc] << "    roc= "<< roc << endl; 
       }
     }else{
       cout << "problem with Roc UB  Roc="<< roc;
       cout << "   nlevel = " << n << endl;
     }
     
     n=findLevels(hADROC[roc],6,l,2);
     if(n!=6){
       cout <<" roc " << roc << "  retry level finding " << endl;
       n=findLevels(hADROC[roc],6,l,2);
     }
	  
     if(n==6){
       fROC[roc][0]=(int) (0.5*(lubb[0]+l[0]));
       fROC[roc][1]=(int) (0.5*(l[0]+l[1]));
       fROC[roc][2]=(int) (0.5*(l[1]+l[2]));
       fROC[roc][3]=(int) (0.5*(l[2]+l[3]));
       fROC[roc][4]=(int) (0.5*(l[3]+l[4]));
       fROC[roc][5]=(int) (0.5*(l[4]+l[5]));
       fROC[roc][6]=(int) (0.5*(l[5]+LHMax));
     }else{
       cout << "problem with Roc levels  Roc="<< roc;
		 cout << "   nlevel = " << n << endl;
     }
   }
}



// ----------------------------------------------------------------------
void BinaryFileReader::updateLevels() {
  if(fLevelMode==1){
    Levels();
    int buf[10];
    buf[0]=fUbTBM;  for(int i=0; i<5; i++){ buf[i+1]=fTBM[i];}
    fcfg->updatea(Form("%s.tbmLevels",fcfgtag), 6, buf,"%4d");
    for(int roc=0; roc< fNROC; roc++){
      buf[0]=fUbROC[roc];  for(int i=0; i<7; i++){ buf[i+1]=fROC[roc][i];}
      fcfg->updatea(Form("%s.rocLevels",fcfgtag), roc, 8, buf,"%4d");
    }
  }
}

// ----------------------------------------------------------------------

int BinaryFileReader::decode(int adc, int nlevel, int* level){
  for(int i=0; i<nlevel; i++){
    if ( (adc>=level[i])&&(adc<level[i+1]) ) return i;
  }
  //cout << "level " << adc << " not found" << endl;
  //cout << level[0] << " " << level[1] << " " << level[2] << " " << level[3] << endl;
  if(adc<level[0]) return 0;
  return nlevel;
}

// ----------------------------------------------------------------------

int BinaryFileReader::getTBMTrigger(){
  int value=  
     (decode(fTBMHeader[4],4,fTBM)<<6)
    +(decode(fTBMHeader[5],4,fTBM)<<4)
    +(decode(fTBMHeader[6],4,fTBM)<<2)
    +(decode(fTBMHeader[7],4,fTBM));
  return value;
}

// ----------------------------------------------------------------------

int BinaryFileReader::getTBMStatus(){
  int value=  (decode(fTBMTrailer[4],4,fTBM)<<6)
    +(decode(fTBMTrailer[5],4,fTBM)<<4)
    +(decode(fTBMTrailer[6],4,fTBM)<<2)
             +(decode(fTBMTrailer[7],4,fTBM));
  return value;
}

/*****************************************************************/
//                                        from Daneks DecodeRawPacket
// Convert dcol&pix to col&row  
// Decodeing from "Weber" pixel addresses to rows for PSI46
// dcol = 0 - 25
// pix = 2 - 161 zigzag pattern.
// colAdd = 1-52   ! col&row start from 1
// rowAdd = 1-53
bool BinaryFileReader::convertDcolToCol(const int dcol,const int pix,
                                       int & colAdd, int & rowAdd) const
{
  const int ROCNUMDCOLS=26;
  const int ROCNUMCOLS=52;
  const int ROCNUMROWS=80;
  const int printWarning=1;
  if(dcol<0||dcol>=ROCNUMDCOLS||pix<2||pix>161)
	 {
		if(printWarning){
		  //		  cout << msgId() <<"wrong dcol or pix in user_decode "<<dcol<<" "<<pix<<endl;
		  cout << "wrong dcol or pix in user_decode "<<dcol<<" "<<pix<<endl;
		}
		rowAdd = -1;     // dummy row Address
		colAdd = -1;     // dummy col Address
		return false;
	 }
  
  // First find if we are in the first or 2nd col of a dcol.
  int colEvenOdd = pix%2;  // module(2), 0-1st sol, 1-2nd col.
  colAdd = dcol * 2 + colEvenOdd; // col address, starts from 0
  rowAdd = abs( int(pix/2) - 80);  // row addres, starts from 0
  if( colAdd<0 || colAdd>ROCNUMCOLS || rowAdd<0 || rowAdd>ROCNUMROWS )
	 {
		if(printWarning)
		  {
			 cout <<"wrong col or row in user_decode "<<colAdd<<" "<<rowAdd<<endl;
			 cout << "wrong dcol or pix in user_decode "<<dcol<<" "<<pix<<endl;
		  }
		rowAdd = -1;    // dummy row Address
		colAdd = -1;    // dummy col Address
		return false;
	 }
  return true;
}


// ----------------------------------------------------------------------
/*
void BinaryFileReader::getHitsObsolete(int *buf){
  int k=0;
  for(int roc=0; roc<fNROC; roc++){
    int j=fOffs[roc]+3;
    for(int i=0; i<fHitROC[roc]; i++){
      //for(int o=0; o<6; o++){cout << fData[j+o] << " ";}; cout << endl;
      //for(int o=0; o<5; o++){cout << decode(fData[j+o],6,fROC[roc]) << " ";}; cout << endl;
      buf[k  ]=roc;
      int dcol =  decode(fData[j  ],6,fROC[roc])*6
  	              +decode(fData[j+1],6,fROC[roc]);
      int pix  =  decode(fData[j+2  ],6,fROC[roc])*6*6
  	         +decode(fData[j+3],6,fROC[roc])*6
  	         +decode(fData[j+4],6,fROC[roc]);
		convertDcolToCol(dcol,pix, buf[k+1], buf[k+2]);
      if(buf[k+1]==-1) fnInvalidAddress++;
      buf[k+3] = fData[j+5];
      j+=6;
      k+=4;
    }
  }
}

*/

// ----------------------------------------------------------------------
// basic conversion from row/column to roc/module coordinates

/*
float BinaryFileReader::colToX(int col){
  const float x0=-54*0.0150/2.;
  if(col==0) return x0;
  if(col==51) return 52*0.0150+x0;
  return col*0.0150+0.0075+x0;
}

float BinaryFileReader::rowToY(int row){
  const float y0=-81*0.0100/2.;
  if(row==79) return 80*0.0100+y0;
  return row*0.0100+y0;
}


void BinaryFileReader::toLocal(pixel& p){
// conversion of roc coordinates to module coordinates 
    if(fIsModule){
      if(p.roc<8){
        p.row=159-p.rowROC;
        p.col=p.roc*52 + p.colROC;
		  p.xy[0]=p.roc*54*0.0150 + colToX(p.colROC);
		  p.xy[1]= 160*0.0100- rowToY(p.rowROC);
      }else{//roc=8..16
        p.row=p.rowROC;
        p.col=(16-p.roc)*52-p.colROC-1;
		  p.xy[0]=(16-p.roc)*54*0.0150-colToX(p.colROC);
		  p.xy[1]= rowToY(p.rowROC);
      }
    }else{// single roc
      p.row=p.rowROC;
      p.col=p.colROC;
		p.xy[0]=colToX(p.colROC);
		p.xy[1]=rowToY(p.rowROC);
    }
}

*/
// ----------------------------------------------------------------------

/*
vector<cluster> BinaryFileReader::getHits(){
  // returns clusters with local coordinates and layer IDs	added 
  //  the layermap is set at construction time
	

  // decodePixels should have been called before to fill pixel buffer pb 


  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  
  vector<cluster> v;
  if(fNHit==0) return v;
  int* gone = new int[fNHit];
  int* layer = new int[fNHit];  
  for(int i=0; i<fNHit; i++){
	 gone[i]=0;
	 // layer[i]=fLayerMap[pb[i].roc];
	 //	 printf("getHits> %d  %d  %d  %d  %d \n",i,pb[i].roc, layer[i], pb[i].col, pb[i].row);
  }
  int seed=0;
  while(seed<fNHit){
    // start a new cluster
    cluster c;
    c.vpix.push_back(pb[seed]); gone[seed]=1;
    c.charge=0.; c.size=0; c.col=0; c.row=0;
	 //c.xy[0]=0;
	 //c.xy[1]=0.;
	 //c.layer=layer[seed];
    // let it grow as much as possible
    int growing;
    do{
      growing=0;
      for(int i=0; i<fNHit; i++){
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
    for(  vector<pixel>::iterator p=c.vpix.begin();  p!=c.vpix.end();  p++){
      double Qpix=p->anaVcal;
      c.charge+=Qpix;
      c.col+=(*p).col*Qpix;
      c.row+=(*p).row*Qpix;
      //c.xy[0]+=(*p).xy[0]*Qpix;
      //c.xy[1]+=(*p).xy[1]*Qpix;
      if((*p).anaVcal>fAnaMin){nBig++;}
    } 
    c.size=c.vpix.size();
	 if(!(c.charge==0)){
		c.col=c.col/c.charge;
		c.row=c.row/c.charge;
		//c.xy[0]/=c.charge;
		//c.xy[1]/=c.charge;
	 }else{
		c.col=(*c.vpix.begin()).col;
		c.row=(*c.vpix.begin()).row;
		//c.xy[0]=(*c.vpix.begin()).xy[0];
		//c.xy[1]=(*c.vpix.begin()).xy[1];
		cout << "BinaryFileReader::GetHits>  cluster with zero charge" << endl;
	 }
	 if(nBig>0){
		v.push_back(c);
	 }
    //look for a new seed
    while((++seed<fNHit)&&(gone[seed]));
  }
  // nothing left,  return clusters
	 delete layer;
	 delete gone;
  return v;
}
*/
// ----------------------------------------------------------------------
void BinaryFileReader::decodePixels(){
  int k=0;
  for(int roc=0; roc<fNROC; roc++){
    int j=fOffs[roc]+3;
    for(int i=0; i<fHitROC[roc]; i++){
      pb[k].roc=roc;//here roc (according to token path)
      pb[k].ana = fData[j+5];
      pb[k].row=-1;
      pb[k].col=-1;
      int dcol =  decode(fData[j  ],6,fROC[roc])*6
  	              +decode(fData[j+1],6,fROC[roc]);
      int pix  =  decode(fData[j+2],6,fROC[roc])*6*6
  	              +decode(fData[j+3],6,fROC[roc])*6
  	              +decode(fData[j+4],6,fROC[roc]);

      //if(!convertDcolToCol( dcol,pix, pb[k].colROC, pb[k].rowROC )){
      if(!convertDcolToCol( dcol,pix, pb[k].col, pb[k].row)){
		  fnInvalidAddress++;
		}
      

      
      if(fPHcal!=NULL){
                  //cout << "Calibration called\n";
		  //if ((pb[k].colROC != -1) && (pb[k].rowROC != -1)){ 
		  if ((pb[k].col != -1) && (pb[k].row != -1)){
		   	 //pb[k].anaVcal = fPHcal->GetVcal(pb[k].ana, roc, pb[k].colROC, pb[k].rowROC);	  
			 pb[k].anaVcal = fPHcal->GetVcal(pb[k].ana, roc, pb[k].col, pb[k].row);
			 hPHVcalROC[roc]->Fill(pb[k].anaVcal);
		  }
      }else{
        //cout << "standard cal called\n";
	pb[k].anaVcal=(pb[k].ana-fAnaOffset[roc])*fAnaSlope[roc];
	if(pb[k].anaVcal<=0){
	  pb[k].anaVcal=1;
        }
	
	
	//cout << fAnaOffset[roc] <<  " " << fAnaSlope[roc] 
	//     <<" " << pb[k].ana << " " << pb[k].anaVcal << endl;
	
	hPHVcalROC[roc]->Fill(pb[k].anaVcal);
      }
      j+=6;
      k++;
    }
  }

// not needed as only one set of coordianted used (single chip only
/*
  // convert to local/module coordinates
  for (int i=0; i<fNHit; i++){
	 if(fIsModule){// is never true as we have only 1 chip
		fRocGeometry.getModLocal(pb[i].roc, pb[i].colROC, pb[i].rowROC, pb[i].xy);
		if(pb[i].roc<8){
        pb[i].row=159-pb[i].rowROC;
        pb[i].col=pb[i].roc*52 + pb[i].colROC;
      }else{
        pb[i].row=pb[i].rowROC;
        pb[i].col=(16-pb[i].roc)*52-pb[i].colROC-1;
      }
	 }else{
		pb[i].col=pb[i].colROC;
		pb[i].row=pb[i].rowROC;
		fRocGeometry.getRocLocal(pb[i].colROC, pb[i].rowROC, pb[i].xy);
	 }
  }
*/
}





// ----------------------------------------------------------------------
void BinaryFileReader::fillPixelMaps(){
  // fill maps
  if(fTrigType&kCalInject){fnCalInjectHistogrammed++;}
  int n2020=0;
  int ndcol10=0;
  for (int i=0; i<fNHit; i++){
	
	//pb[i].col=pb[i].colROC;
	//pb[i].row=pb[i].rowROC;
	
	//cout << "================================\n";
	//cout << "Fill Roc map"<<endl;
	
	// in the following always changed colROC --> col und rowROC -> row
	
	hRocMap[pb[i].roc]->Fill(pb[i].col, pb[i].row);
	
	// are those needed???
	if(fTrigType&kCalInject)       hRocMapCal[pb[i].roc]->Fill(pb[i].col, pb[i].row);
	if(fTrigType&kExternalTrigger) hRocMapExt[pb[i].roc]->Fill(pb[i].col, pb[i].row);
	if(fTrigType&kInternalTrigger) hRocMapInt[pb[i].roc]->Fill(pb[i].col, pb[i].row);
	
	
	
  }
  
  
  
  if((fTrigType&kCalInject)&&(n2020==0)){
	 hGino->Fill(ndcol10);
  }
}

// ----------------------------------------------------------------------
int BinaryFileReader::readRecord() {
  /* read the next record, i.e. header or data, and keep the
	  trigger stack up-to-date
	  histogram update is called for data records,
	  no address decoding or clustering
  */

  // slurp in data until the next Event header (0x8X) is found
  nextBinaryHeader();
  
  // decode and classify  header
  int words=0;
  if(fBufferSize>=3){
    words  = decodeBinaryData();
    if ((fHeader & kData)&&(fBadTrailer==0)&&(fTruncated==0)){
      updateHistos();
    }else if ((fHeader & kData)&&(fBadTrailer==1)){
      fnBadTrailer++;
    }else if ((fHeader & kData)&&(fTruncated==1)){
      fnTruncated++;
    }
  }else{
    fnCorrupt++;
  }


  // stats
  fnRecord++;
  if( fHeader & (kInternalTrigger | kExternalTrigger) ){	 fnTrig++;  }
  if( fHeader & kInternalTrigger ){	 fnTrigInternal++;  }
  if( fHeader & kExternalTrigger ){	 fnTrigExternal++;  }
  if(( fHeader & kData )&&(fBadTrailer==0) ){
	 fnData++;
	 if (fNHit>0) fnDataWithHits++;
  }
  if (fHeader & kReset)     { fnReset++; }
  if (fHeader & kOvflw)     { fnOvflw++; }
  if (fHeader & kInfiniteRO){ fnInfiniteRO++; }
  

  // get trigger time and type for data events
  if(fHeader & kData ){

    // no synchronisation in our case
    /*
    //synchronise trigger counters at the very beginning 
    //probably not necessary
    if (trigger[0].tbmTrigger < 0){ 
    
    	while( (!trigger.empty()) && (trigger[0].reset) ) trigger.erase(trigger.begin());
	if( trigger.empty() ){
	  cout << msgId() << " resynchronization failed!" << endl;
	}else{
	  fResync=0;
	  fSyncOk=1;
	  fnResync++;
	  fTBMTrigger=getTBMTrigger();
	  for(unsigned int k=0; k<trigger.size(); k++){
	    trigger[k].tbmTrigger=fTBMTrigger;
	    fTBMTrigger++; if(fTBMTrigger==256) fTBMTrigger=0;
	  }
	  cout << msgId() <<  "synchronized         tbm counter=" << fTBMTrigger << endl;
        }
    
    }
    */  
  
    // now take the next trigger from the stack
    if(!trigger.empty()){
	fTrigTime=trigger.front().timeStamp;
	if((fTrigTime>fTime)||(fTime-fTrigTime>4000000LL)){
	  cout  <<  fInputFileName
		<< ":   trigger timestamp=" << fTrigTime 
		<< "    data Timestamp=" << fTime << endl;
	  fTrigTime=0;
	  fTrigType=0;
	}
	fTrigType=trigger[0].type;
	
	// only 1 roc no synchronisation -> trigger counting not necesarry
	/*
	// checks if trigger counts of TB and TBM are 
	// equal
	if( !(trigger[0].tbmTrigger==getTBMTrigger()) ){
	  cout << msgId() << "TBM trigger mismatch,   expected " 
	       << trigger[0].tbmTrigger 
	       << "    got " << getTBMTrigger() << endl;
	  //fSyncOk=0;//no synchronisation
	}
	*/
	trigger.erase(trigger.begin());
	// sanity checks on trigger stack
	if (trigger.size()>0){ fNTS5++; } else { fNTS5=0; }
	if( fNTS5>10 ){
	  cout <<  msgId()
	       << ": trigger stack above 0 for > 10 events, resetting " << endl;
	  fNTS5=0;
	  while(!trigger.empty()) {trigger.erase(trigger.begin());}
	}
    }else{
	cout << msgId() 
	     << "Data without pending trigger" << endl;
	fTrigTime=0;
	fTrigType=0;
    }
    
    if(fTrigType&kCalInject) fnCalInject++;
    
  }else{
    // not data
    fTrigTime=0;
    fTrigType=0;
    // limit the stack size
    if(trigger.size()>33) trigger.erase(trigger.begin());
  }

  
  return words;
}



// ----------------------------------------------------------------------
void BinaryFileReader::printHighEffPixels(int nevent, float thresh, TH2I** h){
  for(int roc=0; roc<fNROC; roc++){
	 for(int col=0; col<52; col++){
		for(int row=0; row<80; row++){
		  // root histo bin numbering starts at 1, not at 0
		  int nhit=(int) h[roc]->GetBinContent(col+1, row+1);
		  if(h[roc]->GetBinContent(col+1, row+1)>int(thresh*nevent)){
		    //			 printf("%2d  %2d  %2d   %6f        ( %8d/%8d ) \n",
		    printf("roc%d  %2d  %2d   %6f        ( %8d/%8d ) \n",
				roc,col,row,
				h[roc]->GetBinContent(col+1, row+1)/nevent,
				nhit,nevent);
		  }
		}
	 }
  }
  cout << nevent << " Cal-Inject events" << endl;
}


// ----------------------------------------------------------------------
void BinaryFileReader::printVcalPeaks(const char* tag){
  if(fnCalInjectHistogrammed==0) return;
  float l[100];
  float ph[16];
  for(int roc=0; roc<fNROC; roc++){
	 int n=findLevels(hPHROC[roc], 1, l, 2);
	 if(n==1){
		ph[roc]=l[0];
	 }else{
		ph[roc]=-1000;
	 }
  }

  cout << "@@@" << tag << " ";
  for(int roc=0; roc<fNROC; roc++){	 cout << ph[roc] << " "; }
  cout << endl;
}



// ----------------------------------------------------------------------
int BinaryFileReader::readGoodDataEvent() {
  int stat=readDataEvent();
  while(stat>0){
    if (stat==1) return 1;
    stat=readDataEvent();
  }
  return 0;
}


// ----------------------------------------------------------------------
int BinaryFileReader::readDataEvent() {
  /* get the next data event, return 0 if no more data is available,
     return 1 for good events, >1 for bad events (truncated etc)
  */
  bool good,data;

  // get records until the next data record
  do{
  
    //cout << "Read Record"<<endl;
    readRecord();
    data=(fHeader & kData)==kData;
    good=data && (fBadTrailer==0)  && (fNoTokenPass==0) 
      && (fTruncated==0) ;//&& (fTrigType>0) 
    /*
    cout << msgId() << "readDataEvent  " << data << " " << good 
	 << " badTrailer(ok=0)=" << fBadTrailer 
	 << " noTokenPass(ok=0)=" << fNoTokenPass
	 << " truncated(ok=0)=" << fTruncated
	 << " type(ok>0)=" << fTrigType
	 << endl;
    */
  }while( (!data) && (fEOF==0) );
  
  if ( fEOF !=0 ) {
    // no more data, don't call me again
    return 0;  // 
  }else if ( data ){
    if( good ){
      // good data record
      //cout<<"decode pixels \n";
      decodePixels();
      //cout<<"fill histograms\n";
      fillPixelMaps();

      return 1;
    }else{
      // bad data record
      return 3;
    }
  }else{
    cout << msgId() << "readDataEvent : no clue how we got here" << endl;
    return 0; 
  }
}




// ----------------------------------------------------------------------
void BinaryFileReader::printRunSummary() { 
  float dt=(fTmax-fTmin)*25e-9;
  cout << "----------------------------------------------------------" << endl;
  cout << "time " << dt << " seconds" << endl;
  cout << "triggers: " << fnTrig << "   overflows: " << fnOvflw << "  corrupt: " << fnCorrupt << "   good  readouts: " << fnData << "   w/hits: " << fnDataWithHits << endl;
  cout << "invald pixel addresses " << fnInvalidAddress << endl;
  cout << "reset rate  " << fnReset/dt << " Hz" << endl;
  printf("trigger rates      %8f Hz total       %2f Hz internal       %2f Hz external\n",
			fnTrig/dt, fnTrigInternal/dt, fnTrigExternal/dt);
  cout << "infinite readouts " << fnInfiniteRO << "    "  <<  float(fnInfiniteRO)/float(16*fnTrig) << " per trigger/roc" << endl;
  cout << "Readout truncated " << fnTruncated << endl;
  cout << "Bad trailer " << fnBadTrailer << endl;
  cout << "Raw efficicency (readout with hits / all ) " 
		 << ((double)fnDataWithHits)/((double)fnData) << endl;
  if(fnCalInjectHistogrammed>0){
    //	 cout << "CalInject efficiencies "<<  endl 
    //			<< "roc col row   raw efficiency"<<endl;
	 cout << "CalInject efficiencies "<<  endl 
			<< "       col row  raw efficiency"<<endl;
	 printHighEffPixels(fnCalInjectHistogrammed, 0.3, hRocMapCal);
  }
  cout << "----------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
void BinaryFileReader::printRunSummary2() { 
  float dt=(fTmax-fTmin)*25e-9;
  cout << "@@@time " << dt << " sec" << endl;
  cout << "@@@reset " << fnReset/dt << " 1/sec" << endl;
  cout << "@@@trigger " << fnTrig/dt << " 1/sec" << endl;
  
}

// ----------------------------------------------------------------------
void BinaryFileReader::printPixel(int col, int row) { 
  cout  << col << " " << row;
  for(int roc=0; roc<fNROC; roc++){
    double nhit=hRocMapCal[roc]->GetBinContent(col+1, row+1);
    printf(" %5f",nhit/double(fnCalInjectHistogrammed));
  } 
  cout << endl;
}









