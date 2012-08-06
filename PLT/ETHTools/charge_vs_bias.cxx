#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string.h>
#include <TSystem.h>
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TKey.h"
#include "TStyle.h"

using namespace std;

struct sample{
  char name[50];
  double  phi;  // listOfSamples.tab
  char particle[10];
  int nRuns;
  int trimVcal;
  double temp;
  vector<int> runNr; // rangeXYZ.txt
  vector<int> bias;
  vector<double> Ileak;
  char structure[10];
  vector<double> Q; // charge  from sprectra_XYZ.root
  vector<double> sQ;// error
};


int main(int argc, char **argv){
  
  bool convert = false, createSpectra = false, copyFiles = false, noise = false;
  
  // -- command line arguments
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-c")) {
      convert=true;
      cout <<"Convert to tree\n";
    }else if (!strcmp(argv[i],"-s")) {
      createSpectra=true;
      cout <<"Create Spectra\n";
    }else if (!strcmp(argv[i],"-f")) {
      copyFiles=true;
      cout <<"Copy config.dat and pixelMask.dat\n";
    }else if (!strcmp(argv[i],"-n")) {
      noise = true;
      cout << "Analyze noise\n";
    }else if (!strcmp(argv[i],"-help")) {
      cout <<"./charge_vs_bias -s -c -f\n";
      cout <<"  -c: convert to tree\n";
      cout <<"  -s: create spectra\n";
      cout <<"  -f: copy config.dat and pixelMask.dat\n";
      cout <<"  -n: analyze noise\n";
      cout <<"  -help: only show this message\n";
      return 0;
    }
  }
  
  
  vector<sample> s;
  sample currentSample;
  char currentName[50], currentParticle[10], currentStructure[10];
  double currentPhi, currentTemp;
  int currentTrim;
  
  // Open the file in which the sample name, fluence, particles, trim Vcal, temp, and structure are stored
  ifstream fileInSample;
  fileInSample.open("listOfSamples_joaquin.tab");
  if (fileInSample.bad()){
    cout << "!!!!!!!!!  ----> Could not open file -- terminate\n";
    return 0;
  }  
  cout << "Reading sample names from file "<<endl;
  while(fileInSample.good()){
    fileInSample>>currentName>>currentPhi>>currentParticle>>currentTrim>>currentTemp>>currentStructure;  
    strcpy(currentSample.name, currentName);
    currentSample.phi = currentPhi;
    strcpy(currentSample.particle, currentParticle);
    currentSample.trimVcal = currentTrim;
    currentSample.temp = currentTemp;
    strcpy(currentSample.structure, currentStructure);    
    //cout << currentSample.name<<" "<<currentSample.phi<<" "<<currentSample.particle<<" "<<currentSample.trimVcal<<endl;
    s.push_back(currentSample);
  }
  fileInSample.close();
  cout << "Number of Samples "<< s.size() << endl;
  
  // for Q vs V
  TGraphErrors *g[s.size()];
  double vBias[200], vSbias[200]={0.};
  double vQ[200],vSq[200];
  
  // for Ileak vs V
  TGraph *f[s.size()];
  double vBiasIV[200], vIleak[200], vIleak_corr[200];
  
  
  // for Av. Noise vs. V
  TMultiGraph *h = new TMultiGraph();    
  TLegend legNoise(.7,.7,.95,.95);
  int dot1counter=0, gap30counter=0, gap30_2counter=0, dot1counter2=0, gap30counter2=0, gap30_2counter2=0;
  bool firstLineNoise[7]={true,true,true,true,true,true,true};
  
  // for # Noisy Pix vs. V
  TMultiGraph *npix = new TMultiGraph();
  
  
  // for Q vs Phi numbers of samples by hand
  double pSig0[50], pPhi0[50];
  double pSigPi[50], pPhiPi[50];
  double pSigP600[50], pPhiP600[50];
  double pSigP800[50], pPhiP800[50];  
  double pSigP1000[50], pPhiP1000[50];  
  
  int nSig0 = 0, nSigPi =0, nSigP600=0, nSigP800=0, nSigP1000=0;
  
  // loop over all samples
  
  
  
  for (unsigned int sampleCounter = 0; sampleCounter < s.size(); sampleCounter++){
    char rangeFileName[200];
    //sprintf(rangeFileName,"/data/disk1/rohe/source_test_2008/rawdata/runlists/range_%s.txt",
    //        s.at(sampleCounter).name);
    //    sprintf(rangeFileName,"/home/l_tester/sensorTest/output/CCEStudy/singleROC_%s/range_%s.txt",
    sprintf(rangeFileName,"/home/l_tester/slc5/Output/Efficiency_meaurements/ROC_%s/range_sample.txt",
            s.at(sampleCounter).name);	    
    
    //cout << rangeFileName<<endl;	     
    
    int currentBias, currentRun;
    double currentI;
    
    // open this file and read
    ifstream fileInRange;
    fileInRange.open(rangeFileName);
    if (fileInRange.bad()){
      cout << "!!!!!!!!!  ----> Could not open file "<<rangeFileName<<endl;;
      continue;
    }      
    while(fileInRange.good()){
      fileInRange>>currentBias>>currentRun>>currentI;
      cout << currentBias<<" "<<currentRun<<endl;
      s.at(sampleCounter).bias.push_back(currentBias);
      s.at(sampleCounter).runNr.push_back(currentRun);
      s.at(sampleCounter).Ileak.push_back(currentI);     
    }
    fileInRange.close();
    s.at(sampleCounter).nRuns = s.at(sampleCounter).bias.size();
    //cout << "sample "<< s.at(sampleCounter).name;
    //cout << " number of runs "<< s.at(sampleCounter).nRuns << endl;
    
    // now process all the files
    for (int runCounter =0; runCounter < s.at(sampleCounter).nRuns; runCounter++){
      // now process all the files
      char command1[200], command2[200], command3[200], command4[200];
      
      sprintf(command1,"cp config.dat /home/l_tester/log/bt05r%06d/", s.at(sampleCounter).runNr.at(runCounter));
      if (copyFiles) system(command1);
      
      sprintf(command2,"cp pixelMask.dat /home/l_tester/log/bt05r%06d/", s.at(sampleCounter).runNr.at(runCounter));
      if (copyFiles) system(command2);
      
      //int trimVal = 60;
      //if (s.at(sampleCounter).phi > 20) trimVal = 30;
      
      sprintf(command3,"./convert_to_tree -l -r %d -t %d", s.at(sampleCounter).runNr.at(runCounter), s.at(sampleCounter).trimVcal);
      if (convert){cout << command3 << endl; system(command3);}
      
      sprintf(command4,"./pulseHeightMain %d", s.at(sampleCounter).runNr.at(runCounter));
      if (createSpectra){cout << command4 << endl; system(command4);}
      
      // now open the "spectra file and get out the pulse height
      TFile *rootFile;
      char rootFileName[200];
      sprintf(rootFileName,
              "/home/l_tester/log/bt05r%06d/spectra_%06d.root",
	      s.at(sampleCounter).runNr.at(runCounter),s.at(sampleCounter).runNr.at(runCounter));
      rootFile = new TFile(rootFileName,"read");
      char parHisName[50], errHisName[50];
      TH2D *parHis, *errHis;
      sprintf(parHisName,"hFitPar_%d", s.at(sampleCounter).runNr.at(runCounter));
      sprintf(errHisName,"hFitErr_%d", s.at(sampleCounter).runNr.at(runCounter));     
      
      parHis = (TH2D *) rootFile->Get(parHisName);
      errHis = (TH2D *) rootFile->Get(errHisName);      
      
      s.at(sampleCounter).Q.push_back(parHis->GetBinContent(2,2));
      s.at(sampleCounter).sQ.push_back(errHis->GetBinContent(2,2));
      
      rootFile->Close();
      delete(rootFile);
      
    }// end loop over all runs of a sample
    
    // create a graph (Q vs V) for the current sample
    // loop over the runs
    int countProblems = 0, biasTooLow = 0, usedRuns = 0;
    
    for (int counter=0; counter < s.at(sampleCounter).nRuns; counter++){
      
      bool use = true;
      
      if (s.at(sampleCounter).sQ.at(counter) > 0.05*s.at(sampleCounter).Q.at(counter)){
	// fit probably failed
        cout <<"-------> look at run "<< s.at(sampleCounter).runNr.at(counter)<<endl;
	countProblems++;
	use = false;
      }else if (counter>1){ 
        if (s.at(sampleCounter).Q.at(counter) < 0.8 * s.at(sampleCounter).Q.at(counter-1)){
	  // noise peak fitted?
          cout <<"-------> look at run "<< s.at(sampleCounter).runNr.at(counter)<<" signal became smaller with bias"<<endl;      
	  countProblems++;
	  use = false;	
	}
      }
      if (s.at(sampleCounter).runNr.at(counter) == 4300){
        cout <<"-------> manually decided not to use run "<< s.at(sampleCounter).runNr.at(counter)<<endl;
	countProblems++;
	use = false;      
      }
      if ((s.at(sampleCounter).bias.at(counter) < 100.)&&(s.at(sampleCounter).phi > 0)){
	cout <<"-------> not using run "<<s.at(sampleCounter).runNr.at(counter)<<" bias only "<<s.at(sampleCounter).bias.at(counter)<< endl;
	biasTooLow++;
	use = false;
      }
      if (use){ 
        vBias[usedRuns] = double(s.at(sampleCounter).bias.at(counter));
        vQ[usedRuns] = s.at(sampleCounter).Q.at(counter) * .065;// conversion Vcal -> k electrons
        vSq[usedRuns] = s.at(sampleCounter).sQ.at(counter) * .065;
	usedRuns++;
      }
      
      vBiasIV[counter] = s.at(sampleCounter).bias.at(counter);
      vIleak[counter] = s.at(sampleCounter).Ileak.at(counter);
      
      //correct Ileak for temperature
      double T1 = 263.15, T2;  // T1 = 2.63.15 K (-10C)
      double E_g = 1.12, k = .00008617;  // E_g = 1.12 eV, k = 8.617*10^-5 eV/K
      T2 = s.at(sampleCounter).temp+273.15;  // T2 = measured temp, converted from C to K
      vIleak_corr[counter] = vIleak[counter]*(pow(T1,2)/pow(T2,2))*exp((E_g/(2*k))*((1/T2)-(1/T1)));  //current adjusted for temp; see page 43 of "Pixel Detectors: From Fundamentals to Applications"
                                                                                                      // L. Rossi, P. Fischer, T. Rohe, N. Wermes
      //check that calculation is correct
      if(counter == 0){
	cout << "For the first run of the chip:" << endl;
	cout << "Measured leakage current = " << vIleak[counter] << endl;
	cout << "Corrected leakage current = " << vIleak_corr[counter] << endl;
      }
      
      // for the Q vs Phi plot
      // Phi = 0
      if ((strcmp("oo",s.at(sampleCounter).particle) ==0)&&(s.at(sampleCounter).bias.at(counter) <= 250.)){
        pSig0[nSig0]= s.at(sampleCounter).Q.at(counter) * .065;// conversion Vcal -> k electrons; 
	pPhi0[nSig0]= s.at(sampleCounter).phi;     
      }
      if ((strcmp("pi",s.at(sampleCounter).particle) ==0)&&(s.at(sampleCounter).bias.at(counter) <= 605.)){
        pSigPi[nSigPi]= s.at(sampleCounter).Q.at(counter) * .065; 
	pPhiPi[nSigPi]= s.at(sampleCounter).phi;     
      }
      if (strcmp("pr",s.at(sampleCounter).particle) ==0){
        if (s.at(sampleCounter).bias.at(counter) <= 605.){
          pSigP600[nSigP600]= s.at(sampleCounter).Q.at(counter) * .065; 
	  pPhiP600[nSigP600]= s.at(sampleCounter).phi;     
	}
	if (s.at(sampleCounter).phi >20.){
	  if(s.at(sampleCounter).bias.at(counter) <= 805.){
            pSigP800[nSigP800]= s.at(sampleCounter).Q.at(counter) * .065; 
	    pPhiP800[nSigP800]= s.at(sampleCounter).phi;     
	  } 
	  if(s.at(sampleCounter).bias.at(counter) <= 1005.){
            pSigP1000[nSigP1000]= s.at(sampleCounter).Q.at(counter) * .065; 
	    pPhiP1000[nSigP1000]= s.at(sampleCounter).phi;     
	  } 	
	}
      }
      
    }
    if ((countProblems > 0)||(biasTooLow >0)){ 
      cout<<s.at(sampleCounter).name<<" has "<<countProblems+biasTooLow<<" runs out of "<<s.at(sampleCounter).nRuns<<endl;
    }
    int length = s.at(sampleCounter).nRuns - countProblems - biasTooLow;
    g[sampleCounter] = new TGraphErrors(length, vBias, vQ, vSbias, vSq);
    f[sampleCounter] = new TGraph(s.at(sampleCounter).nRuns,vBiasIV,vIleak_corr);
    
    // for Q vs Phi plot
    if  (strcmp("oo",s.at(sampleCounter).particle) ==0) nSig0++;
    if  (strcmp("pi",s.at(sampleCounter).particle) ==0) nSigPi++;
    if  (strcmp("pr",s.at(sampleCounter).particle) ==0){
      nSigP600++;
      if (s.at(sampleCounter).phi >20.){nSigP800++; nSigP1000++;}
    }
    
    //analyze noise if indicated
    char command5[200];
    sprintf(command5,"./noise %s", s.at(sampleCounter).name);
    if (noise) {cout << command5 << endl; system(command5);}
    
    //Add noise graph to multigraph of all samples
    TGraphErrors *ngraph;
    char noiseRootFileName[200];
    sprintf(noiseRootFileName, "/home/l_tester/slc5/Output/Efficiency_meaurements/ROC_%s/Noise_scurve.root", s.at(sampleCounter).name);
    TFile *noiseRootFile = new TFile(noiseRootFileName,"read");
    TList* list = noiseRootFile->GetListOfKeys();
    for (int i = 0; i < list->GetSize(); i++)
      {
	TKey *key = (TKey *)(list->At(i));
	TObject* obj = key->ReadObj();
	if (obj->InheritsFrom("TGraphErrors"))
	  {
	    ngraph = (TGraphErrors *)obj;
	  }
      }
    TGraph* np = (TGraph*)noiseRootFile->Get( "NumOfNoisyPixels");
    
       
    //Sort noise by chip structure
    if (strcmp(s.at(sampleCounter).structure, "dot1") == 0)
      {
	ngraph->SetMarkerStyle(20);
	np->SetMarkerStyle(20);
	ngraph->SetMarkerColor(2);
	np->SetMarkerColor(2);
	if(dot1counter == 0)
	  {
	    legNoise.AddEntry(ngraph,"Dot 1","p");
	  }
	dot1counter++;
      }
    
    else if (strcmp(s.at(sampleCounter).structure, "gap30") == 0)
      {
	ngraph->SetMarkerStyle(21);
	np->SetMarkerStyle(21);
	ngraph->SetMarkerColor(3);
	np->SetMarkerColor(3);
	if(gap30counter == 0)
	  {
	    legNoise.AddEntry(ngraph,"Gap 30","p");
	  }
	gap30counter++;
      }
    
    else if (strcmp(s.at(sampleCounter).structure, "gap30-2") == 0)
      {
	ngraph->SetMarkerStyle(22);
	np->SetMarkerStyle(22);
	ngraph->SetMarkerColor(4);
	np->SetMarkerColor(4);
	if(gap30_2counter == 0)
	  {
	    legNoise.AddEntry(ngraph,"Gap 30-2","p");
	  }
	gap30_2counter++;
      }
    assert(0);
 
    
    // Sort chips by fluence
    if ((s.at(sampleCounter).phi < .1)&&(s.at(sampleCounter).phi > -.1)) 
      {
	ngraph->SetLineColor(12);
	ngraph->SetLineWidth(2);
	ngraph->SetLineStyle(1);
	np->SetLineColor(12);
	np->SetLineWidth(2);
	np->SetLineStyle(1);
	if (firstLineNoise[0]){firstLineNoise[0]=false;legNoise.AddEntry(ngraph,"Not irradiated","l");} 
      }
    else if (((s.at(sampleCounter).phi - 4.3) < .05)&&((s.at(sampleCounter).phi - 4.3) > -.05))
      {	  
	ngraph->SetLineColor(42); 
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(1);
	np->SetLineColor(42); 
	np->SetLineWidth(1);
	np->SetLineStyle(1);
	if (firstLineNoise[1]){firstLineNoise[1]=false;legNoise.AddEntry(ngraph,"4.3E14, Pions","l");}       
      }  
    else if (((s.at(sampleCounter).phi - 6.2) < .05)&&((s.at(sampleCounter).phi - 6.2) > -.05)) 
      {	 
	ngraph->SetLineColor(30);
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(2);
	np->SetLineColor(30);
	np->SetLineWidth(1);
	np->SetLineStyle(2);
	if (firstLineNoise[2]){firstLineNoise[2]=false;legNoise.AddEntry(ngraph,"6.2E14, Pions","l");}       
      }  
    else if (((s.at(sampleCounter).phi - 6.1) < .05)&&((s.at(sampleCounter).phi - 6.1) > -.05)) 
      {	 
	ngraph->SetLineColor(41);
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(1);
	np->SetLineColor(41);
	np->SetLineWidth(1);
	np->SetLineStyle(1);
	if (firstLineNoise[3]){firstLineNoise[3]=false;legNoise.AddEntry(ngraph,"6.1E14, Protons","l");}
      }      
    else if (((s.at(sampleCounter).phi - 11.0) < .05)&&((s.at(sampleCounter).phi - 11.0) > -.05)) 
      {
	ngraph->SetLineColor(40);
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(2);
	np->SetLineColor(40);
	np->SetLineWidth(1);
	np->SetLineStyle(2);
	if (firstLineNoise[4]){firstLineNoise[4]=false;legNoise.AddEntry(ngraph,"1.1E15, Protons","l");}
      }   
    else if (((s.at(sampleCounter).phi - 28.0) < .05)&&((s.at(sampleCounter).phi - 28.0) > -.05)) 
      {
	ngraph->SetLineColor(50); 
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(4);
	np->SetLineColor(46); 
	np->SetLineWidth(1);
	np->SetLineStyle(4);
	if (firstLineNoise[5]){firstLineNoise[5]=false;legNoise.AddEntry(ngraph,"2.8E15, Protons","l");}
      }
    else if (((s.at(sampleCounter).phi - 51.0) < .05)&&((s.at(sampleCounter).phi - 51.0) > -.05)) 
      {
	ngraph->SetLineColor(46); 
	ngraph->SetLineWidth(1);
	ngraph->SetLineStyle(9);
	np->SetLineColor(46); 
	np->SetLineWidth(1);
	np->SetLineStyle(9);
	if (firstLineNoise[6]){firstLineNoise[6]=false;legNoise.AddEntry(ngraph,"5.1E15, Protons","l");}
      }
    ngraph->SetMarkerSize(0.8);
    np->SetMarkerSize(0.8);
    h->Add(ngraph);
    npix->Add(np);
    
    
  }// end loop over all samples


  TFile all("AllChipPlots.root", "RECREATE");

  
  // plot Graphs charge vs bias
  TCanvas c1;
  c1.cd();
  c1.SetFillColor(10);
  TH2D his("his","",10,0,1150,10,0,27);
  his.SetStats(0);
  his.GetXaxis()->SetTitle("Voltage [V]");
  his.GetYaxis()->SetTitle("Signal [ke^{-}]");  
  his.Draw();
  
  TLegend leg(.6,.55,.85,.8);
  leg.SetFillColor(0);
  
  bool firstLine[7]={true,true,true,true,true,true,true};
  for (int sampleCounter=0; sampleCounter< s.size(); sampleCounter++){
    
    int markerType = 0, lineColor = 1, lineWidth = 1, lineStyle = 1;
    
    if ((s.at(sampleCounter).phi < .1)&&(s.at(sampleCounter).phi > -.1)) {
      markerType = 20; lineColor = 1; lineWidth=2, lineStyle=1;
      if (firstLine[0]){firstLine[0]=false;leg.AddEntry(g[sampleCounter],"Not irradiated","l");} 
    }
    if (((s.at(sampleCounter).phi - 4.3) < .05)&&((s.at(sampleCounter).phi - 4.3) > -.05)){
      markerType = 24; lineColor = 1; lineWidth=1, lineStyle=1;
      if (firstLine[1]){firstLine[1]=false;leg.AddEntry(g[sampleCounter],"4.3E14, Pions","l");}       
    }  
    if (((s.at(sampleCounter).phi - 6.2) < .05)&&((s.at(sampleCounter).phi - 6.2) > -.05)) {
      markerType = 25; lineColor = 1; lineWidth=1, lineStyle=2;
      if (firstLine[2]){firstLine[2]=false;leg.AddEntry(g[sampleCounter],"6.2E14, Pions","l");}             
    }  
    if (((s.at(sampleCounter).phi - 6.1) < .05)&&((s.at(sampleCounter).phi - 6.1) > -.05)) {
      markerType = 21; lineColor = 14; lineWidth=1, lineStyle=1;
      if (firstLine[3]){firstLine[3]=false;leg.AddEntry(g[sampleCounter],"6.1E14, Protons","l");}             
    }      
    if (((s.at(sampleCounter).phi - 11.0) < .05)&&((s.at(sampleCounter).phi - 11.0) > -.05)) {
      markerType = 22; lineColor = 14;lineWidth=1, lineStyle=2;
      if (firstLine[4]){firstLine[4]=false;leg.AddEntry(g[sampleCounter],"1.1E15, Protons","l");}                   
    }   
    if (((s.at(sampleCounter).phi - 28.0) < .05)&&((s.at(sampleCounter).phi - 28.0) > -.05)) {
      markerType = 28;lineColor = 1; lineWidth=1, lineStyle=4;
      if (firstLine[5]){firstLine[5]=false;leg.AddEntry(g[sampleCounter],"2.8E15, Protons","l");}                         
    } 
    if (((s.at(sampleCounter).phi - 51.0) < .05)&&((s.at(sampleCounter).phi - 51.0) > -.05)) {
      markerType = 26;lineColor = 1; lineWidth=1, lineStyle=9;
      if (firstLine[6]){firstLine[6]=false;leg.AddEntry(g[sampleCounter],"5.1E15, Protons","l");}                         
    } 

    //modifications for comparing several of the same chip between the two cold boxes
    //if (sampleCounter < 3){ markerType = 20;}
    //else { markerType = 24;}
    
    
    g[sampleCounter]->SetMarkerStyle(markerType);
    g[sampleCounter]->SetLineColor(lineColor);
    g[sampleCounter]->SetLineWidth(lineWidth);       
    g[sampleCounter]->SetLineStyle(lineStyle);    
    g[sampleCounter]->Draw("plsame");
    
  }
  leg.Draw();
  c1.Print("charge_vs_bias.eps");  
  c1.Print("charge_vs_bias.pdf");
  //c1.SetName("Charge_vs_Bias");
  c1.Write("Charge_vs_Bias");
  //c1.Print("charge_vs_bias.tif");
  
  // make Q vs Phi -plot
  TCanvas c2;
  c2.cd();
  c2.SetFillColor(10);
  TH2D his2("his","",10,-1,30,10,0,27);
  his2.SetStats(0);
  his2.GetXaxis()->SetTitle("#Phi [10^{14} n_{eq}/cm^{2}]");
  his2.GetYaxis()->SetTitle("Signal [ke^{-}]");  
  his2.Draw();  
  
  TLegend leg2(.6,.6,.85,.85);
  leg2.SetFillColor(0);
  
  TGraph *p_zero;
  p_zero = new TGraph(nSig0,pPhi0,pSig0);
  p_zero->SetMarkerStyle(20);
  p_zero->Draw("psame");
  leg2.AddEntry(p_zero,"Not irradiated, 250V","p");
  
  TGraph *p_pi;
  p_pi = new TGraph(nSigPi,pPhiPi,pSigPi);
  p_pi->SetMarkerStyle(22);
  p_pi->Draw("psame");  
  leg2.AddEntry(p_pi,"Pions, 600V","p");  
  
  TGraph *p_p600;
  p_p600 = new TGraph(nSigP600,pPhiP600,pSigP600);
  p_p600->SetMarkerStyle(24);
  p_p600->Draw("psame");    
  leg2.AddEntry(p_p600,"Protons, 600V","p");    
  
  
  TGraph *p_p800;
  p_p800 = new TGraph(nSigP800,pPhiP800,pSigP800);
  p_p800->SetMarkerStyle(26);
  p_p800->Draw("psame");
  leg2.AddEntry(p_p800,"Protons, 800V","p");         
  
  TGraph *p_p1000;
  p_p1000 = new TGraph(nSigP1000,pPhiP1000,pSigP1000);
  p_p1000->SetMarkerStyle(28);
  p_p1000->Draw("psame");     
  leg2.AddEntry(p_p1000,"Protons, 1000V","p");        
  
  leg2.Draw();
  
  c2.Print("charge_vs_phi.eps");
  c2.Print("charge_vs_phi.pdf");
  c2.Write("Charge_vs_Phi");  
  //c2.Print("charge_vs_phi.tif"); 
  
  
  // plot Graphs ileak vs bias
  TCanvas c3;
  c3.cd();
  c3.SetFillColor(10);
  c3.SetLogy(true);
  TH2D hisIV("hisIV","",10,0,1200,10,0.1,100);
  hisIV.SetStats(0);
  hisIV.GetXaxis()->SetTitle("Voltage [V]");
  hisIV.GetYaxis()->SetTitle("Leakage Current [micro A]");  
  hisIV.Draw();
  
  TLegend legIV(.6,.2,.85,.55);
  legIV.SetFillColor(10);
  
  bool firstLineIV[7]={true,true,true,true,true,true,true};
  for (int sampleCounter=0; sampleCounter< s.size(); sampleCounter++){
    
    int markerType = 0, lineColor = 1, lineWidth = 1, lineStyle = 1;
    
    if ((s.at(sampleCounter).phi < .1)&&(s.at(sampleCounter).phi > -.1)) {
      lineColor = 1; lineWidth=2, lineStyle=1;
      if (firstLineIV[0]){firstLineIV[0]=false;legIV.AddEntry(f[sampleCounter],"Not irradiated","l");} 
    }
    else if (((s.at(sampleCounter).phi - 4.3) < .05)&&((s.at(sampleCounter).phi - 4.3) > -.05)){
      lineColor = 1; lineWidth=1, lineStyle=1;
      if (firstLineIV[1]){firstLineIV[1]=false;legIV.AddEntry(f[sampleCounter],"4.3E14, Pions","l");}       
    }  
    else if (((s.at(sampleCounter).phi - 6.2) < .05)&&((s.at(sampleCounter).phi - 6.2) > -.05)) {
      lineColor = 1; lineWidth=1, lineStyle=2;
      if (firstLineIV[2]){firstLineIV[2]=false;legIV.AddEntry(f[sampleCounter],"6.2E14, Pions","l");}             
    }  
    else if (((s.at(sampleCounter).phi - 6.1) < .05)&&((s.at(sampleCounter).phi - 6.1) > -.05)) {
      lineColor = 14; lineWidth=1, lineStyle=1;
      if (firstLineIV[3]){firstLineIV[3]=false;legIV.AddEntry(f[sampleCounter],"6.1E14, Protons","l");}        
    }      
    else if (((s.at(sampleCounter).phi - 11.0) < .05)&&((s.at(sampleCounter).phi - 11.0) > -.05)) {
      lineColor = 14;lineWidth=1, lineStyle=2;
      if (firstLineIV[4]){firstLineIV[4]=false;legIV.AddEntry(f[sampleCounter],"1.1E15, Protons","l");}           
    }   
    else if (((s.at(sampleCounter).phi - 28.0) < .05)&&((s.at(sampleCounter).phi - 28.0) > -.05)) {
      lineColor = 1; lineWidth=1, lineStyle=4;
      if (firstLineIV[5]){firstLineIV[5]=false;legIV.AddEntry(f[sampleCounter],"2.8E15, Protons","l");}           
    } 
    else if (((s.at(sampleCounter).phi - 51.0) < .05)&&((s.at(sampleCounter).phi - 51.0) > -.05)) {
      lineColor = 1; lineWidth=1, lineStyle=9;
      if (firstLineIV[6]){firstLineIV[6]=false;legIV.AddEntry(f[sampleCounter],"5.1E15, Protons","l");}           
    } 


    
    //Sort chips by structure
    if (strcmp(s.at(sampleCounter).structure, "dot1") == 0)
      {
	markerType = 20;
	f[sampleCounter]->SetMarkerColor(2);
	if(dot1counter2 == 0)
	  {
	    legIV.AddEntry(f[sampleCounter],"Dot 1","p");
	  }
	dot1counter2++;
      }
    
    else if (strcmp(s.at(sampleCounter).structure, "gap30") == 0)
      {
	markerType = 21;
	f[sampleCounter]->SetMarkerColor(3);
	if(gap30counter2 == 0)
	  {
	    legIV.AddEntry(f[sampleCounter],"Gap 30","p");
	  }
	gap30counter2++;
      }
    
    else if (strcmp(s.at(sampleCounter).structure, "gap30-2") == 0)
      {
	markerType = 22;
	f[sampleCounter]->SetMarkerColor(4);
	if(gap30_2counter2 == 0)
	  {
	    legIV.AddEntry(f[sampleCounter],"Gap 30-2","p");
	  }
	gap30_2counter2++;
      }
    
    f[sampleCounter]->SetMarkerStyle(markerType);
    f[sampleCounter]->SetMarkerSize(0.8);
    f[sampleCounter]->SetLineColor(lineColor);
    f[sampleCounter]->SetLineWidth(lineWidth);       
    f[sampleCounter]->SetLineStyle(lineStyle);    
    f[sampleCounter]->Draw("lpsame");
    
  }
  legIV.Draw();
  c3.Print("ileak_vs_bias.eps");  
  c3.Print("ileak_vs_bias.pdf");
  c3.Write("Ileak_vs_Bias");
  
  
  // plot Graphs av.noise vs bias
  TCanvas c4;
  c4.cd();
  c4.SetFillColor(10);
  h->SetTitle("Av. Noise vs. Vbias");
  h->Draw("ALP");
  h->GetXaxis()->SetTitle("Voltage [V]");
  h->GetYaxis()->SetTitle("Average Noise [e]");
  h->Draw("ALP");
  legNoise.Draw();
  c4.Print("noise_vs_bias.eps");
  c4.Print("noise_vs_bias.pdf");
  c4.Write("Noise_vs_Bias");

  // plot Graphs # noisey pixels vs bias
  TCanvas c5;
  c5.cd();
  c5.SetFillColor(10);
  npix->SetTitle("Number of Noisy Pixels vs. Vbias");
  npix->Draw("ALP");
  npix->GetXaxis()->SetTitle("Voltage [V]");
  npix->GetYaxis()->SetTitle("#");
  npix->Draw("ALP");
  legNoise.Draw();
  c5.Print("noisyPix_vs_bias.eps");
  c5.Print("noisyPix_vs_bias.pdf");
  c5.Write("NoisyPix_vs_Bias");
  
  
  all.Close();

}

