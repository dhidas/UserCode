#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <TSystem.h>
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TKey.h"

using namespace std;

struct chip{
  char name[50];
  char structure[10];
  int Vbias;
};

int main(){

  vector<chip> s;
  chip currentSample;
  char currentName[50], currentStructure[10];
  int currentVbias;
  
  bool dot1 = true, gap30 = true, gap30_2 = true, unknown = true;

  // Open the file in which the sample name and structure are stored
  ifstream fileInSample;
  fileInSample.open("noiseCompare.tab");
  if (fileInSample.bad()){
    cout << "!!!!!!!!!  ----> Could not open file noiseCompare.tab -- terminate\n";
    return 0;
  }  
  cout << "Reading sample names from file "<<endl;
  while(fileInSample.good()){
    fileInSample>>currentName>>currentStructure>>currentVbias;  
    strcpy(currentSample.name, currentName);
    strcpy(currentSample.structure, currentStructure);    
    currentSample.Vbias = currentVbias;
    s.push_back(currentSample);
  }
  fileInSample.close();
  cout << "Number of Samples "<< s.size() << endl;


  char rootFileName[200];
  sprintf(rootFileName, "noiseCompare.root");
  TFile *rootFile = new TFile(rootFileName,"RECREATE");

  TCanvas cn("AllHists","All Histograms",1);
  TLegend *leg = new TLegend(0.70,0.70,0.85,0.85);

  for (int i = 0; i < s.size(); i++) {

    char noiseRootFileName[200], histName[200], histName2[200];
    sprintf(histName, "Vcal Noise Dis.%iV", s.at(i).Vbias);
    sprintf(histName2, "%s at %iV", s.at(i).name, s.at(i).Vbias);
    TH1F *noiseHist = new TH1F(histName2, "Noise;Noise (Vcal);Frequency",100,0,5);
    sprintf(noiseRootFileName, "/home/l_tester/sensorTest/output/CCEStudy/singleROC_%s/Noise_scurve.root", s.at(i).name);
    TFile *noiseRootFile = new TFile(noiseRootFileName,"READ");
    if (!noiseRootFile->IsOpen()){
      cout << "Could not open file " << noiseRootFileName << endl;
      break;
    }
    cout << "Reading file " << noiseRootFileName << endl;

    TH1F* h = (TH1F*)noiseRootFile->Get(histName);
    double binContent;
    for (int bin = 0; bin < 100; bin++){
      binContent = h->GetBinContent(bin);
      noiseHist->SetBinContent(bin, binContent);
    }
    
    cout << "got histograms" << endl;
    noiseRootFile->Close();
    rootFile->cd();

    //Draw histogram
    if (strcmp(s.at(i).structure,"dot1")==0){ 
      noiseHist->SetLineStyle(1);
      if (dot1){
	dot1 = false;
	leg->AddEntry(noiseHist, "Dot 1", "l");
      }
    }
    else if (strcmp(s.at(i).structure,"gap30")==0){
      noiseHist->SetLineStyle(2);
      if (gap30){
	gap30 = false;
	leg->AddEntry(noiseHist, "Gap 30", "l");
      }
    }
    else if (strcmp(s.at(i).structure,"gap30-2")==0){ 
      noiseHist->SetLineStyle(6); 
      if (gap30_2){
	gap30_2 = false;
	leg->AddEntry(noiseHist, "Gap 30-2", "l");
      }
    }
    else {
      cout << "Warning: " << s.at(i).name << " Structure not recognized"; 
      noiseHist->SetLineStyle(3); 
      if(unknown){
	unknown = false;
	leg->AddEntry(noiseHist, "Unknown structure", "l");
      }
    }

    noiseHist->SetStats(kFALSE);
    noiseHist->Write();
    cn.cd();
    if (i == 0){ noiseHist->SetMaximum(450); noiseHist->Draw();}
    else{ noiseHist->Draw("same");}
 
  }

  leg->Draw("same");
  cn.Write();
  rootFile->Close();

  cout << "Finished!" << endl;

}
