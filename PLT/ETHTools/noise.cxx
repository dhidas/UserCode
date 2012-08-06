#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TLegend.h"
#include "TKey.h"
#include "TList.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "iostream"
#include <fstream>
#include "TGraph2D.h"
#include "fstream"
#include <string.h>
#include<vector>
#include<sstream>
#include "TPaveText.h"
#include <stdio.h>

using namespace std;

int main (int argc, char* argv[50]) 
{

   char name[50];

   if(argc < 2) {

        // in case non suffient commandline input was given
        cout<<"Usage: ./noise dirName \n";
        return 0;
   }
   strcpy(name, argv[1]);

   int Vbias[35];
   int Run_name[35];
   double Ileak[35];


   ifstream read_in;
   char rangeFile [200];
    sprintf(rangeFile,"/home/l_tester/slc5/Output/Efficiency_measurements/ROC_%s/range_%s.txt",name, name);
   read_in.open(rangeFile);
   int cnt_run=0;

   cout << "Reading file " << rangeFile << endl;
   while(!read_in.eof() && read_in >> Vbias[cnt_run] >> Run_name[cnt_run] >> Ileak[cnt_run])
   {
     cout << "Bias  " << Vbias[cnt_run] <<" Run "<< Run_name[cnt_run] << "   I leak  " << Ileak[cnt_run] << endl;
     cnt_run++;
//      cout<<"la marimonda"<<endl;
   }


   read_in.close();


   char rootFileName[200];
   sprintf(rootFileName, "/home/l_tester/slc5/Output/Efficiency_measurements/ROC_%s/Noise_scurve.root", name);
   TFile rootFile(rootFileName, "RECREATE");
   rootFile.Close();


   double AvNoise[35];
   double AvErr[35];
   int noisy[35];


   // Loop over runs
   for( int b = 0; b < cnt_run; b++)
     {
       TFile rootFile2(rootFileName, "UPDATE");
       
       TCanvas* cn = new TCanvas();
       cn->Divide(1,2);

       TH2F * noiseMap = new TH2F(Form("Noise Map %iV", Vbias[b]), "Noise Map;Row;Column",80,0,80,52,0,52);
       TH1F * NoiseDist = new TH1F(Form("Noise Dis.%iV", Vbias[b]), "Noise;Noise (electrons);Frequency",500,0,500);
       TH1F * VcalNoiseDist = new TH1F(Form("Vcal Noise Dis.%iV", Vbias[b]), "Noise;Noise (Vcal);Frequency",100,0,5);


       char Path[200];
       sprintf(Path, "/home/l_tester/slc5/Output/Efficiency_measurements/ROC_%s/SCurve_C0_Vbias%i.dat", name,Vbias[b]);
       cout << Path << endl;
       ifstream input;
       input.open(Path);
     
       if(!input)
	 {
	   cout << "Cannot find file" << endl;
	 }
       std::string str;
       double Threshold, sigma;
       int row_n, col_n;
      
       while(!input.eof() && input >> str)
	 {
	   if(str == "Sigma")
	     {
	       break;
	     }
	 }
       while(!input.eof() && input >> Threshold >> sigma >> str >> col_n >> row_n)
	 {
	   noiseMap->Fill(row_n,col_n,sigma);
	   NoiseDist->Fill(sigma);
	   VcalNoiseDist->Fill(sigma/65.0);
	 }
       input.close();

       noiseMap->SetOption("surf2");
       cn->cd(1);
       noiseMap->Draw();
       noiseMap->Write();
       cn->cd(2);
       NoiseDist->Draw();
       NoiseDist->Fit("gaus");
       NoiseDist->Draw();
       NoiseDist->Write();
       TF1 *fit = NoiseDist->GetFunction("gaus");
       AvNoise[b] = fit->GetParameter(1);
	 cout<<AvNoise[b]<<" ";
       AvErr[b] = fit->GetParameter(2);
       VcalNoiseDist->Draw();
       VcalNoiseDist->Write();
       rootFile2.Close();

       cout << "Run = " << Run_name[b] << endl;
       // Number of noisy pixels, from spectra.root file
       char spectraFileName[200], maskHistName[200];

	 if (Run_name[b]<1000)sprintf(spectraFileName, "/home/l_tester/log/bt05r000%i/spectra_000%i.root", Run_name[b], Run_name[b]);
	 else sprintf(spectraFileName, "/home/l_tester/log/bt05r00%i/spectra_00%i.root", Run_name[b], Run_name[b]);
       sprintf(maskHistName, "pMask_%i", Run_name[b]);
       TFile *spectraFile = new TFile(spectraFileName,"read");
       TH2* mask = (TH2*)spectraFile->Get(maskHistName);
       int noisyPixels = 0, errCode;

       // loop over pixels
       for (int row = 0; row < 80; row++)
	 {
	   for (int col = 0; col < 52; col++)
	     {
	       errCode = (int)mask->GetBinContent(row,col);
	       errCode /= 2;
	       errCode /= 2;
	       if (errCode % 2 == 1)
		 {
		   noisyPixels++;
		 }
	     }
	 }
       spectraFile->Close();
       noisy[b] = noisyPixels;
     }

   TFile rootFile3(rootFileName, "UPDATE");

   /*for(int j = 0; j < cnt_run; j++)
     {
       cout << "Vbias = " << Vbias[j] << "\t# Noisy = " << noisy[j] << endl;
       }*/

   TGraph * noisyPixelDist = new TGraph(cnt_run,Vbias,noisy); 
   noisyPixelDist->SetTitle("Number of Noisy Pixels");
   noisyPixelDist->Draw("AC*");
   noisyPixelDist->GetXaxis()->SetTitle("Vbias [V]");
   noisyPixelDist->GetYaxis()->SetTitle("#");
   noisyPixelDist->Draw("AC*");
   noisyPixelDist->SetName("NumOfNoisyPixels");
   noisyPixelDist->Write();

   double Vbiasd[35];
   int dbl;
   for (int k=0; k < cnt_run; k++)
     {
       dbl = Vbias[k];
       Vbiasd[k] =(double)dbl;
     }

   TGraphErrors* gr = new TGraphErrors(cnt_run,Vbiasd,AvNoise,0,AvErr);
    gr->SetTitle("Average Noise vs Vbias");
    gr->SetName("Noise_vs_Vbias");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(4);
    TCanvas* c1n = new TCanvas();
    c1n->cd();
    gr->SetMinimum(0.0);
    gr->SetMaximum(300.0);
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("Vbias (V)");
    gr->GetYaxis()->SetTitle("Av. Noise (e^{- })");
    double mean = gr->GetMean(2), stdev = gr->GetRMS(2);
    double rms = sqrt(pow(mean,2.0)+pow(stdev,2.0));
    char line1[60], line2[60], line3[60];
    sprintf(line1,"Mean =  %f", mean);
    sprintf(line2,"RMS = %f", rms);
    sprintf(line3,"Standard Deviation = %f", stdev);
    TPaveText *p1 = new TPaveText(300.0,15000.0,500.0,19000.0);
    TText *t1;
    t1= p1->AddText(line1);
    TText *t2;
    t2= p1->AddText(line2);
    TText *t3;
    t3= p1->AddText(line3);
    p1->Draw();
    cout << "Mean of Noise = " << mean << endl;
    cout << "RMS of Noise = " << rms << endl;
    cout << "Standard deviation of Noise = " << stdev << endl;
    gr->Write();
    gr->Draw("AP*");
    c1n->Write();

    char epsFile[200];
    sprintf(epsFile, "/home/l_tester/slc5/Output/Efficiency_measurements/ROC_%s/AvNoiseVsVbias.root", name);
    c1n->Print(epsFile);
    rootFile3.Close();
}
