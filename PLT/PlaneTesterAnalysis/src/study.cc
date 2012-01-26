
#define study_cxx
#include "study.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include "measurer.C"
#include "TH1F.h"
#include "TH3F.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <cstdlib>
#include "TMath.h"
#include "TAxis.h"
#include "dead_noisy.h"//see if that does it, if not, append the file to the end of this one.
//#include<stdio.h>
//using namespace std;
bool debug = false;
char * input_file;
char * output_file;
char * output_file_coeficients;
char * calibration_coeficients_file;

int runid;
int rocid;
int calibration_run;

TH1F *pix_e;

double caldata[100][100][4];  //calibration data, obtained from .txt file, which is created by vcaltest.C
int ipixel[100][100];// map of [ic,ir]->int index;
//caldata[column][row][0] is the 0th parameter for pol3 fitting and so on...

//  Calibration containers
const int max_events = 10000;
const int ic_max = 100;
const int ir_max = 100;

	int min_col = 11;
	int max_col = 41;
	int n_col = 30;
	
	int min_row = 39;
	int max_row = 82;//make sure these are all up to date. 
	int n_row =  43;

int _adc[ic_max][ir_max][max_events];
int _vcal[ic_max][ir_max][max_events];
int _range[ic_max][ir_max][max_events];
int _ievent[ic_max][ir_max];


double a0[ic_max][ir_max];
double a1[ic_max][ir_max];
double a2[ic_max][ir_max];
double a3[ic_max][ir_max];
double sus[ic_max][ir_max];
int noisy_N[ic_max][ir_max];//changed to int from double
int dead_N[ic_max][ir_max];//changed to int from double
float noisy_E[ic_max][ir_max];
float dead_E[ic_max][ir_max];

int saturated_vcal;
int minimal_vcal;



int study::run(int argc, char** argv ){

	//opens the calibration txt file created by calibrator.C
	if ((argc!=5 && argc!=4) || argc<4){
		printf("\n\n\n Warning! Need two arguments:\n");
		printf("./study analyze $rocid $run_number $calibration_run_number\n");
		printf("./study calibrate $rocid $calibration_run_number\n");
		printf("./study scurve $rocid $scurve_run_number\n");
		printf("Exiting\n\n\n");
		printf("argc = %d\n",argc);
		return 1;
	}
	printf("%s %s %s\n",argv[1],argv[2], argv[3]);
	char * mode = argv[1];
	TFile *f;
	FILE *fcoefs;
	TTree * tree;
	FILE *suspix;
	FILE *ded_N;
	FILE *nsy_N;
	FILE *ded_E;
	FILE *nsy_E;
	if (strcmp(mode,"calibrate") == 0 )
	{

		rocid = atoi(argv[2]);
		runid = atoi(argv[3]);
		input_file = Form("data/vcal-s%d-run%d.root",rocid, runid);	

		output_file = Form("output/%d/calibration/%d/calibration_s%d_run%d.root",rocid,runid, rocid,runid);	
system(Form("mkdir -p output/%d/calibration/%d",rocid,runid));
		output_file_coeficients = Form("output/%d/calibration/%d/calibration_coeficients_s%d_run%d.txt",rocid,runid, rocid, runid);
		printf ("Calibrating ...\n");
		printf("Input file:  %s\n", input_file);
		printf("Output file: %s\n", output_file);
		printf("Coeficients: %s\n", output_file_coeficients);
		fcoefs =fopen(output_file_coeficients,"w");
		suspix=fopen(Form("output/%d/calibration/%d/suspix_s%d_run%d.txt", rocid, runid, rocid, runid),"w");
	}
	else if (strcmp(mode,"scurve") == 0 )
	{
		
		rocid = atoi(argv[2]);
		runid = atoi(argv[3]);
		input_file = Form("data/scurve-s%d-run%d.root",rocid, runid);	
		
		output_file = Form("output/%d/scurve/%d/scurve_s%d_run%d.root",rocid,runid, rocid,runid);	
			
            system(Form("mkdir -p output/%d/scurve/%d/",rocid,runid));
//		output_file_coeficients = Form("output/%d/calibration/%d/calibration_coeficients_s%d_run%d.txt",rocid,runid, rocid, runid);
		printf ("Calibrating ...\n");
		printf("Input file:  %s\n", input_file);
		printf("Output file: %s\n", output_file);
//		printf("Coeficients: %s\n", output_file_coeficients);
//		fcoefs =fopen(output_file_coeficients,"w");
	suspix=fopen(Form("output/%d/scurve/%d/suspix_s%d_run%d.txt", rocid, runid, rocid, runid),"w");
	}	
	else if (strcmp(mode,"analyze") == 0 )
	{
		rocid = atoi(argv[2]);
		runid = atoi(argv[3]);
		calibration_run = atoi(argv[4]);
		input_file = Form("data/sr90-s%d-run%d.root",rocid, runid);	
		output_file = Form("output/%d/analysis/%d/analysis_s%d_run%d.root",rocid,runid, rocid, runid);	
                system(Form("mkdir -p output/%d/analysis/%d/",rocid,runid));
		calibration_coeficients_file = Form("output/%d/calibration/%d/calibration_coeficients_s%d_run%d.txt",rocid,calibration_run,rocid, calibration_run);
		printf ("Analyzing ...\n");
		printf("Input file:  %s\n", input_file);
		printf("Output file: %s\n", output_file);
		printf("Coeficients: %s\n", calibration_coeficients_file);
		
		nsy_N=fopen(Form("output/%d/analysis/%d/noisy_s%d_run%d.txt", rocid, runid, rocid, runid),"w");
		ded_N=fopen(Form("output/%d/analysis/%d/dead_s%d_run%d.txt", rocid, runid, rocid, runid),"w");	
		nsy_E=fopen(Form("output/%d/analysis/%d/energy_noisy_s%d_run%d.txt", rocid, runid, rocid, runid),"w");
		ded_E=fopen(Form("output/%d/analysis/%d/energy_dead_s%d_run%d.txt", rocid, runid, rocid, runid),"w");	

	}
	f = new TFile(input_file); //data/sr90-s16-run2.root
	tree = (TTree*)gDirectory->Get("t");
	
	study t(tree, output_file);

	if (strcmp(mode,"calibrate") == 0 ){
	  t.calibrate();
		for (int ir=0;ir< ir_max;ir++) {
			for (int ic=0;ic<ic_max;ic++) {	
			  if (_adc[ic][ir][0] != 0 )
			    {
			  fprintf(fcoefs,"%d %d %f %f %f\n",ic, ir, a0[ic][ir], a1[ic][ir], a2[ic][ir]); 
			    }
			  if(sus[ic][ir]>0) 
			    {
			      fprintf(suspix,"%d %d \n", ic, ir);
			    }
			}
		}
		fclose(fcoefs);
		fclose(suspix);
	}
	if (strcmp(mode,"scurve") == 0 ){
		t.scurve();
		
		//fwhm->GetEntries();
		//outliers
	

		for (int ir=0;ir< ir_max;ir++) {
		  for (int ic=0;ic<ic_max;ic++) {	
		    if (sus[ic][ir]>0 )
		      {	
			fprintf(suspix,"%d %d \n", ic, ir);
			//printf("%d",sus[ir][ic]);
		      }
		    
		  }
		}
	fclose(suspix);
	}
	
	if (strcmp(mode,"analyze") == 0 ){
		printf("Analyzing\n");
		printf("here");
		ifstream fcal;
		fcal.open(calibration_coeficients_file);
		processcal(fcal);
		printf("Read coeficients\n");
		t.analyze();
		printf("Done\n");
		

      		//float deltaX = (max_row - min_row)/n_row;
		//float deltaY = (max_col - min_col)/n_col;

		float col=0;
		float row=0;
		int bin=0;
	for (int ir=0;ir< ir_max;ir++) {
		  for (int ic=0;ic<ic_max;ic++) {

		    //convert bin number to row, col
		    //  float row = min_row + deltaX*(0.5 + (float)ir);
		    // float col = min_col + deltaY*(0.5 + (float)ic);
		    //these mark the pixel as being at the center of the bin,
		    //which probably won't be on an integer row number
		    
		        if (noisy_N[ic][ir]>0 )
		      {	
			//bin=noisy_N[ic][ir];
			//occupancy_noisy->Draw();
			//printf("%d", bin);
			fprintf(nsy_N,"%d %d \n", ic, ir);
			
			//printf("%d",sus[ir][ic]);
		      }
		    if(dead_N[ic][ir]>0)
		      {	
			bin=dead_N[ic][ir];
			fprintf(ded_N,"%d %d \n", ic, ir);
			//printf("%d",sus[ir][ic]);
			}

		    if (noisy_E[ic][ir]>0 )
		      {	
			bin=noisy_E[ic][ir];
			fprintf(nsy_E,"%d %d \n", ic, ir);
			//printf("%d",sus[ir][ic]);
		      }
		    if(dead_E[ic][ir]>0)
		      {
			bin=dead_E[ic][ir];
			
			fprintf(ded_E,"%d %d \n", ic,ir );
			//printf("%d",sus[ir][ic]);
			}
		  }//end for ic
	}//end for ir

	fclose(nsy_N); 
	fclose(ded_N);
	fclose(nsy_E); 
	fclose(ded_E);
	}

    //	out->Close();
	return 0;
}






	
TH2F *get_efficiency(TH2F *in){
	TH2F * efficiency = (TH2F*) in->Clone();
	efficiency->Reset();	
	efficiency->SetName(Form("%s_efficiency",efficiency->GetName()));
	efficiency->Sumw2();
	for(int ic = 1; ic<=in->GetNbinsX();ic++)
		for(int ir = 1; ir<=in->GetNbinsY();ir++)
		{
			float norm = 0;
			float ieff = 0;
			float ieff_err = 0;
			int cnt = 0;
			double neighbours[8];
			for (int i=-1; i<=1; i++) 
			{
				for (int j=-1; j<=1; j++) 
				{
					if (i==0 && j==0) continue;
					neighbours[cnt] = in->GetBinContent(ic+i, ir+j);
					cnt++;
				}
			}
			
			norm = TMath::Mean(8,neighbours);
			if (norm == 0)
			{
				ieff =0;
				ieff_err = 0;
			}
			else{ 
				ieff = in->GetBinContent(ic,ir)/norm;
				ieff_err = TMath::RMS(8,neighbours)/norm;
			
//			printf ("%.2g -- %2g +/- %.2g\n", in->GetBinContent(ic,ir), ieff, ieff_err);
			efficiency->SetBinContent(ic,ir,ieff);
			efficiency->SetBinError(ic,ir,ieff_err);
			}
		}
	
	return efficiency;
	
}

TH1F * get_index(TH2F *in){
	//	char *name = in->GetTitle()+"1D";
	TH1F *temp = new TH1F("","", 1000, 0, 1000);
	temp->SetName(Form("%s_index",in->GetName()));

	int index = 0;
	for(int ic = 1; ic<=in->GetNbinsX();ic++)
		for(int ir = 1; ir<=in->GetNbinsY();ir++)
		{
			index++;
			temp->SetBinContent(index, in->GetBinContent(ic, ir));
			temp->SetBinError(index, in->GetBinError(ic, ir));
		}
	
	return temp;
}

bool is_noisy(int _ic, int _ir){
	
	bool noise = false;
	noise = 
	(_ic == 0 && _ir==0) ||
	(_ic == 0 && _ir==0);
	return noise;
};

bool is_fiducial(int _ic, int _ir){
	
	bool _edge = false;
	_edge = 
	_ic < 15 ||
	_ic > 35 ||
	_ir < 42 ||
	_ir >= 77;

	return !_edge;
};
void study::processcal(ifstream& file)
{
	int ccc=0;
	int column=-1;
	int row=-1;
	double para1=-1;
	double para2=-1;
	double para3=-1;
	double para4=-1;
	while (!file.fail())
    {      double temp=-1;
		file>>temp;
		if(ccc%5==0) {column=temp;}
		else if(ccc%5==1) {row=temp;}
		else if(ccc%5==2) {para1=temp;}
		else if(ccc%5==3) {para2=temp;}
		else 
		{
			para3=temp;
			caldata[column][row][0]=para1;
			caldata[column][row][1]=para2;
			caldata[column][row][2]=para3;
			//			h_calibrated_region->Fill(column,row);
			
			//	  cout << "column = " << column << "; row = " << row << "  p1 =  " << para1 << "   p2 = " << para2 << " p3 = " << para3 << endl;
			//	  caldata[column][row][3]=para4;	  
		}
		ccc++;
    }
}

double erf(double *x, double *par){
	return 20.*TMath::Erf((x[0]-par[1])/par[2]);
//	return 20*TMath::Erf((x[0]-41)/0.5);
}

Double_t turnon_func(Double_t *x, Double_t *par)
{
	double halfpoint = par[0];
	double slope = par[1]/(2*log(2.));
	double plateau = par[2];
	par[2] = 20;
	
	//double offset = par[3];
	//double plateau = 1.0;
	double offset = 0;
	
	double pt = TMath::Max(x[0],0.0);
	
	double arg = 0;
	//cout << pt <<", "<< halfpoint <<", " << slope <<endl;
	arg = (pt - halfpoint)/(slope);
	double fitval = offset+0.5*plateau*(1+TMath::Erf(arg));
	return fitval;
}

void study::fill_scurve(int ic,int ir){
	TF1 *errf = new TF1("errf",turnon_func, 20, 60, 3);
//	erf->SetParameter(1, 41);
	TH1F *vcal_scurve = new TH1F ("vcal_scurve","", 100, 0.5,100.5);
	vcal_scurve->Sumw2();
	//	ofstream suspix(Form("output/%d/scurve/%d/suspix_s%d_run%d.txt", rocid, runid, rocid, runid));
	//FILE *boo;
	//boo=fopen(Form("output/%d/scurve/%d/boo.txt", rocid, runid),"w");
      	
	// Fill Vcal s cirve
	if (debug) 
		printf("adc_scurve>Fill():\n");
	
	for ( int i = 0; i<_ievent[ic][ir]; i++ ){	
//		if (debug)
//			printf("%d %d\n",i, int(_vcal[ic][ir][i]));
		vcal_scurve->Fill(_vcal[ic][ir][i]);
	}
	vcal_scurve->SetName(Form("ic_%d_ir_%d_rate_vs_vcal_scurve",ic,ir));
	if(debug) 
		printf("\n");
	errf->SetParameters(30,1,20);
	vcal_scurve->Fit("errf","rQ","",20, 100);
	
	double _chi2s=errf->GetChisquare();
	
	double _ndfs= errf->GetNDF();
	double sigma1 =0;
	double sigma2=0;
	//sig1=chi2s->GetRMS();
	//sig2=chi2ndfs->GetRMS();
	
	if(_chi2s>20 || (_chi2s/_ndfs)>1)
	{
	    	chi2s->Fill(-1);
		chi2ndfs->Fill(-1);
	anompix->SetBinContent(anompix->GetXaxis()->FindBin(ic),anompix->GetYaxis()->FindBin(ir), 1);
		chi2_2ds->SetBinContent(chi2_2ds->GetXaxis()->FindBin(ic), chi2_2ds->GetYaxis()->FindBin(ir), -1);
	chi2ndf_2ds->SetBinContent(chi2_2ds->GetXaxis()->FindBin(ic), chi2_2ds->GetYaxis()->FindBin(ir), -1);
	sus[ic][ir]=1;

	}
	else
	  {
	chi2s->Fill(_chi2s);
	chi2ndfs->Fill(_chi2s/_ndfs);
	  
	chi2_2ds->SetBinContent(chi2_2ds->GetXaxis()->FindBin(ic), chi2_2ds->GetYaxis()->FindBin(ir), _chi2s);
	chi2ndf_2ds->SetBinContent(chi2_2ds->GetXaxis()->FindBin(ic), chi2_2ds->GetYaxis()->FindBin(ir), _ndfs);
	  }
	//	printf("adc_mean[vcal_array_bins-1] = %f\n",adc_mean[vcal_array_bins-1]-0.1*(adc_mean[vcal_array_bins-1]-adc_mean[0]));
	a0[ic][ir]=(errf->GetParameter(0));
	//printf("value=%.3g\n", (a0[ic][ir]) *65);
	a1[ic][ir]=(errf->GetParameter(1) );
	
	if(fabs((a0[ic][ir])*65)>300 ||fabs((a1[ic][ir])*65)>100)
	  {
	    halfpoint_h ->Fill(-1);
	    fwhm_h    ->Fill(-1);
	anompix->SetBinContent(anompix->GetXaxis()->FindBin(ic),anompix->GetYaxis()->FindBin(ir), 1);
	halfpoint2d_h->SetBinContent(halfpoint2d_h->GetXaxis()->FindBin(ic),halfpoint2d_h->GetYaxis()->FindBin(ir), (a0[ic][ir]) *65);
	fwhm2d_h    ->SetBinContent(fwhm2d_h->GetXaxis()->FindBin(ic),fwhm2d_h->GetYaxis()->FindBin(ir), (a1[ic][ir]) *65);

	sus[ic][ir]=1;
	  }

	else
	{
	halfpoint_h ->Fill((a0[ic][ir])*65);
	fwhm_h->Fill((a1[ic][ir])*65);
	//printf("%.3g %.3g\n", a0[ic][ir]*65 , a1[ic][ir]);
	halfpoint2d_h->SetBinContent(halfpoint2d_h->GetXaxis()->FindBin(ic),halfpoint2d_h->GetYaxis()->FindBin(ir), a0[ic][ir]*65);
	fwhm2d_h    ->SetBinContent(fwhm2d_h->GetXaxis()->FindBin(ic),fwhm2d_h->GetYaxis()->FindBin(ir),a1[ic][ir]*65);
	//sus[ic][ir]=0;
       	  }
	//	a2[ic][ir]=f1->GetParameter(2);
	vcal_scurve->Write();


	delete vcal_scurve;
	return;
}



void study::fill_hist(int ic,int ir){
	
	//creates a function to fit the hist... pol2 for 2 degree polynomial.. pol3 for 3
	f1 = new TF1("f1","pol2");
	f2 = new TF1("f2","pol2");
	
	double adc_mean[1500];
	double adc_rms[1500];
	double vcal_array[1500];
	double vcale_array[1500];
	//	double a[1500];
	//int adc_cnt[1500];
	double _adc_mean[1500];
	double _adc_rms[1500];
	int _vcal_array[1500];
	//initialization
	//for(int j=0; j<1500; j++)
	//{
        int vcal_i=0;
	int vcal_im1=0;
	
	double a[1500][50];


	for(int i=0;i<1500;i++){
	  // adc[i] = 0;
	  _adc_mean[i]=0;
	  _adc_rms[i] = 15;
	  _vcal_array[i] = 0;
	  // adc[i]=0;
	  // adc_cnt[i]=0;
	}
	//}
	TGraph *gr = new TGraph(_ievent[ic][ir], _adc[ic][ir],_vcal[ic][ir]);
	
	int adc_i=0;
	int j=0;
	//	int vcal=0;
	double mean=0;
	double cmax=0;
	for ( int i = 0; i<_ievent[ic][ir]; i++ ){
		if(_vcal[ic][ir][i]>1499) continue;
		
		_vcal_array[_vcal[ic][ir][i]]+=1;
		//j=_vcal_array[_vcal[ic][ir][i]]-1;
		a[_vcal[ic][ir][i]][_vcal_array[_vcal[ic][ir][i]]-1]=_adc[ic][ir][i];
		
		//printf("%f %d\n",a[_vcal[ic][ir][i]][_vcal_array[_vcal[ic][ir][i]]-1],_vcal[ic][ir][i] );
		//_adc_mean[_vcal[ic][ir][i]] +=  _adc[ic][ir][i];
		//  	_adc_rms[_vcal[ic][ir][i]] = ( _adc[ic][ir][i]);
		//if (debug) 	printf("%d %d\n",int(_vcal[ic][ir][i]), int(_vcal_array[_vcal[ic][ir][i]]));
		//		h2->Fill(adc[ic][ir][i],vcal[ic][ir][i]);
	

 
	}
	

	int _cnt =0;
	
	for(int i=0;i<1500;i++){
		if(_vcal_array[i]==0 ) continue;
		_adc_mean[i]=TMath::Mean(_vcal_array[i],a[i]);
		_adc_rms[i]=TMath::RMS(_vcal_array[i],a[i]);
		//_adc_mean[i]/=_vcal_array[i];
		//adc_mean[i]=_adc_mean[i];
		//adc_rms[i]=_adc_rms[i];
		//	adc_rms[i]=TMath::RMS(adc_cnt[i], adc[i]);
		//_adc_mean[i]=TMath::Mean( _adc[ic][ir]);
		// _adc_rms[i]+=((_adc_mean[i]-_adc_rms[i])*(_adc_mean[i]-_adc_rms[i]))/_vcal_array[i];
		//printf("%d %d \n", _adc_mean[i],  _adc_rms[i]);		
		adc_mean[_cnt] = _adc_mean[i];
		adc_rms[_cnt] = _adc_rms[i];
	   
		vcal_array[_cnt] = i;
		vcale_array[_cnt] = 65*i;
		
		_cnt++;
	}
       
	int _vcal_array_bins = _cnt;
	
	int	cnt = 0;
	int vcal_array_bins = _vcal_array_bins;
	if(debug) printf("2: vcal [] = ");
	if(debug) 
		for (int i =0;i<vcal_array_bins;i++){
			printf("%d ",int(vcal_array[i]));
		}		
	if(debug) printf("\n");
	
	if(debug) 	printf("2: adc  [] = ");
	if(debug) 
		for (int i =0;i<vcal_array_bins;i++){
			printf("%d ",int(adc_mean[i]));
		}	if(debug) printf("\n");
	
	residuals_h = new TH1F("residuals","",100, -50,50); // renamed every pixel
	p =new TGraphErrors(vcal_array_bins,adc_mean,vcal_array,adc_rms,0);
	pe =new TGraphErrors(vcal_array_bins,adc_mean,vcale_array,adc_rms,0);
	
	if(vcal_array_bins < 5) {
		printf("Warning! vcal_array_bins = %d <20\n");
		return;
	}
	
	double low_adc = adc_mean[0];
	double high_adc = adc_mean[vcal_array_bins-1];
	double offset = 0.1* (high_adc-low_adc);
	
	if(debug) 	printf ("range = [%d,%d]\n",int(low_adc),int(high_adc));
	p->Fit("f1","Q","",low_adc, high_adc-offset);//adc_mean[0], 0.5*adc_mean[vcal_array_bins-1]);
	pe->Fit("f2","Q","",low_adc, high_adc-offset);//adc_mean[0], 0.5*adc_mean[vcal_array_bins-1]);
	
	//	printf("adc_mean[vcal_array_bins-1] = %f\n",adc_mean[vcal_array_bins-1]-0.1*(adc_mean[vcal_array_bins-1]-adc_mean[0]));
	a0[ic][ir]=f1->GetParameter(0);
	a1[ic][ir]=f1->GetParameter(1);
	a2[ic][ir]=f1->GetParameter(2);
	
	
	minimal_vcal = vcal_array[0];	
	
	for (int i =0;i<vcal_array_bins;i++){
		
		if (adc_mean[i] <= high_adc-offset) saturated_vcal= vcal_array[i];
		if(adc_mean[i]>=low_adc || adc_mean[i]<=high_adc-offset){
			residuals_h->Fill(vcal_array[i]- f1->Eval(adc_mean[i]));
		}
	}
	
	saturated_vcal_h->Fill(65*saturated_vcal);
	minimal_vcal_h->Fill(65*minimal_vcal);
	double _chi2 = f1->GetChisquare();
	if (_chi2>chi2->GetXaxis()->GetBinCenter(chi2->GetNbinsX()))
		_chi2 = chi2->GetXaxis()->GetBinCenter(chi2->GetNbinsX());
	
	double _ndf = f1->GetNDF();
	chi2->Fill(_chi2);
	chi2_2d->SetBinContent(chi2_2d->GetXaxis()->FindBin(ic), chi2_2d->GetYaxis()->FindBin(ir),_chi2);
	
	if(f1->GetNDF()!=0) {
		chi2ndf->Fill(_chi2/_ndf);
		chi2ndf_2d->SetBinContent(chi2ndf_2d->GetXaxis()->FindBin(ic),chi2ndf_2d->GetYaxis()->FindBin(ir),_chi2/_ndf);
	}
	else {
		chi2ndf->Fill(0);
		chi2ndf_2d->SetBinContent(chi2ndf_2d->GetXaxis()->FindBin(ic),chi2ndf_2d->GetYaxis()->FindBin(ir), 0);
	}




	if(_chi2>700 || (_chi2/_ndf)>30)
	{
	  sus[ic][ir]=1;
	  anompix->SetBinContent(chi2ndf_2d->GetXaxis()->FindBin(ic),chi2ndf_2d->GetYaxis()->FindBin(ir), 1);

	}
	else
	  {
	     sus[ic][ir]=0;
	     anompix->SetBinContent(chi2ndf_2d->GetXaxis()->FindBin(ic),chi2ndf_2d->GetYaxis()->FindBin(ir), 0);
	  }
	
	f1->SetName(Form("ic_%d_ir_%d_fit",ic,ir));
	residuals_h->SetName(Form("ic_%d_ir_%d_residual",ic,ir));	
	gr->SetName(Form("ic_%d_ir_%d_adc_vs_vcal",ic,ir));
	p->SetName(Form("ic_%d_ir_%d_adc_vs_vcal_graph",ic,ir));
	pe->SetName(Form("ic_%d_ir_%d_adc_vs_vcale_graph",ic,ir));
	
	p->Write();
	pe->Write();
	residuals_h->Write();
	gr->Write();
	
	delete p;
	delete residuals_h;
	delete gr;
	delete f1;
	
	return;
}



void study::calibrate(){
	// Book all histograms
	for (int i =0;i<ic_max;i++)
	{
		for (int j =0;j<ir_max;j++)
			for (int e =0;e<max_events;e++)
			{
				_adc[i][j][e] = 0;
				_ievent[i][j]=-1;
			}
	}
	
	
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	if(debug) printf("nentries = %d\n",nentries);
	Long64_t nbytes = 0, nb = 0;
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		//		if(jentry % 100000 ==0) printf("");
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry ==0)
			for(int i = 0;i<5;i++){
			  levelinfo_cut_h->SetBinContent(levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), nentries);
				if(debug) printf("cut[%d] = %f\n",levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), levelinfo_cut[i]);
			}
		
		// Event Level Infos
		{
			eventhits->Fill(nhits);
			
			ub1_h->Fill(ub[0]);
			ub2_h->Fill(ub[1]);
			ub3_h->Fill(ub[2]);
			ubROC1_h->Fill(ubROC1);
			ubdiff->Fill(ub[0]-ub[1]);
			ub_tr1_h->Fill(ub_tr1);
			ub_tr2_h->Fill(ub_tr1);
			tbmdiff->Fill(tbm_status[1]-tbm_status[3]);
			tbm1->Fill(tbm_status[0]);
			tbm2->Fill(tbm_status[1]);
			tbm3->Fill(tbm_status[2]);
			tbm4->Fill(tbm_status[3]);
			
			ub1t_h->Fill(tick,ub[0]);
			ub2t_h->Fill(tick,ub[1]);
			ub3t_h->Fill(tick,ub[2]);
			
			ubROC1t_h->Fill(tick,ubROC1);
			ubdifft->Fill(tick,ub[0]-ub[1]);
			ub_tr1t_h->Fill(tick,ub_tr1);
			ub_tr2t_h->Fill(tick,ub_tr1);
			
			
			//		printf("%d %d\n",tick, ub[0]);
			//		printf("levelinfo = ");
			c1->Fill(levelinfo[0]);
			c0->Fill(levelinfo[1]);
			r2->Fill(levelinfo[2]);
			r1->Fill(levelinfo[3]);
			r0->Fill(levelinfo[4]);
			
			for(int i = 0;i<5;i++){
				levelinfo_h->Fill(levelinfo[i]);
				//			levelinfo_cut_h->Fill(levelinfo_cut[i]);
				levelinfot_h->Fill(tick, levelinfo[i]);
				//			printf("%i ",levelinfo[i]);
			}
		}
		
		if (nhits>1) {printf("More than one hit per charge injection!\nSkipping the entry.\n");continue;}
		
		int nhits_fiducial=0;
		int _ic = ic[0];
		int _ir = ir[0];
		
		h_number_of_hits_map->Fill(_ic,_ir);
		
		++_ievent[_ic][_ir];
		
		_adc[_ic][_ir][_ievent[_ic][_ir]]   = adc[0];
		
		if (range == 1)
			_vcal [_ic][_ir][_ievent[_ic][_ir]]  = 7*vcal;
		else
			_vcal [_ic][_ir][_ievent[_ic][_ir]]  = vcal;
		
		_range[_ic][_ir][_ievent[_ic][_ir]] = range;
	}
	
	// Now, once the calibration data has been written,
	// it is time to determine coeficients
	
	for (int ir=0;ir<100;ir++) {
		
		for (int ic=0;ic<100;ic++) {
//			if(debug) 
//				if(ir!=13) continue;
			
			if(_vcal[ic][ir][0]==0) continue;
			if(debug) printf("vcal(%d,%d) = %d\n",ic,ir,int(_vcal[ic][ir][0]));
			printf("Calibrating pixel (%d, %d): ",ic,ir);
			fill_hist(ic,ir);
			printf("Done\n");
			
		}
	}
	chi2->GetRMS();
	out->cd();
	//	ClusterE->Draw();
	out->Write();
	out->Close();	
}

void study::scurve(){
	// Book all histograms
	for (int i =0;i<ic_max;i++)
	{
		for (int j =0;j<ir_max;j++)
			for (int e =0;e<max_events;e++)
			{
				_adc[i][j][e] = 0;
				_ievent[i][j]=-1;
			}
	}
	
	
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	if(debug) printf("nentries = %d\n",nentries);
	Long64_t nbytes = 0, nb = 0;
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		//		if(jentry % 100000 ==0) printf("");
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry ==0)
			for(int i = 0;i<5;i++){
				levelinfo_cut_h->SetBinContent(levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), nentries);
				if(debug) printf("cut[%d] = %f\n",levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), levelinfo_cut[i]);
			}
		
		// Event Level Infos
		{
			eventhits->Fill(nhits);
			
			ub1_h->Fill(ub[0]);
			ub2_h->Fill(ub[1]);
			ub3_h->Fill(ub[2]);
			ubROC1_h->Fill(ubROC1);
			ubdiff->Fill(ub[0]-ub[1]);
			ub_tr1_h->Fill(ub_tr1);
			ub_tr2_h->Fill(ub_tr1);
			tbmdiff->Fill(tbm_status[1]-tbm_status[3]);
			tbm1->Fill(tbm_status[0]);
			tbm2->Fill(tbm_status[1]);
			tbm3->Fill(tbm_status[2]);
			tbm4->Fill(tbm_status[3]);
			
			ub1t_h->Fill(tick,ub[0]);
			ub2t_h->Fill(tick,ub[1]);
			ub3t_h->Fill(tick,ub[2]);
			
			ubROC1t_h->Fill(tick,ubROC1);
			ubdifft->Fill(tick,ub[0]-ub[1]);
			ub_tr1t_h->Fill(tick,ub_tr1);
			ub_tr2t_h->Fill(tick,ub_tr1);
			
			
			//		printf("%d %d\n",tick, ub[0]);
			//		printf("levelinfo = ");
			c1->Fill(levelinfo[0]);
			c0->Fill(levelinfo[1]);
			r2->Fill(levelinfo[2]);
			r1->Fill(levelinfo[3]);
			r0->Fill(levelinfo[4]);
			
			for(int i = 0;i<5;i++){
				levelinfo_h->Fill(levelinfo[i]);
				//			levelinfo_cut_h->Fill(levelinfo_cut[i]);
				levelinfot_h->Fill(tick, levelinfo[i]);
				//			printf("%i ",levelinfo[i]);
			}
		}
		
		if (nhits>1) {printf("More than one hit per charge injection!\nSkipping the entry.\n");continue;}
		
		int nhits_fiducial=0;
		int _ic = ic[0];
		int _ir = ir[0];
		
		h_number_of_hits_map->Fill(_ic,_ir);
		
		++_ievent[_ic][_ir];
		
		_adc[_ic][_ir][_ievent[_ic][_ir]]   = adc[0];
		
		if (range == 1)
			_vcal [_ic][_ir][_ievent[_ic][_ir]]  = 7*vcal;
		else
			_vcal [_ic][_ir][_ievent[_ic][_ir]]  = vcal;
		
		_range[_ic][_ir][_ievent[_ic][_ir]] = range;
	}
	
	// Now, once the calibration data has been written,
	// it is time to determine coeficients
	
	for (int ir=0;ir<100;ir++) {
		for (int ic=0;ic<100;ic++) {
			if(debug) 
				if(ic!=20 || ir !=40) continue;
			
			if(_vcal[ic][ir][0]==0) continue;
			if(debug) printf("vcal(%d,%d) = %d\n",ic,ir, int(_vcal[ic][ir][0]));
			fill_scurve(ic,ir);
		}
	}
	
	out->cd();
	//	ClusterE->Draw();
	//	out->Write();
	//	out->Close();	
}

void study::analyze(){
 	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;
	if(debug) nentries = 10;
	if(debug) printf("nentries = %d\n",nentries);
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry ==0) {
			for(int i = 0;i<5;i++){
				levelinfo_cut_h->SetBinContent(levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), nentries);
				if(debug)	printf("cut[%d] = %f\n",levelinfo_cut_h->GetXaxis()->FindBin(levelinfo_cut[i]), levelinfo_cut[i]);
			}
		}
		if(debug) printf("nhits = %d\n",nhits);
		
		eventhits->Fill(nhits);
		if(debug) printf("nhits = %d\n",nhits);
		
		ub1_h->Fill(ub[0]);
		ub2_h->Fill(ub[1]);
		ub3_h->Fill(ub[2]);
		ubROC1_h->Fill(ubROC1);
		ubdiff->Fill(ub[0]-ub[2]);
		tbmdiff->Fill(tbm_status[1]-tbm_status[3]);
		tbm1->Fill(tbm_status[0]);
		tbm2->Fill(tbm_status[1]);
		tbm3->Fill(tbm_status[2]);
		tbm4->Fill(tbm_status[3]);
		ub_tr1_h->Fill(ub_tr1);
		ub_tr2_h->Fill(ub_tr1);
		
		ub1t_h->Fill(tick,ub[0]);
		ub2t_h->Fill(tick,ub[1]);
		ub3t_h->Fill(tick,ub[2]);
		
		ubROC1t_h->Fill(tick,ubROC1);
		ubdifft->Fill(tick,ub[0]-ub[1]);
		ub_tr1t_h->Fill(tick,ub_tr1);
		ub_tr2t_h->Fill(tick,ub_tr1);
		
		//		printf("%d %d\n",tick, ub[0]);
		//		printf("levelinfo = ");
		c1->Fill(levelinfo[0]);
		c0->Fill(levelinfo[1]);
		r2->Fill(levelinfo[2]);
		r1->Fill(levelinfo[3]);
		r0->Fill(levelinfo[4]);
		
		for(int i = 0;i<5;i++)
		{
			levelinfo_h->Fill(levelinfo[i]);
			levelinfot_h->Fill(tick, levelinfo[i]);
		}
		
		if(debug) printf("ubROC1 = %d\n",ubROC1);
		if(debug) printf("ub_tr1 = %d\n",ub_tr1);
		if (nhits==0) {printf("No hits!!! in an event!!! Investigate!\n");continue;}
		int nhits_fiducial=0;
		for (int ihit = 0;ihit<nhits;ihit++){
			int i = ic[ihit];//
			int j = ir[ihit];
			double signal = adc[ihit];
			double electrons =  65.*(caldata[i][j][0] + (signal)*caldata[i][j][1] + (signal)*(signal)*caldata[i][j][2]);
			pixel_energy->Fill(i,j, electrons);
			if(electrons==0) {
				//				printf("d1:%d %d %f %f %f\n", i, j,caldata[i][j][0],caldata[i][j][1],caldata[i][j][2]);
				continue;
			}
			h_number_of_hits_map->Fill(i,j);
			h_energy_hits_map->Fill(i,j, electrons);
			if (debug) printf("d2: %i %i %.2g %.2g\n",i,j,signal,electrons);  
			RawADC_h->Fill(signal);
		}
		
		//Find most energetic pixel
		double max_e = 0;
		int max_i = 0;
		for (int ihit = 0;ihit<nhits;ihit++){ //search for the most energetic pix			
			double signal = adc[ihit];
			int i = ic[ihit];
			int j = ir[ihit];
			
			double electrons =  65.*(caldata[i][j][0] + (signal)*caldata[i][j][1] + (signal)*(signal)*caldata[i][j][2]);
			if(electrons<0) {
				if(i == 34 && j == 77) printf("d3:e(adc = %f) = %f %d %d %f %f %f\n",signal, electrons, i, j,caldata[i][j][0],caldata[i][j][1],caldata[i][j][2]);
				continue;
			}
			if (electrons > max_e) {
				max_e = electrons;
				max_i = ihit;
			}
		}



		//Find Dead and Noisy pixels		
		/*	Find_Dead_And_Noisy(h_number_of_hits_map, *occupancy_noisy, *occupancy_dead);//xxxx
	Find_Dead_And_Noisy(h_energy_hits_map, *energetically_noisy, *energetically_dead);

		for(int i=1;i<=h_number_of_hits_map->GetNbinsX();++i)
		  for(int j=1;j<=h_number_of_hits_map->GetNbinsY();++j){
		    occupancy_all->Fill(h_number_of_hits_map->GetBinContent(i,j));
		    dead_N[i][j]=occupancy_dead->GetBinContent(i,j);
		    noisy_N[i][j]=occupancy_noisy->GetBinContent(i,j);
		  }//endfor

		for(int i=1;i<=h_energy_hits_map->GetNbinsX();++i)
		  for(int j=1;j<=h_energy_hits_map->GetNbinsY();++j){
		    dead_E[i][j]=energetically_dead->GetBinContent(i,j);
		    noisy_E[i][j]=energetically_noisy->GetBinContent(i,j);
		  }//endfor
		printf("here");
	 
	// Derive efficiencies
	h_number_of_hits_map_efficiency =   get_efficiency(h_number_of_hits_map);
	h_number_of_hits_map_efficiency1d =  get_index(h_number_of_hits_map_efficiency);
	
	for(int i =1;i<= h_number_of_hits_map_efficiency1d->GetNbinsX();i++){
		float bin = h_number_of_hits_map_efficiency1d_all->Fill(h_number_of_hits_map_efficiency1d->GetBinContent(i));
		//		if(bin!=0.0)
		h_number_of_hits_map_efficiency1d_all->Fill(bin);
	}
	
	h_energy_hits_map_efficiency =   get_efficiency(h_energy_hits_map);
	h_energy_hits_map_efficiency1d =  get_index(h_energy_hits_map_efficiency);

	for(int i =1;i<= h_energy_hits_map_efficiency1d->GetNbinsX();i++){
		float bin = h_energy_hits_map_efficiency1d_all->Fill(h_energy_hits_map_efficiency1d->GetBinContent(i));
//		if(bin!=0)
			h_energy_hits_map_efficiency1d_all->Fill(bin);
	}*/
		
		//	printf("here");


	//	bool pix=false;


		// Fill cluster 
		double cluster_energy = max_e;
		double cluster_adc = adc[max_i];
		
		int cnt = 1;
		for (int ihit = 0;ihit<nhits;ihit++){
			if (ihit == max_i) continue;
			int i = ic[ihit];
			int j = ir[ihit];
			//	if(i==31 && j==43)continue;
			//	if(dead_N[i][j]!=0 || noisy_N[i][j]!=0)continue;
			int    dx = abs(ic[ihit] - ic[max_i]);
			int    dy = abs(ir[ihit] - ir[max_i]);
			//			double dr = sqrt(dx*dx+dy*dy);
			/*	if(ic[max_i]==34 && ir[max_i]==43) {
			  cluster_adc =0; 
			  cluster_energy =0;
			  max_e=0;
			  pix=true;
			  break;}
			if(ic[max_i]==21 && ir[max_i]==60)
			  {
			  cluster_adc =0; 
			  cluster_energy =0;
			  max_e=0;
			  break;}
			if(ic[max_i]==29 && ir[max_i]==50)
			  {
			  cluster_adc =0; 
			  cluster_energy =0;
			  max_e=0;
			  break;}
			*/	
		//     if (dr > 2 ) continue;
			if(abs(ic[max_i]-ic[ihit])>2 || abs(ir[max_i]-ir[ihit])>2) continue;
			if (is_fiducial(i,j)) nhits_fiducial++;
			
			cnt ++;//increments everytime you see a pixel in the vacinity of the most energetic pix.
			double signal = adc[ihit];
			
			//			printf("%d %d %f %f %f\n", i, j,caldata[i][j][0],caldata[i][j][1],caldata[i][j][2]);
			double electrons =  65.*(caldata[i][j][0] + (signal)*caldata[i][j][1] + (signal)*(signal)*caldata[i][j][2]);
			//	if (ic[ihit]==20 && ir[ihit]==69) adcelectron_20_69->Fill(electrons);

			if(electrons==0) {
			  //			printf("d4: %d %d %f %f %f\n", i, j,caldata[i][j][0],caldata[i][j][1],caldata[i][j][2]);
			  continue;
			}
			//printf("electrons = %.2g\n",electrons);
			cluster_adc += signal; 
			cluster_energy += electrons;
		}
		
	       	//if (pix==true){printf("cluster energy %e \n", cluster_energy);}
		clusterhits->Fill(cnt);
		
		if (!is_fiducial(ic[max_i],ir[max_i])) continue;
		if(max_e == 0){ 
		 	printf("%d %d\n",ic[max_i], ir[max_i]);
		 	printf("No hits??? Something is wrong, investigate!\n"); 
			continue;
		} 
	      	//	if(ic[max_i]!=34 && ir[max_i]!=43)
		//  {
		//adcelectron_max->Fill(max_e);
		//}
		clusterhits_fiducial->Fill(nhits_fiducial);
			

		if (cnt == 1){ 

		     
		       adcelectron_one -> Fill(cluster_energy);

			  occupancy_one->Fill(ic[max_i],ir[max_i]);
		
			  //	if(ic[max_i]!=21 || ir[max_i]!=60)
			  // {occupancy_one->Fill(ic[max_i],ir[max_i]);}
			  
			  //if(ic[max_i]!=29 || ir[max_i]!=51)
			  // {occupancy_one->Fill(ic[max_i],ir[max_i]);}	
			adcelectron_time->Fill(tick, cluster_energy);
			clusterET1->Fill(tick, cluster_energy);//clusterET1->Fill(tick, cluster_energy);
			//printf("%d %d %d \n", ic[max_i], ir[max_i], cluster_energy); 
		}
		if (cnt == 2) {
		  // if(ic[max_i]!=34 && ir[max_i]!=43)
		  //  if(((ic[max_i]!=34 && ir[max_i]!=43)||(ic[max_i]!=21 && ir[max_i]!=60))||(ic[max_i]!=29 && ir[max_i]!=51))
			   adcelectron_two -> Fill(cluster_energy);
			 //	  if(ic[max_i]!=21 && ir[max_i]!=60)
			 //	 {adcelectron_two -> Fill(cluster_energy);}
			 //	  if(ic[max_i]!=29 && ir[max_i]!=50)
			 //	 {adcelectron_two -> Fill(cluster_energy);}
			 //ok
			 occupancy_two->Fill(ic[max_i],ir[max_i]);//ok
			 clusterET2->Fill(tick, cluster_energy);
			
			 // if(((ic[max_i]!=34 && ir[max_i]!=43)||(ic[max_i]!=21 && ir[max_i]!=60))||(ic[max_i]!=29 && ir[max_i]!=51))	
			 //if(ic[max_i]!=34 && ir[max_i]!=43)
			 // {
			    ClusterADC->Fill(cluster_adc);
			    ClusterE->Fill(cluster_energy);
			    MaxADC->Fill(adc[max_i]);
			    MaxE->Fill(max_e);
			    // }	
			//clusterET2->Fill(tick, cluster_energy);
			
		}
		if (cnt >= 3) { 
			adcelectron_three -> Fill(cluster_energy); 
			occupancy_three->Fill(ic[max_i],ir[max_i]);
			clusterET3->Fill(tick, cluster_energy);
		}
		
		if (cnt >= 4) { 
			adcelectron_four -> Fill(cluster_energy); 
			//clusterET3->Fill(tick, cluster_energy);
		}
	

		
		//		ClusterADC->Fill(cluster_adc);
		//		ClusterE->Fill(cluster_energy);
	}

	/*
	for(int i=1;i<=h_number_of_hits_map->GetNbinsX();++i)
		for(int j=1;j<=h_number_of_hits_map->GetNbinsY();++j){

		        //this seems to duplicate h_number_of_hits_map
			occupancy_all->Fill(h_number_of_hits_map->GetBinContent(i,j));

			if (h_number_of_hits_map->GetBinContent(i,j) < 40)//if dead
			{	
			  occupancy_dead->SetBinContent(i,j,1);//mark it dead as a th2
			  dead_N[i][j]=1; //mark it dead as an array!!
			}
			else
			  {
			  dead_N[i][j]=0; //else, mark alive as an array
			  }
			if (h_number_of_hits_map->GetBinContent(i,j) > 1100) // if noisy
			{
			  occupancy_noisy->SetBinContent(i,j,1);//mark noisy
			   noisy_N[i][j]=1;//mark array noisy
			}
			else
			  {
			    noisy_N[i][j]=0;//else, mark not noisy
			  }
			//	printf("%d\n", h_number_of_hits_map->GetBinContent(i,j));
			}
	*/
	//Mark dead and noisy pixels

	Find_Dead_And_Noisy(h_number_of_hits_map, *occupancy_noisy, *occupancy_dead);//xxxx
	Find_Dead_And_Noisy(h_energy_hits_map, *energetically_noisy, *energetically_dead);

	int x=0;
	int y=0;
	double bin=0;

		for(int i=1;i<=h_number_of_hits_map->GetNbinsX();++i)
		  for(int j=1;j<=h_number_of_hits_map->GetNbinsY();++j)
		  {
		    occupancy_all->Fill(h_number_of_hits_map->GetBinContent(i,j));
		    bin=occupancy_dead->GetBinContent(i,j);
		    if(bin>0)
		      {
		    x=occupancy_dead->GetXaxis()->GetBinCenter(i);
		    y=occupancy_dead->GetYaxis()->GetBinCenter(j);
		    dead_N[x][y]=1;
		      }
		    
		    
		    bin=occupancy_noisy->GetBinContent(i,j);
		    if(bin>0)
		      {
		    x=occupancy_noisy->GetXaxis()->GetBinCenter(i);
		    y=occupancy_noisy->GetYaxis()->GetBinCenter(j);
		    noisy_N[x][y]=1;
		      }
		  }//endfor

		for(int i=1;i<=h_energy_hits_map->GetNbinsX();++i)
		  for(int j=1;j<=h_energy_hits_map->GetNbinsY();++j){

		    bin=energetically_dead->GetBinContent(i,j);
		    if(bin>0)
		      {  
		    x=energetically_dead->GetXaxis()->GetBinCenter(i);
		    y=energetically_dead->GetYaxis()->GetBinCenter(j);
		    dead_E[x][y]=1;
		      }
		  

		    bin=energetically_noisy->GetBinContent(i,j);
		    if(bin>0)
		      {
		    x=energetically_noisy->GetXaxis()->GetBinCenter(i);
		    y=energetically_noisy->GetYaxis()->GetBinCenter(j);
		    noisy_E[x][y]=1;
		      }
			  }//endfor
		
	 
	// Derive efficiencies
	h_number_of_hits_map_efficiency =   get_efficiency(h_number_of_hits_map);
	h_number_of_hits_map_efficiency1d =  get_index(h_number_of_hits_map_efficiency);
	
	for(int i =1;i<= h_number_of_hits_map_efficiency1d->GetNbinsX();i++){
	  float bin = h_number_of_hits_map_efficiency1d_all->Fill(h_number_of_hits_map_efficiency1d->GetBinContent(i));
		//		if(bin!=0.0)
	  h_number_of_hits_map_efficiency1d_all->Fill(bin);
	}
	
	h_energy_hits_map_efficiency =   get_efficiency(h_energy_hits_map);
	h_energy_hits_map_efficiency1d =  get_index(h_energy_hits_map_efficiency);

	for(int i =1;i<= h_energy_hits_map_efficiency1d->GetNbinsX();i++){
	  float bin = h_energy_hits_map_efficiency1d_all->Fill(h_energy_hits_map_efficiency1d->GetBinContent(i));
	  //		if(bin!=0)
			h_energy_hits_map_efficiency1d_all->Fill(bin);
	}
	/*	
	double eff=0;
	for(int i=1; i<=h_energy_hits_map_efficiency->GetNbinsX(); i++){
	  for(int j=1;j<=h_energy_hits_map_efficiency->GetNbinsY();++j)
		  {
		    
		    eff=h_energy_hits_map_efficiency->GetBinContent(i,j);
		    if(eff<0.5)markloweff->SetBinContent(i,j,1);
		    if(eff>1.5)markhigheff->SetBinContent(i,j,1);

		    
		  }
	}
	
	// Create 1D energy histrograms for 
	for(int i=1;i<=pixel_energy->GetNbinsX();++i)
		{
		  for(int j=1;j<=pixel_energy->GetNbinsY();++j)
		  {
		    x=pixel_energy->Project3D("x")->GetBinCenter(i);
		    y=pixel_energy->Project3D("y")->GetBinCenter(j);
		    if((x<12 ||x>39)||(y<39 || y>79))continue;
		    pixel_energy->GetXaxis()->SetRange(i,i);
		    pixel_energy->GetYaxis()->SetRange(j,j);
		    pix_e=(TH1F*)pixel_energy->Project3D("z");
		    pix_e->SetName(Form("Pix_e_ic%d_ir%d", x,y));
		    pix_e->SetTitle(Form("Pixel Energy Distrubtion ic:%d ir:%d", x,y));
	    //pix_e->Write();
		    
		  }
		}
	*/
		//out->cd();
	//	ClusterE->Draw();
	//out->Write();
	//	out->Close();	
}





