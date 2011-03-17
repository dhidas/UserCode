#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>


#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TFile.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TGraphAsymmErrors.h"

int bg_pars=3;
int mc_pars=6;

///defining fit functions
Double_t Func_bg(Double_t* x, Double_t* par)
{

  //par[2] and par[3] used to fix normalization
  //takes bg_pars+mc_pars + 2 parameters as input
  //return 1*par[0];
  return par[0]*TMath::Landau(x[0],par[1],par[2])+
    par[3]*TMath::Gaus(x[0],par[4],par[5]);

}
Double_t Func_data(Double_t* x, Double_t* par)
{

  //par[2] and par[3] used to fix normalization
  //takes bg_pars+mc_pars + 2 parameters as input
  //return 1*par[0];
  return par[0]*TMath::Landau(x[0],par[1],par[2])+
    par[3]*TMath::Gaus(x[0],par[4],par[5]);

}



void  CMS_significance(){
  //gRandom->SetSeed(std::fmod( (double) time(NULL),(double) 100000.));
  gRandom->SetSeed(1234567);

  ///////////////////////////////////////////////////////////
  //defining outputfile
  TFile *outputfile=new TFile("CMS_significance_allsmear_nobgfix.root","recreate");
  outputfile->mkdir("numevt");
  outputfile->mkdir("not_converged");

  //varibales for saving stuff in each pse
  std::vector <float > num_evt;
  std::vector<float > vmass;
  std::vector <float > num_evt_1sigdown;
  std::vector<float > num_evt_1sigup;
  std::vector <float > num_evt_2sigdown;
  std::vector<float > num_evt_2sigup;
  std::vector <float > num_evt_95;
  std::vector <float > num_evt_95_1sigdown;
  std::vector<float > num_evt_95_1sigup;
  std::vector <float > num_evt_95_2sigdown;
  std::vector<float > num_evt_95_2sigup;


  std::vector <float > vec_median;  
  std::vector <float > vec_95;  
  std::vector <float > vec_median_smear;  
  std::vector <float > vec_95_smear; 

  //define number of pseudoexperiments
  int nps=10000;
  //defining lumi
  float lumi=35;
  //defining expo parameters and errors
  double expo0=5.61166911473381180e+00; double expo0err=0.03*5.61166911473381180e+00;
  double expo1=-6.66717160493066584e-03;   double expo1err=0.05*6.66717160493066584e-03; 

  //defining the mass point we want to look at
  //float Mjjj=380;
  //defining exponential background function
  TF1 func_mcexpo("func_mcexpo", "expo(0)", 0, 800);
  //histograms to save outcome cross section (evt is acutally the xsec)
  TH1D *evt_1sigma=new TH1D("crossection","crossection",400,-300,300);
  TH1D *width=new TH1D("width","width",400,0,200);
  std::vector< float> vec_numobs_all;
  std::vector< float> vec_numobs_1sigma;
  std::vector< float> vec_numobs_2sigma;
  std::vector< float> vec_numobs_cuts;
  std::vector<float> vec_mass_all;
  std::vector<float> vec_mass_cuts;
  //loop for PES
  int count_pes_pass=0;
  int count_pes=0;
  for (float Mjjj=380; Mjjj<=385; Mjjj=Mjjj+10){
    for (int p=0; p<nps; p++){
      count_pes++;
      //defining acceptance as function of Mjjj
      double acceptance=-5.56179969055969892e-03 - 4.01623842089755741e-06 * Mjjj + 2.49580149009780901e-07 * Mjjj * Mjjj;
      //introducing a smear factor according to acceptance uncertainty (either Gauss, Uniform or no smearing)
      //float smearfactor=gRandom->Gaus(1,0.18);
      float smearfactor=1;
      // float smearfactor=gRandom->Uniform(1-0.18,1+0.18);

      //adding a wiggle to the expo parameters (Gaussian, Uniform or no wiggle)
      //I use uniform usuak=lly
      //double expo0_pse=gRandom->Uniform(expo0-expo0err,expo0+expo0err);
      //double expo1_pse=gRandom->Uniform(expo1-expo1err,expo1+expo1err);
      //double expo0_pse=gRandom->Uniform(expo0,expo0);
      //double expo1_pse=gRandom->Uniform(expo1,expo1);
      double expo0_pse=gRandom->Gaus(expo0,expo0err);
      double expo1_pse=gRandom->Gaus(expo1,expo1err);

      //modelling background with

      Double_t funPars_mcexpo[9]={expo0_pse,expo1_pse};
      func_mcexpo.SetParameters(funPars_mcexpo);

      outputfile->cd();
      //func_mcexpo.Write();
      //defining histogram to fill with data background
      char hisID[50], title[50];
      sprintf(hisID,"hBmass%03i",p); 
      TH1D *hBmass = new TH1D(hisID,title,1000,0,1000); 
      //Throw the PSE
      for(int i=0; i<hBmass->GetNbinsX(); i++)
      {//function values are per 10 GeV bin -> we have to divide by 10 when filling in 1 GeV steps
        float value=gRandom->Poisson(func_mcexpo.Eval(i+0.5)/10);
        hBmass->SetBinContent(i, value);
      }
      //changing some properties of the PSE histogram     
      hBmass->Sumw2(); 
      hBmass->SetLineColor(30);
      hBmass->SetLineWidth(2);
      hBmass->Rebin(10);
      hBmass->GetXaxis()->SetRangeUser(160,800);

      //Fittng functions with expo
      TF1 *fit_dataexpo=new TF1("fit_dataexpo", "expo(0)+gaus(2)",200, 800);    
      //fixin the expo values to the nominal values without wiggle
      float err=0.000000001;
      //comment in the next lines for floating bg parameters
      //	fit_dataexpo->SetParameter(0,expo0);	
      //fit_dataexpo->SetParLimits(0,expo0-expo0err,expo0+expo0err);
      //fit_dataexpo->SetParameter(1,expo1);	
      //fit_dataexpo->SetParLimits(1, expo1-expo1err,expo1+expo1err); 
      //comment in for fixed bg parameters
      fit_dataexpo->SetParameter(0,expo0);	
      fit_dataexpo->SetParLimits(0,expo0-err,expo0+err);
      fit_dataexpo->SetParameter(1,expo1);	
      fit_dataexpo->SetParLimits(1, expo1-err,expo1+err); 
      //initizalising the Gaussian part of the fit, allow negativ amplitude and fix mass	
      float amp=15;
      fit_dataexpo->SetParameters(2,amp);
      fit_dataexpo->SetParLimits(2,-1000,1000);
      float mass=380;
      fit_dataexpo->SetParameters(3,Mjjj);
      //mass fixed
      fit_dataexpo->SetParLimits(3,Mjjj-0.1,Mjjj+0.1);
      //fit_dataexpo->SetParLimits(3,200,600);
      //width of the initialize depending on Mjjj
      float mwidth = -999;
      if (Mjjj <= 250){
        mwidth=12.5;
      }
      if ((Mjjj > 250) && (Mjjj < 350)){
        mwidth=15;
      }
      if (Mjjj >= 350){
        mwidth=22.5;
      }

      fit_dataexpo->SetParameters(4,mwidth);

      //set limits for the width depending on Mjjj

      float lowwidth  = -999;
      float highwidth = -999;
      if (Mjjj <= 250){
        lowwidth=10;
        highwidth=15;
      }
      if ((Mjjj > 250) && (Mjjj < 350)){
        lowwidth=10 + 10 * (Mjjj - 250.) / 100.;
        highwidth=15 + 10 * (Mjjj - 250.) / 100.;
      }
      if (Mjjj >= 350){
        lowwidth=20;
        highwidth=25;
      }


      fit_dataexpo->SetParLimits(4,lowwidth,highwidth);
      //fit binned likelihood
      hBmass->Fit(fit_dataexpo,"QL","",170,800);

      outputfile->cd();
      //hBmass->Write();
      //check if the fit converged, if not redo
      if (!gMinuit->fCstatu.Contains("CONVERGED"))
      {
        outputfile->cd("not_converged");
        hBmass->Write();
        hBmass->Fit(fit_dataexpo,"QL","",170,800);
        hBmass->Write();
      }

      amp=fit_dataexpo->GetParameter(2);
      //getting the Gauss integral back (factor of sqrt(2 pi)=2.50663
      float sig_evt=amp*fit_dataexpo->GetParameter(4)*2.50663/10;
      width->Fill(fit_dataexpo->GetParameter(4));



      //	cout<<sig_evt<<"  "<<<<endl;
      //bringing the smearfactor again into the acceptance
      acceptance=acceptance*smearfactor;
      //turn n_evts into a cross section
      vec_numobs_1sigma.push_back(sig_evt/(lumi*acceptance));
      //write out histograms that have measured cross section higher than value in data
      if ((sig_evt/(lumi*acceptance)>=38.6))
      {	outputfile->cd();
        hBmass->Write();
        count_pes_pass++;
        // cout<<smearfactor<<endl;
      }
      evt_1sigma->Fill(sig_evt/(lumi*acceptance));
      std::cout<<"NPEs: "<<p<<"  mass:"<<fit_dataexpo->GetParameter(3)<<"  sig_evt: "<<sig_evt<<" acc: "<<acceptance<<" xsec: "<<sig_evt/(lumi*acceptance)<<"  r:"<<(fit_dataexpo->GetParameter(2)/fit_dataexpo->GetParError(2))<<std::endl;

      //bringing the smearfactor

      delete hBmass; delete fit_dataexpo; 
    }
  }
  //sorting the cross sections
  std::sort(vec_numobs_1sigma.begin(), vec_numobs_1sigma.end());
  //find length of the array
  float num=vec_numobs_1sigma.size();
  float fxsec=vec_numobs_1sigma[(int)(nps*0.5-1)];

  //pushing everything back to save it for a graph
  vec_median.push_back(fxsec);
  num_evt_1sigup.push_back(fabs(vec_numobs_1sigma[floor(0.841*nps)]-vec_numobs_1sigma[floor(nps*0.5-1)]));
  num_evt_1sigdown.push_back(fabs(vec_numobs_1sigma[floor(nps*0.5-1)]-vec_numobs_1sigma[ceil(0.159*nps)]));	
  num_evt_2sigup.push_back(fabs(vec_numobs_1sigma[floor(0.977*nps)]-vec_numobs_1sigma[floor(nps*0.5-1)]));
  num_evt_2sigdown.push_back(fabs(vec_numobs_1sigma[floor(nps*0.5-1)]-vec_numobs_1sigma[ceil(0.023*nps-1)]));


  //as well as the median with smearing


  outputfile->cd();
  evt_1sigma->Write();
  width->Write();

  delete evt_1sigma;       
  vmass.push_back(170); 

  //writing out graphs for smear or unsmear
  //  cout<<vec_median.size()<<"   "<<vec_95.size()<<endl;
  TGraphAsymmErrors* graph_median_1sigma = new TGraphAsymmErrors(vec_median.size(),
      &vmass[0], 
      &vec_median[0],0,0,&num_evt_1sigdown[0],&num_evt_1sigup[0]); 
  outputfile->cd();

  graph_median_1sigma->SetTitle("xsec with insec0 of observed events median _1sigma");
  graph_median_1sigma->SetName("median_1sigma");
  graph_median_1sigma->Write();
  float significance=0;

  TH1D *h_significance = new TH1D("significance","signifiance",2,0,2); 
  h_significance->Fill(1.,(float) count_pes_pass);
  h_significance->Fill(0.,(float) count_pes);
  outputfile->cd();
  h_significance->Write();
  std::cout<<count_pes_pass<<" PES out of "<<count_pes<<" give an excess with r=gauss_amp/gauss_amp_err >= 2.1"<<std::endl;





}







int main (int argc, char* argv[]) {
  CMS_significance();

  return 0;
}
