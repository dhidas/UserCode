//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  6 21:57:28 2009 by ROOT version 5.22/00
// from TTree t/out tree
// found on file: sr90-s21-run16.root
//////////////////////////////////////////////////////////

#ifndef study_h
#define study_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "dead_noisy.h"

class study {
public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain
	TFile *out;
	
	// Declaration of leaf types
	Int_t           tick;
	Int_t           container[204];
	Int_t           vcal;
	Int_t           range;
	Int_t           levelinfo_cut[5];
	Int_t           ub_threshold;
	Int_t           levelinfo[5];
	Int_t           ub[3];
	Int_t           ubROC1;
	Int_t           ub_tr1;
	Int_t           ub_tr2;
	Int_t           tbm_status[4];
	Int_t           nhits;
	Int_t           ic[100];
	Int_t           ir[100];
	Int_t           adc[100];
	
	float min_ic;
	float max_ic;
	int n_ic;
	
	float min_ir;
	float max_ir;
	int n_ir;
	
	
	TH1F *RawADC_h;
	TH1F * Electrons_h;
	TH1F * ClusterE;
	TH1F *ClusterADC;
	TH1F *MaxADC;
	TH1F *MaxE;
	TH3F * contained3d;
	

	TH2F *clusterET1;
	TH2F *clusterET2;
	TH2F *clusterET3;
	
	//hist for bugs finding
	TH1F* ub1_h;
	TH1F* ub2_h;
	TH1F* ub3_h;
	TH1F* ubROC1_h;
	TH1F* ub_tr1_h;
	TH1F* ub_tr2_h;
	TH1F* ubdiff;
	TH1F* tbmdiff;
	TH1F* tbm1;// =new TH1F("tbm1", "tbm1",n_tbm,min_tbm,max_tbm);
	TH1F* tbm2;// =new TH1F("tbm2", "tbm2",n_tbm,min_tbm,max_tbm);
	TH1F* tbm3;// =new TH1F("tbm3", "tbm3",n_tbm,min_tbm,max_tbm);
	TH1F* tbm4;// =new TH1F("tbm4", "tbm4",n_tbm,min_tbm,max_tbm);
	
	TH2F* ub1t_h;// =new TH2F("ub1t", "ultrablack one vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	TH2F* ub2t_h;// =new TH2F("ub2t", "ultrablack two vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	TH2F* ub3t_h;// =new TH2F("ub3t", "ultrablack three vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	TH2F* ubROC1t_h;// =new TH2F("ubROC1t", "ROC1 ultrablack vs time", n_time, min_time, max_time,n_ubroc, min_ubroc, max_ubroc);
	TH2F* ub_tr1t_h;// =new TH2F("ub_tr1t", "1^{st} trailer ultrablack vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	TH2F* ub_tr2t_h;// =new TH2F("ub_tr2t", "2^{nd} trailer ultrablack vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	TH2F* ubdifft;// =new TH2F("ubdifft", "difference in ultrablack level b/t 1st and 2nd ub vs time", n_time, min_time, max_time,100,-50,50);
	
	TH2F* levelinfot_h;// =new TH2F("levelinfo_t", "level info", n_time, min_time, max_time,n_level, min_level, max_level);
	
	TH1F* levelinfo_h;// =new TH1F("levelinfo", "level info",n_level, min_level, max_level);
	TH1F* levelinfo_cut_h;// =new TH1F("levelinfo_cut", "cut on level info",n_level, min_level, max_level);
	
	TH1F* c0;// =new TH1F("c0", "level info: c0",n_level, min_level, max_level);
	TH1F* c1;// =new TH1F("c1", "level info: c1",n_level, min_level, max_level);
	TH1F* r0;// =new TH1F("r0", "level info: r0",n_level, min_level, max_level);
	TH1F* r1;// =new TH1F("r1", "level info: r1",n_level, min_level, max_level);
	TH1F* r2;// =new TH1F("r2", "level info: r2",n_level, min_level, max_level);
	
	//Occupancy map
	TH2F* h_number_of_hits_map;// =new TH2F("h_number_of_hits_map", "num hits per pixel", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);

	TH2F * h_number_of_hits_map_efficiency;
	TH1F *h_number_of_hits_map_efficiency1d;
	TH1F *h_number_of_hits_map_efficiency1d_all;
	
	TH2F * h_energy_hits_map_efficiency;
	TH1F *h_energy_hits_map_efficiency1d;
	TH1F *h_energy_hits_map_efficiency1d_all;
	TH2F *markloweff;
	TH2F *markhigheff;
	
	TH2F* occupancy_one;// =new TH2F("occupancy_one", "occumapncy of one pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	TH2F* occupancy_two;// =new TH2F("occupancy_two", "occumapncy of 2 pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	TH2F* occupancy_three;// =new TH2F("occupancy_three", "occumapncy of >=3 pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	TH2F* occupancy_ratio;
	TH1F* h_occupancy_ratio;
	
	TH2F* occupancy_dead ;//=new TH2F("occupancy_dead", "occumapncy of one pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	TH2F* occupancy_noisy;// =new TH2F("occupancy_noisy", "occumapncy of one pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	TH2F* energetically_dead;
	TH2F* energetically_noisy;
	TH2F* occupancy_abs_dead;
	TH2F* energetically_abs_dead;

	
	// Number of hits: per event, per cluster, per cluster fiducial
	TH1F* eventhits ;//           = new TH1F("eventhits", "number of pixel hits per event",11,-0.5,10.5);
	TH1F* clusterhits;//          = new TH1F("clusterhits", "number of pixel hits per cluster",11,-0.5,10.5);
	TH1F* clusterhits_fiducial;// = new TH1F("clusterhits_fiducial", "number of pixel hits per cluster in fiducial region",11,-0.5,10.5);
	
	TH1F* occupancy_all;//            = new TH1F("occupancy_all", "hit occupancy",20000,0,20000);
	
	//additional vcal to e hist , _one is for event with only one pixel hit, _two is for 2 pixelhits, _three is for 3+
	TH1F* adcelectron_one;//    = new TH1F("adcelectron_one",    "Pulse height 1 pixelhit events", n_e, min_e, max_e);
	TH1F* adcelectron_two ;//   = new TH1F("adcelectron_two",    "Pulse height for 2 pixelhit events", n_e, min_e, max_e);
	TH1F* adcelectron_three ;// = new TH1F("adcelectron_three",  "Pulse height for 3+ pixelhit events", n_e, min_e, max_e);
	TH1F* adcelectron_four;//  = new TH1F("adcelectron_four",  "Pulse height for 3+ pixelhit events", n_e, min_e, max_e);
	TH1F* adcelectron_max ;//= new TH1F("h_max_pulse_height", "Pulse heghts of pixels with max pulse height", n_e, min_e, max_e);

	TH1F* adcelectron_20_69;//    = new TH1F("adcelectron_20_69",    "Pulse height for pixel (20,69)", n_e, min_e, max_e);
	
	TH2F* adcelectron_time;//    = new TH2F("adcelectron_vs_time","Pulse height 2 pixelhit events vs time", n_time, min_time, max_time,  n_e, min_e, max_e);

	
	// Energy deposition map
	TH2F* h_energy_hits_map;// =new TH2F("h_number_of_hits_map", "num hits per pixel", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);

	
	TH2F *h2;
	TH1F *residuals_h;
	TH1F *saturated_vcal_h;
	TH1F *minimal_vcal_h;
	TH1F *chi2;
	TH1F * chi2ndf;
	TH2F *chi2_2d;
	TH2F * chi2ndf_2d;
	TF1 *f1;
	TF1 *f2;
	TGraphErrors *p;
	TGraphErrors *pe;
	
	TH1F *chi2s;
	TH1F * chi2ndfs;
	TH2F *chi2_2ds;
	TH2F * chi2ndf_2ds;

	TH1F * fwhm_h;
	TH1F *halfpoint_h;
	TH2F * fwhm2d_h;
	TH2F *halfpoint2d_h;
	TH2F *anompix;
	
	TH3F *pixel_energy;
	// List of branches
   TBranch        *b_t;   //!

	study(TTree *tree=0, char* out_name="");
	virtual ~study();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     analyze();
	virtual void     scurve();
	virtual void     fill_hist   (int ic,int ir);
	virtual void     fill_scurve (int ic,int ir);
	virtual void     calibrate();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
void processcal(ifstream& file);

  int run (int argc, char** argv);
};

#endif

#ifdef study_cxx
study::study(TTree *tree, char* out_name)
{		
	out = new TFile(out_name,"recreate");
	if (tree == 0) {
		printf("oops, can't find the tree \"t\"\n");
	}
	
	Init(tree);
}

study::~study()
{
   if (!fChain) return;
	delete fChain->GetCurrentFile();
	out->Write();
	out->Close();
}

Int_t study::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t study::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void study::Init(TTree *tree)
{
	int min_adc = 0;
	int max_adc = 4000;
	int n_adc = 100;
	
	int min_e = 0;
	int max_e = 50000;
	int n_e = 100;
	 
	int min_level = 1800;
	int max_level = 3500;
	int n_level = 100;
	
	int min_tbm = 1800;
	int max_tbm = 2500;
	int n_tbm = 300;
	
	int min_ub = 1200;
	int max_ub = 1400;
	int n_ub = 100;
	
	int min_ubroc = 1200;
	int max_ubroc = 1500;
	int n_ubroc = 100;
	
	min_ic = 12;
	max_ic = 42;
	n_ic = 30;
	
	min_ir = 40;
	max_ir = 83;
	n_ir =  43;
	
	int min_time = 0;
	int max_time = 86400;
	int n_time = 20*60*60;// i.e. one bin is 60 seconds
	
	int min_index = 0;
	int max_index = 1000;
	int n_index = 1000;
	
	
	RawADC_h = new TH1F("RawADC_h","pulse height in Raw ADC", 500,2000,4500);
	Electrons_h = new TH1F("Electrons_h", "Pulse height in Electrons",100,0,60000);
	ClusterE = new TH1F("ClusterE", "Pulse height in a 2x2 Clusters", 100,0,60000);
	ClusterADC = new TH1F("ClusterADC","pulse height in Raw ADC",500,2000,4500);
	MaxADC = new TH1F("MaxADC","pulse height in Raw ADC max hit", 500,2000,4500);
	MaxE = new TH1F("MaxE","pulse height E od max hit",100,0,60000);
	contained3d = new TH3F("container3d","",100, 0, 100, 28, 13, 40, 41, 40, 80);
	//hist for bugs finding
	ub1_h =new TH1F("ub1", "ultrablack one", n_ub, min_ub, max_ub);
	ub2_h =new TH1F("ub2", "ultrablack two", n_ub, min_ub, max_ub);
	ub3_h =new TH1F("ub3", "ultrablack three", n_ub, min_ub, max_ub);
	ubROC1_h =new TH1F("ubROC1", "ROC1 ultrablack",n_ubroc, min_ubroc, max_ubroc);
	ub_tr1_h =new TH1F("ub_tr1", "1^{st} trailer ultrablack", n_ub, min_ub, max_ub);
	ub_tr2_h =new TH1F("ub_tr2", "2^{nd} trailer ultrablack", n_ub, min_ub, max_ub);
	ubdiff =new TH1F("ubdiff", "difference in ultrablack level b/t 1st and 2nd ub",100,-50,50);
	tbmdiff =new TH1F("tbmdiff", "difference in tbm status",100,-50,50);
	tbm1 =new TH1F("tbm1", "tbm1",n_tbm,min_tbm,max_tbm);
	tbm2 =new TH1F("tbm2", "tbm2",n_tbm,min_tbm,max_tbm);
	tbm3 =new TH1F("tbm3", "tbm3",n_tbm,min_tbm,max_tbm);
	tbm4 =new TH1F("tbm4", "tbm4",n_tbm,min_tbm,max_tbm);
	
	ub1t_h =new TH2F("ub1t", "ultrablack one vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	ub2t_h =new TH2F("ub2t", "ultrablack two vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	ub3t_h =new TH2F("ub3t", "ultrablack three vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	ubROC1t_h =new TH2F("ubROC1t", "ROC1 ultrablack vs time", n_time, min_time, max_time,n_ubroc, min_ubroc, max_ubroc);
	ub_tr1t_h =new TH2F("ub_tr1t", "1^{st} trailer ultrablack vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	ub_tr2t_h =new TH2F("ub_tr2t", "2^{nd} trailer ultrablack vs time", n_time, min_time, max_time, n_ub, min_ub, max_ub);
	ubdifft =new TH2F("ubdifft", "difference in ultrablack level b/t 1st and 2nd ub vs time", n_time, min_time, max_time,100,-50,50);
	
	levelinfot_h =new TH2F("levelinfo_t", "level info", n_time, min_time, max_time,n_level, min_level, max_level);
	
	levelinfo_h =new TH1F("levelinfo", "level info",n_level, min_level, max_level);
	levelinfo_cut_h =new TH1F("levelinfo_cut", "cut on level info",n_level, min_level, max_level);
	
	c0 =new TH1F("c0", "level info: c0",n_level, min_level, max_level);
	c1 =new TH1F("c1", "level info: c1",n_level, min_level, max_level);
	r0 =new TH1F("r0", "level info: r0",n_level, min_level, max_level);
	r1 =new TH1F("r1", "level info: r1",n_level, min_level, max_level);
	r2 =new TH1F("r2", "level info: r2",n_level, min_level, max_level);
	
	//Occupancy map
	h_number_of_hits_map =new TH2F("h_number_of_hits_map", "num hits per pixel", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	
	// Derive efficiencies
	h_number_of_hits_map_efficiency =  new TH2F("h_number_of_hits_map_efficiency", "", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	h_number_of_hits_map_efficiency1d =  new TH1F("h_number_of_hits_map_efficiency_index","", 1000, 0, 1000);
	h_number_of_hits_map_efficiency1d_all =  new TH1F("h_number_of_hits_map_efficiency_all","", 100, -0.1, 3);

	h_energy_hits_map_efficiency =  new TH2F("h_energy_hits_map_efficiency", "", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	h_energy_hits_map_efficiency1d =  new TH1F("h_energy_hits_map_efficiency_index","", 1000, 0, 1000);
	h_energy_hits_map_efficiency1d_all =  new TH1F("h_energy_hits_map_efficiency_all","", 100, -0.1, 3);
		
      markloweff =  new TH2F("markloweff", "", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
      markhigheff =  new TH2F("markhigheff", "", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	
	occupancy_one =new TH2F("occupancy_one", "occumapncy of one pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	
	occupancy_two =new TH2F("occupancy_two", "occumapncy of 2 pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	occupancy_three =new TH2F("occupancy_three", "occumapncy of >=3 pixel clusters", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	h_occupancy_ratio = new TH1F("h_occupancy_ratio", "occupancy ratio", 60, 0, 3);
	
	occupancy_dead =new TH2F("occupancy_dead", "Pixels with anomolously low occupancy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	occupancy_noisy =new TH2F("occupancy_noisy", "Pixels with anomolously high occupancy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	energetically_dead =new TH2F("energetically_dead", "Pixels with anomolously low energy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	energetically_noisy =new TH2F("energetically_noisy", "Pixels with anomolously high energy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	occupancy_abs_dead =new TH2F("occupancy_abs_dead", "Pixels with anomolously low occupancy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	energetically_abs_dead =new TH2F("energetically_abs_dead", "Pixels with anomolously low energy", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);

	
	// Number of hits: per event, per cluster, per cluster fiducial
	eventhits            = new TH1F("eventhits", "number of pixel hits per event",11,-0.5,10.5);
	clusterhits          = new TH1F("clusterhits", "number of pixel hits per cluster",11,-0.5,10.5);
	clusterhits_fiducial = new TH1F("clusterhits_fiducial", "number of pixel hits per cluster in fiducial region",11,-0.5,10.5);
	
	occupancy_all            = new TH1F("occupancy_all", "hit occupancy",20000,0,20000);
	
	//additional vcal to e hist , _one is for event with only one pixel hit, _two is for 2 pixelhits, _three is for 3+
	adcelectron_one    = new TH1F("adcelectron_one",    "Pulse height 1 pixelhit events", n_e, min_e, max_e);
	adcelectron_two    = new TH1F("adcelectron_two",    "Pulse height for 2 pixelhit events", n_e, min_e, max_e);
	adcelectron_three  = new TH1F("adcelectron_three",  "Pulse height for 3+ pixelhit events", n_e, min_e, max_e);
	adcelectron_four  = new TH1F("adcelectron_four",  "Pulse height for 3+ pixelhit events", n_e, min_e, max_e);
	adcelectron_max = new TH1F("h_max_pulse_height", "Pulse heghts of pixels with max pulse height", n_e, min_e, max_e);
	
//	adcelectron_one->Sumw2();
//	adcelectron_two->Sumw2();
//	adcelectron_three->Sumw2();
	adcelectron_20_69    = new TH1F("adcelectron_20_69",    "Pulse height for pixel (20,69)", n_e, min_e, max_e);
	
	adcelectron_time    = new TH2F("adcelectron_vs_time","Pulse height 1 pixelhit events vs time", n_time, min_time, max_time,  n_e, min_e, max_e);
	
	clusterET1   = new TH2F("clusterET1","Cluster Energy 1 pixel vs time", n_time, min_time, max_time,  n_e, min_e, max_e);
	clusterET1->SetBit(TH2::kCanRebin);
	clusterET2   = new TH2F("clusterET2","Cluster Energy 2 pixels vs time", n_time, min_time, max_time,  n_e, min_e, max_e);
	clusterET2->SetBit(TH2::kCanRebin);
	clusterET3   = new TH2F("clusterET3","Cluster Energy 2 or more pixels vs time", n_time, min_time, max_time,  n_e, min_e, max_e);
	clusterET3->SetBit(TH2::kCanRebin);


	h_energy_hits_map =new TH2F("h_energy_hits_map", "energy deposited in pixels", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);

	saturated_vcal_h = new TH1F("saturated_vcal", "",1000, 0, 65*1000);
	minimal_vcal_h = new TH1F("minimal_vcal", "", 200, 0,65*200);
	chi2 = new TH1F("chi2","chi2 of all fits", 300, 0, 3000);
	chi2ndf = new TH1F("chi2ndf","chi2/n.d.f of all fits", 100, 0, 100);
	chi2_2d = new TH2F("chi2_2d","fit chi2 map", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	chi2ndf_2d = new TH2F("chi2ndf_2d","fit chi2/n.d.f. map", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);

	//fits for scurves
	chi2s = new TH1F("chi2s","chi2s of all fits", 300, 0, 1);
	chi2s->SetBit(TH1::kCanRebin);
	chi2ndfs = new TH1F("chi2ndfs","chi2/n.d.f of all fits", 100, 0, 1);
	chi2ndfs->SetBit(TH1::kCanRebin);
	chi2_2ds = new TH2F("chi2_2ds","fit chi2 map", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	chi2_2ds->SetBit(TH1::kCanRebin);
	chi2ndf_2ds = new TH2F("chi2ndf_2ds","fit chi2/n.d.f. map", n_ic, min_ic, max_ic,n_ir, min_ir, max_ir);
	chi2_2ds->SetBit(TH1::kCanRebin);
	
	
	fwhm_h=new TH1F("fwhm", "FWDM", 50, 0, 1000);
//	fwhm_h->SetBit(TH1::kCanRebin);
	halfpoint_h=new TH1F("halfpoint", "halfpoint", 100, 0, 10000);
//	halfpoint_h->SetBit(TH1::kCanRebin);
	fwhm2d_h=new TH2F("fwhm2d", "FWDM", n_ic, min_ic, max_ic, n_ir, min_ir, max_ir); 
//	fwhm2d_h->SetBit(TH1::kCanRebin);
	halfpoint2d_h=new TH2F("halfpoint2d", "halfpoint", n_ic, min_ic, max_ic, n_ir, min_ir, max_ir);
//	halfpoint2d_h->SetBit(TH1::kCanRebin);
	anompix= new TH2F("anompix","Pixels With Outlying Fits", n_ic, 0.5, 0.6,n_ir, 0.5, 0.6);
//	anompix->SetBit(TH1::kCanRebin);
	

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("t", &tick, &b_t);
   Notify();
}

Bool_t study::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void study::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t study::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef study_cxx
