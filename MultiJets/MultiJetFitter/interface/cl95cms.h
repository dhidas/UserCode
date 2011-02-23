#ifndef GUARD_cl95cms_h
#define GUARD_cl95cms_h

#include <iostream>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TArrow.h"
#include "TCanvas.h"


static Double_t A0, sA, B0, sB, epsilon, MaxSig = 100.;
static Double_t MinLike = 1.e-6, Precision = 1.e-5;
static Int_t N;
static bool lGauss = kFALSE, plot = kTRUE;
static Int_t I = 0;
static Double_t sigma_a = 0., sigma_b = 0., tau_a = 0., tau_b = 0.;

Double_t Likelihood(Double_t *x, Double_t *p);
Double_t Inner(Double_t *x, Double_t *par);
Double_t Outer(Double_t *x, Double_t *p);
Double_t Poisson(Double_t Mu, Int_t n);
Double_t CL95(Double_t ilum, Double_t slum, Double_t eff, Double_t seff, Double_t bck, Double_t sbck, Int_t n, Bool_t gauss = kFALSE, Int_t nuisanceModel = 0);
Double_t CLA(Double_t ilum, Double_t slum, Double_t eff, Double_t seff, Double_t bck, Double_t sbck, Int_t bckint = 0);

#endif
