class TH1;

class LangauFitter{
 private:
  Double_t fp[4], fpe[4];
  TF1 *fitsnr; 
  
 public:
  //  LangauFitter();
  void Fit(TH1* h, Double_t limitLow = 0., Double_t limitHigh = 0.);
  static Double_t langaufun(Double_t *x, Double_t *par);
  Double_t GetParameter(Int_t i){return fp[i];};
  Double_t GetError(Int_t i){return fpe[i];};
  TF1 *GetFitFct(){return fitsnr;};

 private:
  TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
  Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
};
