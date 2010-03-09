////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on:  2007ish
//
////////////////////////////////////////////////////////////////////

#ifndef GUARD_TAnaHist_h
#define GUARD_TAnaHist_h

#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"


class TAnaHist
{
  public:
    TAnaHist () {
      fRootFile = new TFile("TAnaHist.root", "recreate");
      fHistDirName = ".";
      Init();
    };

    TAnaHist (TFile* rootfile, TString name) {
      fRootFile = rootfile;
      fHistDirName = name;
      Init();
    };

    ~TAnaHist () {
      //fHistDir->Write();
    }







    void NewTH1F (TString name, int bins, float min, float max) {
      fMapTH1F[name] = new TH1F(name, name, bins, min, max);
      fMapTH1F[name]->SetDirectory(fHistDir);
      return;
    }

    void NewTH1F (TString name, TH1F *InHist) {
      fMapTH1F[name] = (TH1F*) InHist->Clone(name);
      fMapTH1F[name]->SetDirectory(fHistDir);
      return;
    }

    void NewTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax) {
      fMapTH2D[name] = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
      fMapTH2D[name]->SetDirectory(fHistDir);
      return;
    }

    void AddTH1F (TString name, TH1F *InHist) {
      if (fMapTH1F.find(name) == fMapTH1F.end()) {
        NewTH1F(name, InHist);
      } else {
        fMapTH1F[name]->Add(InHist);
      }
      return;
    }

    int FillTH1F (TString name, float value) {
      return fMapTH1F[name]->Fill(value);
    }

    int FillTH1F (TString name, int bins, float min, float max, float value, float weight=1) {
      if (fMapTH1F.find(name) == fMapTH1F.end()) {
        NewTH1F(name, bins, min, max);
      }
      return fMapTH1F[name]->Fill(value, weight);
    }

    int FillTH1F (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight=1) {
      if (fMapTH1F.find(name) == fMapTH1F.end()) {
        NewTH1F(name, bins, min, max);
        fMapTH1F[name]->SetTitle(title);
        fMapTH1F[name]->SetXTitle(xtitle);
        fMapTH1F[name]->SetYTitle(ytitle);
      }
      return fMapTH1F[name]->Fill(value, weight);
    }

    TH1F* GetTH1F(TString const name) {
      if (fMapTH1F.find(name) == fMapTH1F.end()) {
        NewTH1F(name, 100, -1000, 1000);
        std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
      }
      return fMapTH1F[name];
    }






    void NewTH1D (TString name, int bins, float min, float max) {
      fMapTH1D[name] = new TH1D(name, name, bins, min, max);
      fMapTH1D[name]->SetDirectory(fHistDir);
      fMapTH1D[name]->Sumw2();
      return;
    }

    void NewTH1D (TString name, TH1D *InHist) {
      fMapTH1D[name] = (TH1D*) InHist->Clone(name);
      fMapTH1D[name]->SetDirectory(fHistDir);
      fMapTH1D[name]->Sumw2();
      return;
    }


    void AddTH1D (TString name, TH1D *InHist) {
      if (fMapTH1D.find(name) == fMapTH1D.end()) {
        NewTH1D(name, InHist);
      } else {
        fMapTH1D[name]->Add(InHist);
      }
      return;
    }

    int FillTH1D (TString name, float value) {
      return fMapTH1D[name]->Fill(value);
    }

    int FillTH1D (TString name, int bins, float min, float max, float value, float weight=1) {
      if (fMapTH1D.find(name) == fMapTH1D.end()) {
        NewTH1D(name, bins, min, max);
      }
      return fMapTH1D[name]->Fill(value, weight);
    }

    int FillTH1D (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight=1) {
      if (fMapTH1D.find(name) == fMapTH1D.end()) {
        NewTH1D(name, bins, min, max);
        fMapTH1D[name]->SetTitle(title);
        fMapTH1D[name]->SetXTitle(xtitle);
        fMapTH1D[name]->SetYTitle(ytitle);
      }
      return fMapTH1D[name]->Fill(value, weight);
    }

    TH1D* GetTH1D(TString const name) {
      if (fMapTH1D.find(name) == fMapTH1D.end()) {
        NewTH1D(name, 100, -1000, 1000);
        std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
      }
      return fMapTH1D[name];
    }








    int FillTH2D (TString name, float xvalue, float yvalue) {
      return fMapTH2D[name]->Fill(xvalue, yvalue);
    }

    int FillTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight=1) {
      if (fMapTH2D.find(name) == fMapTH2D.end()) {
        NewTH2D(name, xbins, xmin, xmax, ybins, ymin, ymax);
      }
      return fMapTH2D[name]->Fill(xvalue, yvalue, weight);
    }
    int FillTH2D (TString name, TString title, TString xtitle, TString ytitle, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight=1) {
      if (fMapTH2D.find(name) == fMapTH2D.end()) {
        NewTH2D(name, xbins, xmin, xmax, ybins, ymin, ymax);
        fMapTH2D[name]->SetTitle(title);
        fMapTH2D[name]->SetXTitle(xtitle);
        fMapTH2D[name]->SetYTitle(ytitle);
      }
      return fMapTH2D[name]->Fill(xvalue, yvalue, weight);
    }


    TH2D* GetTH2D(TString const name) {
      if (fMapTH2D.find(name) == fMapTH2D.end()) {
        NewTH2D(name, 100, -1000, 1000, 100, -1000, 1000);
        std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
      }
      return fMapTH2D[name];
    }


  private:
    TFile* fRootFile;
    TDirectory* fHistDir;
    TString fHistDirName;

    std::map<TString, TH1F*> fMapTH1F;
    std::map<TString, TH1D*> fMapTH1D;
    std::map<TString, TH2D*> fMapTH2D;

    void Init () {
      if (fHistDirName == ".") {
        fHistDir = (TDirectory*) fRootFile;
      } else {
        fHistDir = fRootFile->mkdir(fHistDirName);
      }
      return;
    }


};








#endif
