////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Jul 14 11:18:08 PDT 2010
//
////////////////////////////////////////////////////////////////////

#ifndef GUARD_TAnaHist_h
#define GUARD_TAnaHist_h

#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TString.h"


class TAnaHist
{
  public:
    TAnaHist ();

    TAnaHist (TFile* rootfile, TString name);

    ~TAnaHist ();







    void NewTH1F (TString name, int bins, float min, float max);

    void NewTH1F (TString name, TH1F *InHist);

    void NewTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax);

    void NewTH3D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax);

    void AddTH1F (TString name, TH1F *InHist);

    int FillTH1F (TString name, float value);

    int FillTH1F (TString name, int bins, float min, float max, float value, float weight=1);

    int FillTH1F (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight=1);

    TH1F* GetTH1F(TString const name);






    void NewTH1D (TString name, int bins, float min, float max);

    void NewTH1D (TString name, TH1D *InHist);


    void AddTH1D (TString name, TH1D *InHist);

    int FillTH1D (TString name, float value);

    int FillTH1D (TString name, int bins, float min, float max, float value, float weight=1);

    int FillTH1D (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight=1);

    TH1D* GetTH1D(TString const name);








    int FillTH2D (TString name, float xvalue, float yvalue);

    int FillTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight=1);
    int FillTH2D (TString name, TString title, TString xtitle, TString ytitle, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight=1);


    TH2D* GetTH2D(TString const name);


    int FillTH3D (TString name, float xvalue, float yvalue, float zvalue);

    int FillTH3D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax, float xvalue, float yvalue, float zvalue, float weight=1);
    int FillTH3D (TString name, TString title, TString xtitle, TString ytitle, TString ztitle, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax, float xvalue, float yvalue, float zvalue, float weight=1);


    TH3D* GetTH3D(TString const name);


  private:
    TFile* fRootFile;
    TDirectory* fHistDir;
    TString fHistDirName;

    std::map<TString, TH1F*> fMapTH1F;
    std::map<TString, TH1D*> fMapTH1D;
    std::map<TString, TH2D*> fMapTH2D;
    std::map<TString, TH3D*> fMapTH3D;

    void Init ();


};








#endif
