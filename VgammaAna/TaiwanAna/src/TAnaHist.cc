////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Jul 14 11:18:08 PDT 2010
//
////////////////////////////////////////////////////////////////////

#include "VgammaAna/TaiwanAna/interface/TAnaHist.h"

#include <iostream>



TAnaHist::TAnaHist ()
{
  fRootFile = new TFile("TAnaHist.root", "recreate");
  fHistDirName = ".";
  Init();
}

TAnaHist::TAnaHist (TFile* rootfile, TString name)
{
  fRootFile = rootfile;
  fHistDirName = name;
  Init();
}

TAnaHist::~TAnaHist ()
{
  //fHistDir->Write();
}






void TAnaHist::NewTH1F (TString name, int bins, float min, float max)
{
  fMapTH1F[name] = new TH1F(name, name, bins, min, max);
  fMapTH1F[name]->SetDirectory(fHistDir);
  return;
}

void TAnaHist::NewTH1F (TString name, TH1F *InHist)
{
  fMapTH1F[name] = (TH1F*) InHist->Clone(name);
  fMapTH1F[name]->SetDirectory(fHistDir);
  return;
}

void TAnaHist::NewTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax)
{
  fMapTH2D[name] = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
  fMapTH2D[name]->SetDirectory(fHistDir);
  return;
}

void TAnaHist::NewTH3D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax)
{
  fMapTH3D[name] = new TH3D(name, name, xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax);
  fMapTH3D[name]->SetDirectory(fHistDir);
  return;
}

void TAnaHist::AddTH1F (TString name, TH1F *InHist)
{
  if (fMapTH1F.find(name) == fMapTH1F.end()) {
    NewTH1F(name, InHist);
  } else {
    fMapTH1F[name]->Add(InHist);
  }
  return;
}

int TAnaHist::FillTH1F (TString name, float value)
{
  return fMapTH1F[name]->Fill(value);
}

int TAnaHist::FillTH1F (TString name, int bins, float min, float max, float value, float weight)
{
  if (fMapTH1F.find(name) == fMapTH1F.end()) {
    NewTH1F(name, bins, min, max);
  }
  return fMapTH1F[name]->Fill(value, weight);
}

int TAnaHist::FillTH1F (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight)
{
  if (fMapTH1F.find(name) == fMapTH1F.end()) {
    NewTH1F(name, bins, min, max);
    fMapTH1F[name]->SetTitle(title);
    fMapTH1F[name]->SetXTitle(xtitle);
    fMapTH1F[name]->SetYTitle(ytitle);
  }
  return fMapTH1F[name]->Fill(value, weight);
}

TH1F* TAnaHist::GetTH1F(TString const name)
{
  if (fMapTH1F.find(name) == fMapTH1F.end()) {
    NewTH1F(name, 100, -1000, 1000);
    std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
  }
  return fMapTH1F[name];
}






void TAnaHist::NewTH1D (TString name, int bins, float min, float max)
{
  fMapTH1D[name] = new TH1D(name, name, bins, min, max);
  fMapTH1D[name]->SetDirectory(fHistDir);
  fMapTH1D[name]->Sumw2();
  return;
}

void TAnaHist::NewTH1D (TString name, TH1D *InHist)
{
  fMapTH1D[name] = (TH1D*) InHist->Clone(name);
  fMapTH1D[name]->SetDirectory(fHistDir);
  fMapTH1D[name]->Sumw2();
  return;
}


void TAnaHist::AddTH1D (TString name, TH1D *InHist)
{
  if (fMapTH1D.find(name) == fMapTH1D.end()) {
    NewTH1D(name, InHist);
  } else {
    fMapTH1D[name]->Add(InHist);
  }
  return;
}

int TAnaHist::FillTH1D (TString name, float value)
{
  return fMapTH1D[name]->Fill(value);
}

int TAnaHist::FillTH1D (TString name, int bins, float min, float max, float value, float weight)
{
  if (fMapTH1D.find(name) == fMapTH1D.end()) {
    NewTH1D(name, bins, min, max);
  }
  return fMapTH1D[name]->Fill(value, weight);
}

int TAnaHist::FillTH1D (TString name, TString title, TString xtitle, TString ytitle, int bins, float min, float max, float value, float weight)
{
  if (fMapTH1D.find(name) == fMapTH1D.end()) {
    NewTH1D(name, bins, min, max);
    fMapTH1D[name]->SetTitle(title);
    fMapTH1D[name]->SetXTitle(xtitle);
    fMapTH1D[name]->SetYTitle(ytitle);
  }
  return fMapTH1D[name]->Fill(value, weight);
}

TH1D* TAnaHist::GetTH1D(TString const name)
{
  if (fMapTH1D.find(name) == fMapTH1D.end()) {
    NewTH1D(name, 100, -1000, 1000);
    std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
  }
  return fMapTH1D[name];
}








int TAnaHist::FillTH2D (TString name, float xvalue, float yvalue)
{
  return fMapTH2D[name]->Fill(xvalue, yvalue);
}

int TAnaHist::FillTH2D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight)
{
  if (fMapTH2D.find(name) == fMapTH2D.end()) {
    NewTH2D(name, xbins, xmin, xmax, ybins, ymin, ymax);
  }
  return fMapTH2D[name]->Fill(xvalue, yvalue, weight);
}

int TAnaHist::FillTH2D (TString name, TString title, TString xtitle, TString ytitle, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, float xvalue, float yvalue, float weight)
{
  if (fMapTH2D.find(name) == fMapTH2D.end()) {
    NewTH2D(name, xbins, xmin, xmax, ybins, ymin, ymax);
    fMapTH2D[name]->SetTitle(title);
    fMapTH2D[name]->SetXTitle(xtitle);
    fMapTH2D[name]->SetYTitle(ytitle);
  }
  return fMapTH2D[name]->Fill(xvalue, yvalue, weight);
}


TH2D* TAnaHist::GetTH2D(TString const name)
{
  if (fMapTH2D.find(name) == fMapTH2D.end()) {
    NewTH2D(name, 100, -1000, 1000, 100, -1000, 1000);
    std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
  }
  return fMapTH2D[name];
}


int TAnaHist::FillTH3D (TString name, float xvalue, float yvalue, float zvalue)
{
  return fMapTH3D[name]->Fill(xvalue, yvalue, zvalue);
}

int TAnaHist::FillTH3D (TString name, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax, float xvalue, float yvalue, float zvalue, float weight)
{
  if (fMapTH3D.find(name) == fMapTH3D.end()) {
    NewTH3D(name, xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax);
  }
  return fMapTH3D[name]->Fill(xvalue, yvalue, zvalue, weight);
}

int TAnaHist::FillTH3D (TString name, TString title, TString xtitle, TString ytitle, TString ztitle, int xbins, float xmin, float xmax, int ybins, float ymin, float ymax, int zbins, float zmin, float zmax, float xvalue, float yvalue, float zvalue, float weight)
{
  if (fMapTH3D.find(name) == fMapTH3D.end()) {
    NewTH3D(name, xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax);
    fMapTH3D[name]->SetTitle(title);
    fMapTH3D[name]->SetXTitle(xtitle);
    fMapTH3D[name]->SetYTitle(ytitle);
    fMapTH3D[name]->SetZTitle(ztitle);
  }
  return fMapTH3D[name]->Fill(xvalue, yvalue, zvalue, weight);
}


TH3D* TAnaHist::GetTH3D(TString const name)
{
  if (fMapTH3D.find(name) == fMapTH3D.end()) {
    NewTH3D(name, 100, -1000, 1000, 100, -1000, 1000, 100, -1000, 1000);
    std::cerr << "WARNING: TAnaHist possible error in histogram name " << name << std::endl;
  }
  return fMapTH3D[name];
}




void TAnaHist::Init ()
{
  if (fHistDirName == ".") {
    fHistDir = (TDirectory*) fRootFile;
  } else {
    fHistDir = fRootFile->mkdir(fHistDirName);
  }
  return;
}
