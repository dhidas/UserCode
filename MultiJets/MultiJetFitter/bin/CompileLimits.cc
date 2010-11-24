////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Oct 14 16:53:54 EDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TStyle.h"


float Median (std::vector<float> const& Vec)
{
  size_t const n = Vec.size();
  float Ret;
  if (n % 2 == 0) {
     Ret = Vec[n/2];
  } else {
    Ret = Vec[n/2 + 1];
  }

  return Ret;
}


float Quantile (std::vector<float> const& Vec, float Q)
{
  // Good enough for large stats

  size_t const n   = (size_t) Vec.size();
  int ntoadd =  (Q * (float) n);
  if (Q > 0) {
    //std::cout << Q << "  " << ntoadd << std::endl;
    return Vec[ntoadd];
  }
  //std::cout << Q << "  " << n << "  " << n-1+ntoadd<< std::endl;
  return Vec[n-1+ntoadd];
}




int CompileLimits (TString const DataFileName, std::vector<TString>&  PEFileNames)
{
  // some variabes we'll read from the file
  int ipe;
  double pvalue, gmean, gsigma;
  TString line;

  std::vector<double> vPValue;
  std::vector<double> vMean;

  // Loop over all input filees
  std::vector<float> Masses;
  std::vector< std::vector<float> > PE;
  std::vector<float> Data;
  size_t NMasses;
  std::string tmp;
  float dtmp;
  for (size_t ifile = 0; ifile != PEFileNames.size(); ++ifile) {

    // open this file
    std::ifstream InFile(PEFileNames[ifile].Data());
    if (!InFile.is_open()) {
      std::cerr << "ERROR: cannot open input file: " << PEFileNames[ifile] << std::endl;
      continue;
    }
    std::cout << "Reading file: " << PEFileNames[ifile] << std::endl;


    if (ifile == 0) {
      float mass;
      std::istringstream InStream;
      std::getline(InFile, tmp);
      InStream.str(tmp);
      while (InStream >> mass) {
        std::cout << "Mass found: " << mass << std::endl;
        Masses.push_back(mass);
      }
      NMasses = Masses.size();
      std::cout << "NMasses: " << NMasses << std::endl;
      PE.resize(NMasses);
      for (size_t i = 0; i != PE.size(); ++i) PE[i].reserve(100000);
    } else {
      std::getline(InFile, tmp);
    }
    while (!InFile.eof()) {
      for (size_t imass = 0; imass < NMasses; ++imass) {
        InFile >> dtmp;
        PE[imass].push_back(dtmp);
      }
    }

  }

  // open the data file and read vals
  std::ifstream InFile(DataFileName.Data());
  if (!InFile.is_open()) {
    std::cerr << "ERROR: cannot open input file: " << DataFileName << std::endl;
  } else {
    std::cout << "Reading Data file: " << DataFileName << std::endl;
  }
  std::getline(InFile, tmp);
  for (size_t imass = 0; imass < NMasses; ++imass) {
    InFile >> dtmp;
    std::cout << "Data limit for " << Masses[imass] << "  = " << dtmp << std::endl;
    Data.push_back(dtmp);
  }

  // sort the values
  for (size_t imass = 0; imass < NMasses; ++imass) {
    std::cout << "Sorting PEs for Mass " << Masses[imass] << std::endl;
    std::sort(PE[imass].begin(), PE[imass].end());
  }

  // Make and fill histogram
  std::vector<TH1F*> Hist;
  Hist.resize(NMasses);
  char BUFF[200];
  for (size_t imass = 0; imass < NMasses; ++imass) {
    std::cout << "Filling histogram for Mass: " << Masses[imass] << std::endl;
    sprintf(BUFF, "LimitsFor_%i", (int) Masses[imass]);
    Hist[imass] = new TH1F(BUFF, BUFF, 100, PE[imass].front(), Quantile(PE[imass], 0.99));
    std::cout << Masses[imass] << "  " << PE[imass].front() << "  " << PE[imass].back() << std::endl;

    for (size_t ipe = 0; ipe != PE[imass].size(); ++ipe) {
      Hist[imass]->Fill(PE[imass][ipe]);
    }
  }

  // make a line for data
  std::vector<TLine*> Line;
  std::vector<TLine*> MLinePE;
  std::vector<TLine*> M2LinePE;
  std::vector<TLine*> M1LinePE;
  std::vector<TLine*> P1LinePE;
  std::vector<TLine*> P2LinePE;
  Line.resize(NMasses);
  MLinePE.resize(NMasses);
  M2LinePE.resize(NMasses);
  M1LinePE.resize(NMasses);
  P1LinePE.resize(NMasses);
  P2LinePE.resize(NMasses);
  for (size_t imass = 0; imass < NMasses; ++imass) {
    Line[imass] = new TLine(Data[imass], 0, Data[imass], Hist[imass]->GetMaximum());
    Line[imass]->SetLineColor(1);
    Line[imass]->SetLineWidth(2);
    float const MedianPE = Median(PE[imass]);
    MLinePE[imass] = new TLine(MedianPE, 0, MedianPE, Hist[imass]->GetMaximum());
    MLinePE[imass]->SetLineColor(2);
    MLinePE[imass]->SetLineWidth(2);
    MLinePE[imass]->SetLineStyle(2);
    float const M2PE = Quantile(PE[imass],  0.95);
    float const M1PE = Quantile(PE[imass],  0.67);
    float const P1PE = Quantile(PE[imass], -0.67);
    float const P2PE = Quantile(PE[imass], -0.95);
    M2LinePE[imass] = new TLine(M2PE, 0, M2PE, Hist[imass]->GetMaximum());
    M2LinePE[imass]->SetLineColor(3);
    M2LinePE[imass]->SetLineWidth(1);
    M2LinePE[imass]->SetLineStyle(3);
    P2LinePE[imass] = new TLine(P2PE, 0, P2PE, Hist[imass]->GetMaximum());
    P2LinePE[imass]->SetLineColor(3);
    P2LinePE[imass]->SetLineWidth(1);
    P2LinePE[imass]->SetLineStyle(3);
    M1LinePE[imass] = new TLine(M1PE, 0, M1PE, Hist[imass]->GetMaximum());
    M1LinePE[imass]->SetLineColor(4);
    M1LinePE[imass]->SetLineWidth(1);
    M1LinePE[imass]->SetLineStyle(2);
    P1LinePE[imass] = new TLine(P1PE, 0, P1PE, Hist[imass]->GetMaximum());
    P1LinePE[imass]->SetLineColor(4);
    P1LinePE[imass]->SetLineWidth(1);
    P1LinePE[imass]->SetLineStyle(2);
  }

  // Make canvas, plot histogram, and save canvas
  gStyle->SetOptStat(1111111111);  
  for (size_t imass = 0; imass < NMasses; ++imass) {
    TCanvas Can;
    Can.cd();
    Hist[imass]->Draw("hist");
    Line[imass]->Draw("same");
    MLinePE[imass]->Draw("same");
    M2LinePE[imass]->Draw("same");
    M1LinePE[imass]->Draw("same");
    P1LinePE[imass]->Draw("same");
    P2LinePE[imass]->Draw("same");
    sprintf(BUFF, "LimitsPlot_%i.eps", (int) Masses[imass]);
    Can.SaveAs(BUFF);
  }


  FILE*  OutFile = fopen("Limits.dat", "w");
  if (OutFile == NULL) {
    std::cerr << "ERROR: cannot open output file" << std::endl;
    exit(1);
  }
  fprintf(OutFile, "Test\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Masses[i]);              fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Quantile(PE[i], -0.95)); fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Quantile(PE[i], -0.68)); fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Median(PE[i]));          fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Quantile(PE[i],  0.68)); fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Quantile(PE[i],  0.95)); fprintf(OutFile, "\n");
  for (size_t i = 0; i != NMasses; ++i)   fprintf(OutFile, "%10.3f ", Data[i]);
  fprintf(OutFile, "\n1");
  fclose(OutFile);


  for (size_t i = 0; i != NMasses; ++i) std::cout << Masses[i] << "  MQ: " << Median(PE[i]) << "  " << Quantile(PE[i],0.5) << "  " << Quantile(PE[i],-0.5) << std::endl;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [DataFile.dat] [PEFile.dat]s" << std::endl;
    return 1;
  }

  TString const DataFileName = argv[1];
  std::vector<TString> PEFiles;
  for (int i = 2; i < argc; ++i) {
    PEFiles.push_back(argv[i]);
  }

  CompileLimits(DataFileName, PEFiles);

  return 0;
}
