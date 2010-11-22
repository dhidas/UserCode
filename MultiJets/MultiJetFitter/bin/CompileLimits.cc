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




int CompileLimits (TString const OutFileName, TString const DataFileName, std::vector<TString>&  PEFileNames)
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
        std::cout << mass << std::endl;
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
    //std::reverse(PE[imass].begin(), PE[imass].end());
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
  Line.resize(NMasses);
  for (size_t imass = 0; imass < NMasses; ++imass) {
    Line[imass] = new TLine(Data[imass], 0, Data[imass], Hist[imass]->GetMaximum());
    Line[imass]->SetLineColor(2);
    Line[imass]->SetLineWidth(2);
  }

  // Make canvas, plot histogram, and save canvas
  gStyle->SetOptStat(1111111111);  
  for (size_t imass = 0; imass < NMasses; ++imass) {
    TCanvas Can;
    Can.cd();
    Hist[imass]->Draw("hist");
    Line[imass]->Draw("same");
    sprintf(BUFF, "LimitsPlot_%i.eps", (int) Masses[imass]);
    Can.SaveAs(BUFF);
  }


  FILE*  OutFile = fopen("Limits.dat", "w");
  if (OutFile == NULL) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
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


  for (size_t i = 0; i != NMasses; ++i) std::cout << "MQ: " << Median(PE[i]) << "  " << Quantile(PE[i],0.5) << "  " << Quantile(PE[i],-0.5) << std::endl;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " [Outname.eps] [DataFile.dat] [PEFile.dat]s" << std::endl;
    return 1;
  }

  TString const OutFileName = argv[1];
  TString const DataFileName = argv[2];
  std::vector<TString> PEFiles;
  for (int i = 3; i < argc; ++i) {
    PEFiles.push_back(argv[i]);
  }

  CompileLimits(OutFileName, DataFileName, PEFiles);

  return 0;
}
