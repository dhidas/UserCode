////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Fri Jun 10 13:50:29 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraph.h"


float AcceptanceFromFit (TH1F* Hist, int const NGenerated, TString const Name)
{
    //TF1 Func("landgaus", "[0] * TMath::Landau(x, [1], [2], 1) + [3] * TMath::Gaus(x, [4], [5], 1)", 0, 800);
    //Func.SetParameter(0, 100);
    //Func.SetParameter(1, 170);
    //Func.SetParameter(2, 100);
    //Func.SetParameter(3, 40);
    //Func.SetParameter(4, 170);
    //Func.SetParameter(5, 12);
    //Func.SetParLimits(1, 160, 180);
    //Func.SetParLimits(2, 100, 130);
    //Func.SetParLimits(4, 160, 180);
    //Func.SetParLimits(5, 10, 15);

    TF1 Func2("Other", "[0] * TMath::Gaus(x, [1], [2], 1) + [3] * x*x*x + [4] * x*x + [5] * x + [6]", 160, 180);
    //Func2.SetParameter(0, 100);
    //Func2.SetParameter(1, 172);
    //Func2.SetParameter(2, 14);
    Func2.SetParameter(0,  1.68307e+04);
    Func2.SetParameter(1,  1.69423e+02);
    Func2.SetParameter(2,  1.46436e+01);
    Func2.SetParameter(3,  6.31311e-05);
    Func2.SetParameter(4, -5.93896e-02);
    Func2.SetParameter(5,  1.65220e+01);
    Func2.SetParameter(6, -1.06987e+03);
    Func2.SetParLimits(0, 0, 9000000);
    Func2.SetParLimits(1, 150, 200);
    Func2.SetParLimits(2, 10, 20);

    Hist->Fit("Other", "MLPQ", "", 100, 400);


    TF1 FuncBG("BG", "[3] * x*x*x + [4] * x*x + [5] * x + [6]", 100, 400);
    FuncBG.SetParameter(3, Func2.GetParameter(3));
    FuncBG.SetParameter(4, Func2.GetParameter(4));
    FuncBG.SetParameter(5, Func2.GetParameter(5));
    FuncBG.SetParameter(6, Func2.GetParameter(6));
    FuncBG.SetLineColor(2);
    FuncBG.SetLineStyle(2);



    TCanvas Can;
    Can.cd();
    Hist->Draw();
    FuncBG.Draw("same");
    TString const HistName = Hist->GetName();
    TString SaveName = Name + "_" + HistName( HistName.Last('/')+1, HistName.Length() - HistName.Last('/') - 1) + ".gif";
    Can.SaveAs(SaveName, "Q");

    return (Func2.GetParameter(0) / Hist->GetBinWidth(1)) / NGenerated;
}



int RoughAcceptance (TString const InFileName, int const NGenerated, float const CrossSection, float const Luminosity)
{
  // Open input file
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open file: " << InFileName << std::endl;
    throw;
  }

  TString Name = InFileName.Contains("/") ? InFileName( InFileName.Last('/')+1, InFileName.Length() - InFileName.Last('/') - 6) : InFileName;

  FILE* fOut = fopen( (Name+".out").Data(), "w");

  std::set< std::pair<int, float> > AccSet;
  for (int icut = 100; icut != 105; icut += 5) {
    TString const HistName = TString::Format("LeptonPlusJets/TriJetMass_NJetGE04_dP%03i", icut);
    TH1F* Hist = (TH1F*) InFile.Get(HistName);

    float const Acceptance = AcceptanceFromFit(Hist, NGenerated, Name);
    float const NExpected = CrossSection * Luminosity * Acceptance;
    //printf("%12.3E %12.3E %12.3E\n", CrossSection, Luminosity, Acceptance);
    fprintf(fOut, "%15s %3i %9.3E %9.3E\n", Name.Data(), icut, Acceptance, NExpected);

    AccSet.insert( std::make_pair<int, float>(icut, Acceptance) );

  }
  fclose(fOut);

  TGraph Graph( (int) AccSet.size() );
  int Point = 0;
  for (std::set< std::pair<int, float> >::iterator it = AccSet.begin(); it != AccSet.end(); ++it) {
    Graph.SetPoint( Point, (float) it->first, it->second );
    ++Point;
  }

  TCanvas Can;
  Can.cd();
  Graph.SetTitle(Name);
  Graph.GetXaxis()->SetTitle("Diag cut");
  Graph.GetYaxis()->SetTitle("Acceptance");
  Graph.Draw("AC*");
  Can.SaveAs("Acc_"+Name+".gif");


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TString const InFileName = argv[1];

  //RoughAcceptance(InFileName, 1);

  TString Files[9] = {
    "~/TestPlots/ZprimeM500PAT.root",
    "~/TestPlots/ZprimeM750PAT.root",
    "~/TestPlots/ZprimeM1000PAT.root",
    "~/TestPlots/ZprimeM1250PAT.root",
    "~/TestPlots/ZprimeM1500PAT.root",
    "~/TestPlots/ZprimeM2000PAT.root",
    "~/TestPlots/ZprimeM3000PAT.root",
    "~/TestPlots/ZprimeM4000PAT.root",
    "~/TestPlots/TTJets.root"
  };

  int NGen[9] = {
    213112,
    207843,
    108713,
    225142,
    191228,
    234087,
    127665,
    131049,
    615699
  };

  float CrossSection[9] = {
    55,
    15,
     4,
     2,
     1,
     0.1,
     0.01,
     0.001,
   157.5
  };

  for (int i = 0; i != 9; ++i) {
    RoughAcceptance(Files[i], NGen[i], CrossSection[i], 1000);
  }


  

  return 0;
}
