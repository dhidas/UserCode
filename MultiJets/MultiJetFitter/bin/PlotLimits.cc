#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>



#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveLabel.h"
#include "TLine.h"
#include "TROOT.h"


void DrawLimits (std::vector<TString> const&);

int main (int argc, char* argv[])
{
  if (argc <= 1) {
    std::cerr << "Usage: " << argv[0] << " [FileName]s" << std::endl;
    return 1;
  }


  // Add the filenames given to the vector
  std::vector<TString> FileNames;
  for (int i=1; i != argc; ++i) {
    std::cout << "Adding file: " << argv[i] << std::endl;
    FileNames.push_back(argv[i]);
  }

  DrawLimits(FileNames);

  return 0;
}



void DrawLimits (std::vector<TString> const& FileNames)
{
  // Check the size of the input vector
  if (FileNames.size() == 0) {
    std::cerr << "ERROR: No files given" << std::endl;
    return;
  }

  // Should we draw the observed line?
  bool const DrawObserved = true;

  // Set the default style to plain
  gROOT->SetStyle("Plain");

  // Make the canvas!
  TCanvas HiggsLimits("Limits","Limits",200,10,700,500);
  HiggsLimits.SetGrid();
  if (false) {
    HiggsLimits.DrawFrame(150,0.1,450,190);
    HiggsLimits.SetLogy(0);
  } else {
    HiggsLimits.DrawFrame(150,0,450,190);
    HiggsLimits.SetLogy(0);
  }
  HiggsLimits.cd();

  // Setup the Legend
  TLegend MyLegend(0.66,0.6,0.88,0.88, "");
  MyLegend.SetBorderSize(0);
  MyLegend.SetFillColor(0);


  // TGraphs for shading.  Made later.  Want to add them near last to legend...
  TGraph *UseShade1S;
  TGraph *UseShade2S;
  std::string ShadeLabel;


  // Loop over all filenames
  for (size_t ifile=0; ifile != FileNames.size(); ++ifile) {

    // Open the file and check it opened
    std::ifstream InFile(FileNames[ifile].Data());
    if (!InFile.is_open()) {
      std::cerr << "ERROR: File cannot be opened: " << FileNames[ifile] << std::endl;
      return;
    }
    std::cout << "Reading file: " << FileNames[ifile] << std::endl;


    // Let's attempt to read this file line by line

    // Get the Label (first line)
    std::string MyLabel;
    std::getline(InFile, MyLabel);
    std::cout << "Label for this file: " << MyLabel << std::endl;

    // For the rest of the lines
    std::string InLine;
    std::istringstream InStreamLine;

    // Read the masses lines and count how many there are
    std::getline(InFile, InLine);
    std::cout << "Masses for this file:\n   " << InLine << std::endl;
    InStreamLine.clear();
    InStreamLine.str(InLine);
    int NMass = 0;
    for (float mass; InStreamLine >> mass; ++NMass) {
    }
    std::cout << MyLabel << " NMass = " << NMass << std::endl;

    // Now set the Mass vector (x-axis)
    float Mass[NMass];
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> Mass[imass];
      std::cout << MyLabel << " Mass: " << Mass[imass] << std::endl;
      
    }


    // Get the -2 sigma line
    float M2Sigma[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> M2Sigma[imass];
      std::cout << MyLabel << " M2Sigma[" << imass << "]: " << M2Sigma[imass] << std::endl;
    }

    // Get the -1 sigma line
    float M1Sigma[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> M1Sigma[imass];
      std::cout << MyLabel << " M1Sigma[" << imass << "]: " << M1Sigma[imass] << std::endl;
    }

    // Get the Median line
    float Median[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> Median[imass];
      std::cout << MyLabel << " Median[" << imass << "]: " << Median[imass] << std::endl;
    }

    // Get the +1 sigma line
    float P1Sigma[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> P1Sigma[imass];
      std::cout << MyLabel << " P1Sigma[" << imass << "]: " << P1Sigma[imass] << std::endl;
    }

    // Get the +2 sigma line
    float P2Sigma[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> P2Sigma[imass];
      std::cout << MyLabel << " P2Sigma[" << imass << "]: " << P2Sigma[imass] << std::endl;
    }

    // Get the Observed line
    float Observed[NMass];
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    for (int imass=0; imass != NMass; ++imass) {
      InStreamLine >> Observed[imass];
      if (DrawObserved) {
        std::cout << MyLabel << " Observed[" << imass << "]: " << Observed[imass] << std::endl;
      }
    }

    // Get the Color to use
    int Color;
    std::getline(InFile, InLine);
    InStreamLine.clear();
    InStreamLine.str(InLine);
    InStreamLine >> Color;
    std::cout << MyLabel << " will be Color = " << Color << std::endl;


    // Make the TGraphs we want
    TGraph *grMedian = new TGraph(NMass, Mass, Median);
    TGraph *grObserved = new TGraph(NMass, Mass, Observed);


    // Make the shade graphs.
    TGraph *grShade1S = new TGraph(2*NMass);
    TGraph *grShade2S = new TGraph(2*NMass);
    for (int i=0; i < NMass; ++i) {
      grShade1S->SetPoint(i, Mass[i], P1Sigma[i]);
      grShade1S->SetPoint(NMass+i, Mass[NMass-i-1], M1Sigma[NMass-i-1]);
      grShade2S->SetPoint(i, Mass[i], P2Sigma[i]);
      grShade2S->SetPoint(NMass+i, Mass[NMass-i-1], M2Sigma[NMass-i-1]);
    }
    grShade1S->SetFillColor(3);
    grShade2S->SetFillColor(5);

    // Draw the shades only for first file
    if (ifile == 0) {
      std::cout << "Will draw +/- 1.2 sigma for: " << MyLabel << std::endl;
      UseShade2S = grShade2S;
      UseShade1S = grShade1S;
      UseShade2S->Draw("f");
      UseShade1S->Draw("f");
      ShadeLabel = MyLabel;
    }

    // Draw the median
    grMedian->SetLineWidth(3);
    grMedian->SetFillColor(0);
    grMedian->SetLineStyle(4);
    grMedian->SetMarkerColor(Color);
    grMedian->SetLineColor(Color);
    grMedian->SetMarkerStyle(0);
    grMedian->Draw("l");
    MyLegend.AddEntry(grMedian, (MyLabel+" Expectation").c_str(), "l");

    // Draw the Observed
    if (DrawObserved) {
      grObserved->SetLineWidth(3);
      grObserved->SetFillColor(0);
      grObserved->SetLineStyle(1);
      grObserved->SetMarkerColor(Color);
      grObserved->SetLineColor(Color);
      grObserved->SetMarkerStyle(0);
      grObserved->Draw("l");
      MyLegend.AddEntry(grObserved, (MyLabel+" Observed").c_str(), "l");
    }

  }

  // Draw the sigma shaded bands and add to legend
  MyLegend.AddEntry(UseShade1S, (ShadeLabel+" \\pm 1\\sigma").c_str(), "f");
  MyLegend.AddEntry(UseShade2S, (ShadeLabel+" \\pm 2\\sigma").c_str(), "f");



  // Add the SM line
  TLine *SMLine = new TLine(110, 1, 200, 1);
  SMLine->SetLineWidth(3);
  SMLine->SetLineStyle(1);
  SMLine->SetLineColor(4);
  //SMLine->Draw();
  //MyLegend.AddEntry(SMLine, "SM", "l");

  TPaveLabel *HMassLabel = new TPaveLabel();
  HMassLabel->SetLabel("Mass [GeV/c^{2}]");
  HMassLabel->SetX1NDC(0.60);
  HMassLabel->SetX2NDC(0.90);
  HMassLabel->SetY1NDC(0.01);
  HMassLabel->SetY2NDC(0.05);
  HMassLabel->SetFillColor(0);
  HMassLabel->SetTextSize();
  HMassLabel->Draw("same");

  TPaveLabel *YLabel = new TPaveLabel();
  YLabel->SetLabel("Expected 95\% C.L.# of Events");
  YLabel->SetX1NDC(0.00);
  YLabel->SetX2NDC(0.05);
  YLabel->SetY1NDC(0.10);
  YLabel->SetY2NDC(0.90);
  YLabel->SetFillColor(0);
  YLabel->SetTextAngle(90);
  YLabel->SetTextAlign(22);
  YLabel->SetTextSize(0.07);
  YLabel->Draw("same");

  TPaveLabel *Title = new TPaveLabel();
  Title->SetLabel("Higgs cross section limits");
  Title->SetX1NDC(0.25);
  Title->SetX2NDC(0.75);
  Title->SetY1NDC(0.92);
  Title->SetY2NDC(0.99);
  Title->SetFillColor(0);
  Title->SetTextAlign(22);
  Title->SetTextSize();
  //Title->Draw("same");

  TPaveLabel *RunII = new TPaveLabel();
  //RunII->SetLabel("CDF Run II Preliminary");
  RunII->SetLabel("CDF Run II UnBlessed");
  RunII->SetX1NDC(0.1);
  RunII->SetX2NDC(0.5);
  RunII->SetY1NDC(0.92);
  RunII->SetY2NDC(0.99);
  RunII->SetFillColor(0);
  //RunII->SetTextAlign(22);
  RunII->SetTextSize(0.5);
  //RunII->Draw("same");

  TPaveLabel *Lumi = new TPaveLabel();
  Lumi->SetLabel("#int L = 1.9 fb^{-1}");
  Lumi->SetX1NDC(0.5);
  Lumi->SetX2NDC(0.9);
  Lumi->SetY1NDC(0.92);
  Lumi->SetY2NDC(0.99);
  Lumi->SetFillColor(0);
  Lumi->SetTextAlign(31);
  Lumi->SetTextSize(0.44);
  //Lumi->Draw("same");



  // Draw the legend and canvas.  Save to file
  MyLegend.Draw();
  HiggsLimits.Draw();
  HiggsLimits.SaveAs("GraphLimits.eps");

  return;
}



