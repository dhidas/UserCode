////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Oct 13 03:38:37 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>

#include "StandardHypoTestInvDemo.h"

#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"


int AddHTIR (int const SignalMass, std::vector<TString>& InFileNames)
{
  // Cumulative result
  RooStats::HypoTestInverterResult* MyResult = 0x0;

  // Grab the result from each file and Add it to MyResult
  for (size_t ifile = 0; ifile != InFileNames.size(); ++ifile) {
    TFile InFile(InFileNames[ifile], "read");
    if (!InFile.IsOpen()) {
      std::cerr << "ERROR: cannot open input file: " << InFileNames[ifile] << std::endl;
      exit(1);
    }

    RooStats::HypoTestInverterResult* ThisResult = (RooStats::HypoTestInverterResult*) InFile.Get("result_xs");

    // Check if this object exists.. if not, it might not be done running, and so on
    if (!ThisResult) {
      continue;
    }

    if (MyResult == 0x0) {
      MyResult = ThisResult;
    } else {
      MyResult->Add(*ThisResult);
    }
  }

  // Open output file for limits
  std::ofstream OutFile(TString::Format("Limits_CLs_%i.dat", (int) SignalMass).Data());
  if (!OutFile.is_open()) {
    std::cerr << "ERROR: cannot open output file Limits_CLs_*.dat for mass: " << SignalMass << std::endl;
    throw;
  }

  // Setup output root file based on mass
  TFile OutRootFile(TString::Format("MultiJetsCLs_%i.root", (int) SignalMass), "recreate");
  if (!OutRootFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output root file for mass: " << SignalMass << std::endl;
    throw;
  }
  OutRootFile.cd();

  // Number of entries in result
  const int NEntries = MyResult->ArraySize();

  // Just some names
  const char* TypeName = "Hybrid";
  const char* ResultName = MyResult->GetName();
  TString PlotTitle = TString::Format("%s CL Scan for workspace %s", TypeName, ResultName);

  // Grab the result plot
  RooStats::HypoTestInverterPlot *Plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", PlotTitle, MyResult);
  TCanvas CanCLb("CLb_2CL", "CLb_2CL");
  CanCLb.cd();
  Plot->SetTitle(TString::Format("Hybrid CL Scan for M_{jjj} = %i", (int) SignalMass));
  Plot->Draw("CLb 2CL");  // plot all and Clb
  CanCLb.SaveAs(TString::Format("CLb2L_%i.eps", (int) SignalMass));
  OutRootFile.cd();
  CanCLb.Write();

  // Draw the sampling distributions
  TCanvas CanHTI("HTI_Result", "HTI_Result");
  CanHTI.Divide(2, (int) TMath::Ceil(NEntries/2));
  for (int i = 0; i < NEntries; ++i) {
    CanHTI.cd(i + 1);
    RooStats::SamplingDistPlot * SamplingPlot = Plot->MakeTestStatPlot(i);
    SamplingPlot->SetLogYaxis(true);
    delete SamplingPlot;
  }
  CanHTI.SaveAs(TString::Format("HTI_Result_%i.eps", (int) SignalMass));
  OutRootFile.cd();
  CanHTI.Write();

  // Print the limits
  printf(" expected limit (-2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-2));
  printf(" expected limit (-1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(-1));
  printf(" expected limit (median) %12.3E\n", MyResult->GetExpectedUpperLimit(0) );
  printf(" expected limit (+1 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(1) );
  printf(" expected limit (+2 sig) %12.3E\n", MyResult->GetExpectedUpperLimit(2) );
  printf(" observed limit          %12.3E +/- %12.3E\n", MyResult->UpperLimit(), MyResult->UpperLimitEstimatedError()); 


  // Write results and close file
  OutFile << SignalMass << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-2) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit(-1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 0) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 1) << std::endl;
  OutFile << MyResult->GetExpectedUpperLimit( 2) << std::endl;
  OutFile << MyResult->UpperLimit() << std::endl;
  OutFile.close();

  // Close root file
  OutRootFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [FileName]s" << std::endl;
    return 1;
  }

  // Map for filenames
  std::map<int, std::vector<TString> > MassFileMap;

  // put each file with it's mass buddies
  for (int i = 1; i != argc; ++i) {
    TString Name = argv[i];
    int const Mass = atoi( Name(Name.Last('_') + 1, Name.Last('.') - Name.Last('_') - 1).Data() );
    MassFileMap[Mass].push_back(Name);
  }

  // For each mass calculate a limit
  for (std::map<int, std::vector<TString> >::iterator it = MassFileMap.begin(); it != MassFileMap.end(); ++it) {
    AddHTIR(it->first, it->second);
  }


  return 0;
}
