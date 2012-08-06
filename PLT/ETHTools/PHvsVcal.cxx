//////////////////////////////////////////////////
// Last modified by Jennifer Sibille 12.05.2010 //
// Directories are hard-coded                   //
//////////////////////////////////////////////////

#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "iostream"
#include <fstream>
#include <string.h>
#include <vector>
#include <sstream>
#include <stdio.h>

using namespace std;

int main (int argc, char* argv[50]) 
{

  char dir[200];

  if (argc < 2)
    {
      //in case non sufficient command line input was given
      cout << "Usage: ./PHvsVcal dirName \n";
      return 0;
    }

  strcpy(dir, argv[1]);

  double PH[10], Vcal[10], fitParam[4];
  char pix[5], pixF[5];
  int row, rowF, col, colF, pcol=0;


  //Set Vcal values (fixed)

  //Low range
  Vcal[0] = 50.0;
  Vcal[1] = 100.0;
  Vcal[2] = 150.0;
  Vcal[3] = 200.0;
  Vcal[4] = 250.0;
  //High range
  Vcal[5] = 7.0 * 30.0;
  Vcal[6] = 7.0 * 50.0;
  Vcal[7] = 7.0 * 70.0;
  Vcal[8] = 7.0 * 90.0;
  Vcal[9] = 7.0 * 200.0;
 

  //Read in file
  char infile[200], fitInfile[200];
  string input;
  sprintf(infile, "/home/l_tester/slc5/Output/Efficiency_measurements/%s/phCalibration_C0.dat", dir);
  sprintf(fitInfile, "/home/l_tester/slc5/Output/Efficiency_measurements/%s/phCalibrationFitTan_C0.dat", dir);
  ifstream read_in, fitReadIn;
  read_in.open(infile);
  fitReadIn.open(fitInfile);
  cout << "Reading files: \n" << infile << " \n " << fitInfile << endl;

  char rootFileName[200];
  char grName[100];
  sprintf(rootFileName, "/home/l_tester/slc5/Output/Efficiency_measurements/%s/PHvsVcal.root", dir);
  TFile rootFile(rootFileName, "RECREATE");

  for (int j = 0; j < 3; j++)
    {
      read_in.ignore(100,'\n');
    }

  for (int k = 0; k < 3; k++)
    {
      fitReadIn.ignore(100,'\n');
    }

   while(!read_in.eof())
  {
      for (int i = 0; i < 10; i++)
	{
	  read_in >> input;
	  if(input.compare("N/A")==0)
	    {
	      PH[i] = -500;
	    }
	  else
	    {
	      stringstream(input) >> PH[i];
	    }
	}
      read_in >> pix >> col >> row;
      fitReadIn >> fitParam[0] >> fitParam[1] >> fitParam[2] >> fitParam[3] >> pixF >> colF >> rowF;
      if(col != pcol)
	{
	  cout << "Col " << pcol << endl;
	}
      //cout << col << " " << row << " " << PH[0] << endl;


      //Check that pixels match between fit file and ph file
      if (col==colF && row==rowF)
	{
	  //Make graph
	  TGraph* gr = new TGraph(10, Vcal, PH);
	  gr->SetMarkerStyle(20);
	  gr->SetMarkerColor(4);
	  sprintf(grName, "C%iR%i", col, row);
	  TCanvas* c1 = new TCanvas(grName);
	  c1->cd();
	  gr->Draw("PA");
	  gr->GetXaxis()->SetTitle("VCal");
	  gr->GetYaxis()->SetTitle("PH");
	  gr->SetTitle("Pulse Height vs. Vcal");
	  gr->Draw("PA");

	  //Draw fit on graph
	  TF1* fit = new TF1("fit","[3]+[2]*tanh([0]*x-[1])",0,1600);
	  fit->SetParameters( fitParam[0], fitParam[1], fitParam[2], fitParam[3]);
	  fit->Draw("lsame");
	  c1->Write();
	  
	  //gr->Write();
	  //gr->Clear();
	  pcol = col;
	}
      else 
	{
	  cout << "ERROR: pixel address does not match in PH and Fit files" << endl;
	  break;
	}

  }
  cout << "Done" << endl;
  read_in.close();
  rootFile.Close();
}
