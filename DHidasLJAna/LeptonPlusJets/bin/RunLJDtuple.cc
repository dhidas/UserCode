////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon May  2 04:36:37 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "TString.h"
#include "TFile.h"

#include "DHidasLJAna/LeptonPlusJets/interface/Dtuple.h"
#include "DHidasLJAna/LeptonPlusJets/interface/TAnaHist.h"
#include "DHidasLJAna/LeptonPlusJets/interface/DHRunTracker.h"




void PlotLeptonPlusJets (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  static TAnaHist Hist(&OutFile, "LeptonPlusJets");

  size_t const NJets = Ev.Jet.size();

  if (NJets < 3 || Ev.Lep.size() != 1) {
    return;
  }


  static char NAME[400];
  static char SIGN;

  TString const LepType = Dtuple::LeptonEventType(Ev);


  // Make a plot of the Jet Pts
  for (size_t i = 0; i < NJets; ++i) {
    Hist.FillTH1D("JetPt", 100, 0, 800, Ev.Jet[i].Pt());
    Hist.FillTH1D("JetPt_"+LepType, 100, 0, 800, Ev.Jet[i].Pt());
    sprintf(NAME, "JetPt%02i", (int) i);
    Hist.FillTH1D(NAME, 100, 0, 800, Ev.Jet[i].Pt());
    sprintf(NAME, "JetPt%02i_", (int) i);
    Hist.FillTH1D(NAME+LepType, 100, 0, 800, Ev.Jet[i].Pt());
  }

  Hist.FillTH1D("METMag", 100, 0, 800, Ev.MET.Mod());
  sprintf(NAME, "METMag_NJetEQ%02i_", (int) NJets);
  Hist.FillTH1D(NAME, 100, 0, 800, Ev.MET.Mod());
  Hist.FillTH1D("METMag_"+LepType, 100, 0, 800, Ev.MET.Mod());

  for (size_t i = 0; i < NJets - 2; ++i) {
    for (size_t j = i+1; j < NJets - 1; ++j) {
      for (size_t k = j+1; k < NJets; ++k) {
        float const Mass      = (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).M();
        float const SumPtJets = Ev.Jet[i].Pt() + Ev.Jet[j].Pt() + Ev.Jet[k].Pt();
        Hist.FillTH2D("TriJetSumPt_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
        Hist.FillTH1D("TriJetMass", 100, 0, 800, Mass);
        Hist.FillTH1D("LeptonPt", 100, 0, 800, Ev.Lep[0].Pt());

        Hist.FillTH2D("TriJetSumPt_vs_Mass_"+LepType, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
        Hist.FillTH1D("TriJetMass_"+LepType, 100, 0, 800, Mass);
        Hist.FillTH1D("LeptonPt_"+LepType, 100, 0, 800, Ev.Lep[0].Pt());

        Hist.FillTH1D("TripletPt", 100, 0, 800, (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).Pt());
        Hist.FillTH2D("TriJetMas_vs_TripletPt", 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);

        Hist.FillTH2D("METMag_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, Ev.MET.Mod(), Mass);

        for (size_t ijets = 3; ijets < 10; ++ijets) {
          if (NJets == ijets) {
            sprintf(NAME, "TriJetSumPt_vs_Mass_NJetEQ%02i", (int) ijets);
            Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
            sprintf(NAME, "TriJetMass_NJetEQ%02i", (int) ijets);
            Hist.FillTH1D(NAME, 100, 0, 800, Mass);
          } 
        }

        for (size_t ijets = 3; ijets < 10; ++ijets) {
          if (NJets >= ijets) {
            sprintf(NAME, "TriJetSumPt_vs_Mass_NJetGE%02i", (int) ijets);
            Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
            sprintf(NAME, "TriJetMass_NJetGE%02i", (int) ijets);
            Hist.FillTH1D(NAME, 100, 0, 800, Mass);
          } 
        }


        // Ones passing some cut!
        for (float cut = -100; cut <= 200; cut += 5) {
          if (cut >= 0) {
            SIGN = 'P';
          } else {
            SIGN = 'N';
          }

          if (Mass < SumPtJets - cut) {
            //sprintf(NAME, "TriJetSumPt_vs_Mass_d%c%03i", SIGN, abs((int) cut));
            //Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
            sprintf(NAME, "TriJetMass_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Mass);

            //sprintf(NAME, "TriJetSumPt_vs_Mass_NJet%02i_d%c%03i", (int) NJets, SIGN, abs((int) cut));
            //Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
            sprintf(NAME, "TriJetMass_NJetEQ%02i_d%c%03i", (int) NJets, SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Mass);

            // Lepton Pt
            sprintf(NAME, "LeptonPt_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Ev.Lep[0].Pt());

            // Jet Pt plots for the 3 jets
            sprintf(NAME, "JetPt_TripletJet0_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Ev.Jet[i].Pt());
            sprintf(NAME, "JetPt_TripletJet1_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Ev.Jet[j].Pt());
            sprintf(NAME, "JetPt_TripletJet2_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Ev.Jet[k].Pt());

            // SystemPt
            sprintf(NAME, "TripletPt_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).Pt());

            // MET
            sprintf(NAME, "METMag_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH1D(NAME, 100, 0, 800, Ev.MET.Mod());

            // Met vs Mass
            sprintf(NAME, "METMag_vs_Mass_d%c%03i", SIGN, abs((int) cut));
            Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, Ev.MET.Mod(), Mass);


            for (size_t ijets = 3; ijets < NJets; ++ijets) {
              if (NJets >= ijets) {
                //sprintf(NAME, "TriJetSumPt_vs_Mass_NJetGE%02i_d%c%03i", (int) ijets, SIGN, abs((int) cut));
                //Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
                sprintf(NAME, "TriJetMass_NJetGE%02i_d%c%03i", (int) ijets, SIGN, abs((int) cut));
                Hist.FillTH1D(NAME, 100, 0, 800, Mass);
              } 
            }
          }
        }


      }
    }
  }


  return;
}





void PlotDileptonEvents (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  static TAnaHist Hist(&OutFile, "DileptonEvents");

  if (Ev.NLeptons != 2) {
    return;
  }

  TString const DilType = Dtuple::LeptonEventType(Ev);

  Hist.FillTH1F("DileptonMass", "Dilepton Mass", "Dilepton Mass (GeV)", "Events", 250, 0, 500, (Ev.Lep[0]+Ev.Lep[1]).M());
  Hist.FillTH1F("LeptonPt0", "Leading Lepton P_{T}", "Lepton P_{T} (GeV)", "Events", 250, 0, 500, Ev.Lep[0].Pt());
  Hist.FillTH1F("LeptonPt1", "Sub-Leading Lepton P_{T}", "Lepton P_{T} (GeV)", "Events", 250, 0, 500, Ev.Lep[1].Pt());

  Hist.FillTH1F("DileptonMass_"+DilType, "Dilepton Mass "+DilType, "Dilepton Mass (GeV)", "Events", 250, 0, 500, (Ev.Lep[0]+Ev.Lep[1]).M());
  Hist.FillTH1F("LeptonPt0_"+DilType, "Leading Lepton P_{T} "+DilType, "Lepton P_{T} (GeV)", "Events", 250, 0, 500, Ev.Lep[0].Pt());
  Hist.FillTH1F("LeptonPt1_"+DilType, "Sub-Leading Lepton P_{T} "+DilType, "Lepton P_{T} (GeV)", "Events", 250, 0, 500, Ev.Lep[1].Pt());
  return;
}




int RunLJDtuple (TString const OutFileName, std::vector<TString> const& InFileNames)
{
  // Open OutFile
  TFile OutFile(OutFileName, "create");
  if (!OutFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
    throw;
  }

  // Keep track of runs
  DHRunTracker RT;


  // Grab small ntuple.  Everything is contained in "Ev"
  Dtuple DT(InFileNames);
  Dtuple::SimpleEvent& Ev = DT.GetEvt();
  for (long long ientry = 0; DT.GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 10000 == 0) {
      std::cout << "Processing: " << ientry << std::endl;
    }

    // Duplicate event check
    if (RT.IsDuplicate(Ev.Run, Ev.Event)) {
      std::cerr << "Duplicate run found and being skipped: " << Ev.Run << "  " << Ev.Event << std::endl;
      continue;
    }

    // Analysis starts here!
    PlotDileptonEvents(Ev, OutFile);
    PlotLeptonPlusJets(Ev, OutFile);




    DT.ClearDtuple();
  }

  OutFile.Write();
  OutFile.Close();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName] [InFileName]s" << std::endl;
    return 1;
  }

  TString const OutFileName = argv[1];
  std::vector<TString> InFileNames;
  for (int i = 2; i < argc; ++i) {
    InFileNames.push_back(argv[i]);
  }

  RunLJDtuple(OutFileName, InFileNames);

  return 0;
}
