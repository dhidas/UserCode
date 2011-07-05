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
#include "DHidasLJAna/LeptonPlusJets/interface/DHidasJSON.h"




void PlotSingleLeptonEvents (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  if (Ev.Lep.size() != 1) {
    return;
  }

  static TAnaHist Hist(&OutFile, "SingleLeptonEvents");

  TString const LepType = Dtuple::LeptonEventType(Ev);

  Hist.FillTH1D("LeptonPt", 100, 0, 600, Ev.Lep[0].Pt());
  Hist.FillTH1D("LeptonPt_"+LepType, 100, 0, 600, Ev.Lep[0].Pt());

  return;
}


void PlotGammaJetJet (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  static TAnaHist Hist(&OutFile, "PlotGammaJetJet");

  size_t const NJets = Ev.Jet.size();
  size_t const NLeptons = Ev.Lep.size();
  size_t const NPhotons = Ev.Pho.size();

  if (NLeptons != 1 || NPhotons != 1 || NJets < 2) {
    return;
  }

  TString const LepType = Dtuple::LeptonEventType(Ev);

  Hist.FillTH1D("PhotonPt", 100, 0, 800, Ev.Pho[0].Pt());
  Hist.FillTH1D("PhotonPt_"+LepType, 100, 0, 800, Ev.Pho[0].Pt());
  Hist.FillTH1D("LeptonPt", 100, 0, 800, Ev.Lep[0].Pt());
  Hist.FillTH1D("LeptonPt_"+LepType, 100, 0, 800, Ev.Lep[0].Pt());

  if (NJets == 2) {
    Hist.FillTH1D("GammaJJ", 100, 0, 800, (Ev.Pho[0] + Ev.Jet[0] + Ev.Jet[1]).M());
  }

  for (size_t i = 0; i < Ev.Jet.size() - 1; ++i) {
    for (size_t j = i + 1; j < Ev.Jet.size(); ++j) {
      Hist.FillTH1D("GammaJJ_All", 100, 0, 800, (Ev.Pho[i] + Ev.Jet[j] + Ev.Jet[1]).M());
    }
  }


  return;
}


void PlotLeptonPlusJets (Dtuple::SimpleEvent& Ev, TFile& OutFile)
{
  static TAnaHist Hist(&OutFile, "LeptonPlusJets");

  size_t const NJets = Ev.Jet.size();

  TString const LepType = Dtuple::LeptonEventType(Ev);

  /*
  if (NJets == 0 && Ev.Lep.size() == 1 && LepType == "m" && Ev.Lep[0].Pt() > 25) {
    std::cout << "ASD" << std::endl;
  }
  */
  /*
  if (NJets >= 3 && Ev.Lep.size() == 1 && LepType == "m" && Ev.Lep[0].Pt() > 25) {
    printf("%12i  %9.1f  %3i", Ev.Event, Ev.Lep[0].Pt(), (int) Ev.Jet.size());
    for (size_t i = 0; i != Ev.Jet.size(); ++i) {
      printf("  %9.1f", Ev.Jet[i].Pt());
    }
    printf("\n");
  }
  */

  if (NJets < 4 || Ev.Lep.size() != 1) {
    return;
  }

  if (LepType == "e" && Ev.Lep[0].Pt() < 45) {
    return;
  } else if (LepType == "m" && Ev.Lep[0].Pt() < 30) {
    return;
  } else if (LepType != "e" && LepType != "m") {
    std::cerr << "I don't understand this type: PlotLeptonPlusJets(): " << LepType << std::endl;
    return;
  }
 


  static char NAME[400];
  static char SIGN;



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

  Hist.FillTH1D("LeptonPt_"+LepType, 100, 0, 800, Ev.Lep[0].Pt());


  for (size_t i = 0; i < NJets - 2; ++i) {
    for (size_t j = i+1; j < NJets - 1; ++j) {
      for (size_t k = j+1; k < NJets; ++k) {
        float const Mass      = (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).M();
        float const SumPtJets = Ev.Jet[i].Pt() + Ev.Jet[j].Pt() + Ev.Jet[k].Pt();
        float const VecSumPtJets = (Ev.Jet[i] + Ev.Jet[j] + Ev.Jet[k]).Pt();

        //if (Mass > 168 && Mass < 178) {
        //  if (Mass < SumPtJets - 80) {
        //    Hist.FillTH1D("TripletPt_TopMass_cut80", 100, 0, 800, (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).Pt());
        //  }
        //}

        Hist.FillTH2D("TriJetSumPt_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
        Hist.FillTH1D("TriJetMass", 100, 0, 800, Mass);
        Hist.FillTH1D("LeptonPt", 100, 0, 800, Ev.Lep[0].Pt());

        //Hist.FillTH2D("TriJetSumPt_vs_Mass_"+LepType, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
        Hist.FillTH1D("TriJetMass_"+LepType, 100, 0, 800, Mass);

        Hist.FillTH1D("TripletPt", 100, 0, 800, (Ev.Jet[i]+Ev.Jet[j]+Ev.Jet[k]).Pt());
        Hist.FillTH2D("TriJetMas_vs_TripletPt", 1000, 0, 1000, 1000, 0, 1000, VecSumPtJets, Mass);

        //Hist.FillTH2D("METMag_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, Ev.MET.Mod(), Mass);

        Hist.FillTH2D("TriVecJetSumPt_vs_Mass", 1000, 0, 1000, 1000, 0, 1000, VecSumPtJets, Mass);

        for (float JetPtCut = 20; JetPtCut < 45; JetPtCut += 10) {
          // Cut on 4-th jet Pt
          if (Ev.Jet[3].Pt() < JetPtCut) {
            continue;
          }
          if (Ev.Jet[k].Pt() >= JetPtCut) {
            sprintf(NAME, "TriJetSumPt_vs_Mass_JetPtCut%03i", (int) JetPtCut);
            //Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, SumPtJets, Mass);
            sprintf(NAME, "TriJetMass_JetPtCut%03i", (int) JetPtCut);
            Hist.FillTH1D(NAME, 100, 0, 800, Mass);
          }
        }


        for (size_t ijets = 4; ijets < 10; ++ijets) {
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

            for (float JetPtCut = 20; JetPtCut < 45; JetPtCut += 10) {
              if (Ev.Jet[3].Pt() < JetPtCut) {
                continue;
              }
              if (Ev.Jet[k].Pt() >= JetPtCut) {
                sprintf(NAME, "TriJetMass_d%c%03is_JetPtCut%03i", SIGN, abs((int) cut), (int) JetPtCut);
              Hist.FillTH1D(NAME, 100, 0, 800, Mass);
              }
            }

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
            //Hist.FillTH2D(NAME, 1000, 0, 1000, 1000, 0, 1000, Ev.MET.Mod(), Mass);


            for (size_t ijets = 4; ijets < NJets; ++ijets) {
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

  for (size_t i = 0; i != Ev.Lep.size(); ++i) {
    if ((*Ev.LeptonType)[i] == Dtuple::kElectron && Ev.Lep[i].Pt() < 45) {
      return;
    } else if ((*Ev.LeptonType)[i] == Dtuple::kMuon && Ev.Lep[i].Pt() < 30) {
      return;
    }
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




int RunLJDtuple (bool const IsData, TString const OutFileName, std::vector<TString> const& InFileNames)
{
  if (IsData) {
    std::cout << "I see that this is DATA" << std::endl;
  } else {
    std::cout << "This looks like MC..." << std::endl;
  }

  // Open OutFile
  TFile OutFile(OutFileName, "create");
  if (!OutFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output file: " << OutFileName << std::endl;
    throw;
  }

  DHidasJSON JSON("json/June27thGoldJSON.txt", IsData);

  // Keep track of runs
  DHRunTracker RT;


  // Grab small ntuple.  Everything is contained in "Ev"
  Dtuple DT(InFileNames);
  Dtuple::SimpleEvent& Ev = DT.GetEvt();
  for (long long ientry = 0; DT.GetEntry(ientry) > 0; ++ientry) {
    if (ientry % 10000 == 0) {
      std::cout << "Processing: " << ientry << std::endl;
    }

    // Check lumi section
    if (!JSON.IsGoodLumiSection(Ev.Run, Ev.LumiSection)) {
      printf("JSON Skipping Run %9i  LumiSection %9i\n", Ev.Run, Ev.LumiSection);
      continue;
    }

    // Duplicate event check
    if (RT.IsDuplicate(Ev.Run, Ev.Event)) {
      std::cerr << "Duplicate run found and being skipped: " << Ev.Run << "  " << Ev.Event << std::endl;
      continue;
    }

    // Analysis starts here!
    PlotDileptonEvents(Ev, OutFile);
    PlotLeptonPlusJets(Ev, OutFile);
    PlotGammaJetJet(Ev, OutFile);
    PlotSingleLeptonEvents(Ev, OutFile);





    DT.ClearDtuple();
  }

  OutFile.Write();
  OutFile.Close();
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [IsData] [OutFileName] [InFileName]s" << std::endl;
    return 1;
  }

  bool const IsData = atoi(argv[1]) == 1 ? true : false;
  TString const OutFileName = argv[2];
  std::vector<TString> InFileNames;
  for (int i = 3; i < argc; ++i) {
    InFileNames.push_back(argv[i]);
  }

  RunLJDtuple(IsData, OutFileName, InFileNames);

  return 0;
}
