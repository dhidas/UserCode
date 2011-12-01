/*
 * Run.cc
 *
 *  Created on: Mar 11, 2010
 *      Author: lkreczko
 */
//#include "cms_bristol_ana_v3.hh"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "Analysis.h"
#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <fstream>

#include "Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.C"
#include "Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.C"
#include "Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.C"

using namespace ROOT;
using namespace std;
using namespace BAT;

void setUpOnce() {
    //needed to proper link vector<float> etc.
    gROOT->ProcessLine("#include <vector>");
    //prevent automatic ownership of ROOT objects
    TH1F::AddDirectory(false);
    //ignore ROOT errors (temporaly due to different nTuple content)
    gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
}

int main(int argc, char **argv) {
    setUpOnce();
    TStopwatch watch;
    watch.Start();

    int const Section = argc > 1 ? atoi(argv[1]) : -1;
    TString const ListName = argc > 2 ? argv[2]  : "";


    if (argc >= 4) {
    	if (strcmp(argv[3], "up") == 0) {
    		Jet::correctDirection = JetCorrDirection::PLUS;
    		cout << "Adding JECUnc to jet p4" << endl;
    	} else if (strcmp(argv[3], "down") == 0) {
    		Jet::correctDirection = JetCorrDirection::MINUS;
    		cout << "Subtracting JECUnc to jet p4" << endl;
    	}
    }
    Analysis::useJetAlgorithm(JetAlgorithm::PF2PAT);
    Analysis::useElectronAlgorithm(ElectronAlgorithm::ParticleFlow);
    // Analysis::useElectronAlgorithm(ElectronAlgorithm::Calo);
		// These are selectecPatElectrons, as Claudia uses

    Analysis::useMuonAlgorithm(MuonAlgorithm::ParticleFlow);
    // Analysis::useMuonAlgorithm(MuonAlgorithm::Default);

    Analysis::useMETAlgorithm(METAlgorithm::ParticleFlowMET);
    // Analysis::luminosity = 1091.0; // luminosity()/1000000.;
    // Analysis::luminosity = 5000.0; // luminosity()/1000000.;
    //Analysis::luminosity = luminosityReReco() / 1000000.;
    Analysis::luminosity = luminosity() / 1000000.;
    cout<<"Luminosity: "<<Analysis::luminosity<<" pb"<<endl;
    // Analysis::luminosity = 13.63;
    // Analysis::luminosity = 36.145;

    Analysis::useCustomConversionTagger(false);
    Analysis::usePFIsolation(true);
    // Analysis::usePFIsolation(false);	// For Calo electrons

    boost::scoped_ptr<Analysis> myAnalysis(new Analysis());
    myAnalysis->setUsedNeutrinoSelectionForTopPairReconstruction(NeutrinoSelectionCriterion::pzClosestToLepton);
    myAnalysis->setUsedTTbarReconstructionCriterion(TTbarReconstructionCriterion::chi2);

    if (Section != -1) {
      myAnalysis->histMan.Section = Section;
    }

    //myAnalysis->addInputFile("/cms/se/store/user/clseitz/OSUNtuples/ttbar/TTJets_TuneZ2_7TeV-madgraph-tauola_srappocc-ttbsm_v9_Summer11-PU_S4_START42_V11-v1-bf57a985b107a689982b667a3f2f23c7_USER_CMSSW_4_2_4_Ov10.6_LQnTuple_34_1_Pzg.root");
    std::vector<TString> FileList;
    if (ListName == "") {
      myAnalysis->addInputFile("/cms/se/store/user/clseitz/OSUNtuples/ttbar/TTJets_TuneZ2_7TeV*.root");
    } else {
      std::ifstream InFile(ListName.Data());
      if (!InFile.is_open()) {
        std::cerr << "ERROR cannot open input list" << std::endl;
        return 2;
      }

      int NFiles = 0;
      for (TString Line; Line.ReadLine(InFile); ++NFiles) {
        FileList.push_back(Line);
      }

      size_t const NSections = 20;
      size_t const NPerSection = FileList.size() / NSections + 1;
      size_t const Start = Section * NPerSection;
      size_t const End   = (Section + 1) * NPerSection;
      for (size_t ifile = Start; ifile < End; ++ifile) {
        if (ifile >= FileList.size()) {
          std::cout << "Reached end of input file list" << std::endl;
          break;
        }
        std::cout << "Adding file: " << FileList[ifile] << std::endl;
        myAnalysis->addInputFile(FileList[ifile]);
      }
    }

    cout << "starting analysis" << endl;
    myAnalysis->analyze(Section, ListName);
    watch.Stop();
    watch.Print();

    return 0;
}
