// **************************************************************
// ***** IntroAnalysis.cpp
// ***** Introductory analysis using OAT
// ***** Contains most of the code for a complete analysis
// ***** Important sections are marked with comment blocks.
// ***** Other comments explain important steps.
// ***** Related files:  IntroAnalysis.h for declarations,
// ***** IntroAnalysisUtils.cpp for utility functions.
// **************************************************************


#include <cmath>
#include <math.h>
#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>

#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"

#include "IntroAnalysis.h"
#include "../interface/EventCounter.h"
#include "../interface/Printers/EventTablePrinter.h"

#include "Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON.C"

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

		// **************************************************************
		// *********** Set up **********
		// **************************************************************

    setUpOnce();
    TStopwatch watch;
    watch.Start();

    // *** Set configuration values

    // Set luminosity in the data ntuples
    IntroAnalysis::luminosity = luminosity();

    // Set jet algorithm
    IntroAnalysis::useJetAlgorithm(JetAlgorithm::PF2PAT);

    // Set electron algorithm
    IntroAnalysis::useElectronAlgorithm(ElectronAlgorithm::ParticleFlow);

    // Set muon algorithm
    IntroAnalysis::useMuonAlgorithm(MuonAlgorithm::ParticleFlow);

    // Set MET algorithm
    IntroAnalysis::useMETAlgorithm(METAlgorithm::ParticleFlowMET);

    // Set boolean options

    // Set whether to use the custom conversion tagger
    IntroAnalysis::useCustomConversionTagger(false);

    // Set whether to use particle-flow isolation
    IntroAnalysis::usePFIsolation(true);

    // Instantiate analysis object
    boost::scoped_ptr<IntroAnalysis> myAnalysis(new IntroAnalysis());
    //  myAnalysis->setMaximalNumberOfEvents(100);

    // Set criterion for neutrino selection
    myAnalysis->setUsedNeutrinoSelectionForTopPairReconstruction(NeutrinoSelectionCriterion::chi2);

    // Set method of ttbar reconstruction
    myAnalysis->setUsedTTbarReconstructionCriterion(TTbarReconstructionCriterion::TopMassDifference);

		// **************************************************************
    // ****** End of set-up section ***********
		// **************************************************************


		// **************************************************************
    // ***** Input section:  Set location of input ntuples
		// **************************************************************

    // DATA
    myAnalysis->addInputFile("/data/data/ElectronHad_Run2011A-PromptReco-v2_AOD_CMSSW_4_1_4_Ov1.2_LQnTuple*.root");

    // MC samples
    // TTJets, TToBLNu_t-channel, TToBLnu_tw-channel, WJets, DYJets
    myAnalysis->addInputFile(
			"/data/mc//TTJets_TuneD6T_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_AODSIM_CMSSW_4_1_4_Ov1.1_LQnTuple*.root");
    // QCD BC to E:  20-30, 30-80, 80-170
    // QCD EM Enriched:  20-30, 30-80, 80-170
    // QCD Photon + Jets:  40-100, 100-200, 200

		// **************************************************************
    // ***** End of input section
		// **************************************************************

    cout << "starting analysis" << endl;

    myAnalysis->analyze();	// Do analysis:  See below

    watch.Stop();
    watch.Print();

    return 0;
}


void IntroAnalysis::analyze() {
    createHistograms();
    cout << "detected samples:" << endl;
    for (unsigned int sample = 0; sample < DataType::NUMBER_OF_DATA_TYPES; ++sample) {
        if (eventReader->getSeenDatatypes()[sample])
            cout << DataType::names[sample] << endl;
    }

		// **************************************************************
    // ***** Main event loop
		// **************************************************************
    while (eventReader->hasNextEvent()) {
        initiateEvent();
        printNumberOfProccessedEventsEvery(100000);
        inspectEvents();
        doTTbarCutFlow(); // Counts events that pass each cut
        doDiElectronAnalysis();  // Special additional analysis

        doTTBarAnalysis();  // Main analaysis function -- defined below

        // Additional analyses
        doNotePlots();
        doQCDStudy();
        doJetAnalysis();
        if (currentEvent.getDataType() == DataType::ttbar)
        	doMCttbarReconstruction();

    }
    printInterestingEvents();
    printSummary();
}

void IntroAnalysis::initiateEvent() {
    currentEvent = eventReader->getNextEvent();

    if( !isGood(currentEvent.runnumber(), currentEvent.lumiblock()) ) return;

    // **************************************************************
    // ***** Instantiate reconstructed candidates
		// **************************************************************
    ttbarCandidate = TopPairCandidateTutorial(currentEvent);
    tPrimeCandidate = ToplikeCandidate(currentEvent);
    loneTopsNoMassConstr = TopNoMassConstraint(currentEvent);

    weight = weights.getWeight(currentEvent.getDataType());
    if(!currentEvent.isRealData()){
        weight *= weights.reweightPileUp(currentEvent.numberOfGeneratedPileUpVertices());
    }

    histMan.setCurrentDataType(ttbarCandidate.getDataType());
    histMan.setCurrentJetBin(currentEvent.GoodJets().size());
    histMan.setCurrentBJetBin(currentEvent.GoodBJets().size());
}


// **************************************************************
// ***** Perform main ttbar analysis
// **************************************************************


// **************************************************************
// ***** Select events
// **************************************************************

// *** Main event selection function.  See below it for functions for each
// selection step.

bool TopPairCandidateTutorial::passesSelectionStep(enum TTbarEPlusJetsSelection::Step step) const {
    switch (step) {
    case TTbarEPlusJetsSelection::FilterOutScraping:
        return passesScrapingFilter();
    case TTbarEPlusJetsSelection::HighLevelTrigger:
        return passesHighLevelTrigger();
    case TTbarEPlusJetsSelection::GoodPrimaryvertex:
        return hasOneGoodPrimaryVertex();
    case TTbarEPlusJetsSelection::OneIsolatedElectron:
        return hasOnlyOneGoodIsolatedElectron();
    case TTbarEPlusJetsSelection::ConversionRejection:
        return isolatedElectronDoesNotComeFromConversion();
    case TTbarEPlusJetsSelection::ConversionFinder:
        return isolatedElectronNotTaggedAsFromConversion();
    case TTbarEPlusJetsSelection::LooseMuonVeto:
        return hasNoIsolatedMuon();
    case TTbarEPlusJetsSelection::AtLeastOneGoodJets:
        return hasAtLeastOneGoodJet();
    case TTbarEPlusJetsSelection::AtLeastTwoGoodJets:
        return hasAtLeastTwoGoodJets();
    case TTbarEPlusJetsSelection::AtLeastThreeGoodJets:
        return hasAtLeastThreeGoodJets();
    case TTbarEPlusJetsSelection::AtLeastFourGoodJets:
        return hasAtLeastFourGoodJets();
    case TTbarEPlusJetsSelection::Zveto:
        return isNotAZBosonEvent();
    default:
        return false;
    }
}

/*
For details about the following selection steps, search as indicated
through the code reference file.

hasOneGoodPrimaryVertex(); -- See Event::PrimaryVertex() and Vertex::isGood()

hasOnlyOneGoodIsolatedElectron(); -- See Electron::isGood() and
	Electron::relativeIsolation()

isolatedElectronDoesNotComeFromConversion(); -- See Electron::isFromConversion()

isolatedElectronNotTaggedAsFromConversion(); --
	See Electron::isTaggedAsConversion()

hasNoIsolatedMuon(); -- See Muon::relativeIsolation()

hasAtLeastOneGoodJet();
hasAtLeastTwoGoodJets();
hasAtLeastThreeGoodJets();
hasAtLeastFourGoodJets(); -- See Jet::isGood()
*/


bool TopPairCandidateTutorial::passesScrapingFilter() const {
    if (tracks.size() > 10) {
        if (numberOfHighPurityTracks / (1.0 * tracks.size()) > 0.25)
            return true;
        else
            return false;
    } else
        return isBeamScraping == false;
}


bool TopPairCandidateTutorial::passesHighLevelTrigger() const {
    if (isRealData()) {
        if (runNumber < 140041)
            return HLT(HLTriggers::HLT_Ele10_LW_L1R);
        else if (runNumber >= 140041 && runNumber <= 143962)
            return HLT(HLTriggers::HLT_Ele15_SW_L1R);
        else if (runNumber > 143962 && runNumber <= 146427)
            return HLT(HLTriggers::HLT_Ele15_SW_CaloEleId_L1R);
        else if (runNumber > 146427 && runNumber <= 147116)
            return HLT(HLTriggers::HLT_Ele17_SW_CaloEleId_L1R);
        else if (runNumber > 147116 && runNumber <= 148818)
            return HLT(HLTriggers::HLT_Ele17_SW_TightEleId_L1R);
        else if (runNumber >= 148819 && runNumber < 149181)
            return HLT(HLTriggers::HLT_Ele22_SW_TighterEleId_L1R);
            //return HLT(HLTriggers::HLT_Ele22_SW_TighterEleId_L1R_v2);
        else if(runNumber >= 149181 && runNumber < 160000)
            return HLT(HLTriggers::HLT_Ele22_SW_TighterEleId_L1R);
            //return HLT(HLTriggers::HLT_Ele22_SW_TighterEleId_L1R_v3);
        else if(runNumber > 160000)
            return HLT(HLTriggers::HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_TrkIdT_QuadCentralJet30) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_QuadCentralJet30) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP) || HLT(
                    HLTriggers::HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20);
        else
            return false;
    }
    else // do not use HLT for MC
        return true;
}


bool TopPairCandidateTutorial::isNotAZBosonEvent() const {
    float invariantMass = 0;
    bool isZEvent = false;
    ElectronPointer isoElectron;

    if (Event::usePFIsolation && goodPFIsolatedElectrons.size() > 0)
        isoElectron = goodPFIsolatedElectrons.front();
    else if (goodIsolatedElectrons.size() > 0)
        isoElectron = goodIsolatedElectrons.front();

    if (isoElectron != NULL && allElectrons.size() > 1) {
        for (unsigned int index = 0; index < allElectrons.size(); ++index) {
            const ElectronPointer looseElectron = allElectrons.at(index);
            bool passLooseIso = false;

            if (Event::usePFIsolation)
                passLooseIso = looseElectron->isLoose() && looseElectron->pfIsolation() < 1.;
            else
                passLooseIso = looseElectron->isLoose() && looseElectron->relativeIsolation() < 1.;

            if (passLooseIso)
                invariantMass = isoElectron->invariantMass(looseElectron);
            else
                invariantMass = 0;

            bool passesLowerLimit = invariantMass > 76;
            bool passesUpperLimit = invariantMass < 106;
            if (passesLowerLimit && passesUpperLimit)
                isZEvent = true;
        }

    }


    return isZEvent == false;
}


// The following three functions are just interfaces to the main selection
// function above.

bool TopPairCandidateTutorial::passesFullTTbarEPlusJetSelection() const {
    unsigned int newstep = (int) TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS - 1;
    return passesSelectionStepUpTo((TTbarEPlusJetsSelection::Step) newstep);
}

bool TopPairCandidateTutorial::passesSelectionStepUpTo(enum TTbarEPlusJetsSelection::Step step) const {
    if (step == TTbarEPlusJetsSelection::FilterOutScraping)
        return passesSelectionStep(step);
    else {
        unsigned int newstep = (int) step - 1;
        return passesSelectionStep(step) && passesSelectionStepUpTo((TTbarEPlusJetsSelection::Step) newstep);
    }
}

bool TopPairCandidateTutorial::passesNMinus1(enum TTbarEPlusJetsSelection::Step omitted) const {
    bool passes(true);

    for (unsigned int cut = 0; cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
        if (cut == (unsigned int) omitted)
            continue;
        passes = passes && passesSelectionStep((TTbarEPlusJetsSelection::Step) cut);
    }
    return passes;
}

// **************************************************************
// ***** End of event selection section
// **************************************************************

// **************************************************************
// ***** ttbar candidate reconstruction
// **************************************************************

void TopPairCandidateTutorial::reconstructTTbar(ElectronPointer electron) {
    if (goodJets.size() < 4)
      throw ReconstructionException("Not enough jets available to reconstruct top event using Mass Equality method.");
    electronFromW = electron;
    selectedNeutrino = 0;
    currentSelectedNeutrino = 0;
    reconstructNeutrinos();
    double chosen_TopMassDifference(9999999.);
    double chosen_Chi2Total(9999999.);

    for (unsigned short hadBindex = 0; hadBindex < goodJets.size(); ++hadBindex) {
            for (unsigned short lepBindex = 0; lepBindex < goodJets.size(); ++lepBindex) {
                if (lepBindex == hadBindex)
                    continue;
                for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
                    if (jet1Index == lepBindex || jet1Index == hadBindex)
                        continue;
                    for (unsigned short jet2Index = 0; jet2Index < goodJets.size(); ++jet2Index) {
                        if (jet2Index == jet1Index || jet2Index == lepBindex || jet2Index == hadBindex)
                            continue;
                        hadronicBJet = goodJets.at(hadBindex);
                        leptonicBJet = goodJets.at(lepBindex);
                        jet1FromW = goodJets.at(jet1Index);
                        jet2FromW = goodJets.at(jet2Index);

                        leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *electronFromW));
                        leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *electronFromW));
                        hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
                        leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
                        leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
                        hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
                        fillHypotheses();
                        selectNeutrinoSolution();
                        double TopMassDifference = calculateTopMassDifference(currentSelectedNeutrino);
                        double chi2 = getTotalChi2(currentSelectedNeutrino);
                        switch (usedTTbarReconstruction) {
                        case TTbarReconstructionCriterion::TopMassDifference:
                        	if (TopMassDifference < chosen_TopMassDifference) {
                        		hadronicBIndex = hadBindex;
                        		leptonicBIndex = lepBindex;
                        		jet1FromWIndex = jet1Index;
                        		jet2FromWIndex = jet2Index;
                        		chosen_TopMassDifference = TopMassDifference;
                        		selectedNeutrino = currentSelectedNeutrino;
                        	}
                        	break;

                        case TTbarReconstructionCriterion::chi2:
                        	if (chi2 < chosen_Chi2Total) {
                        		hadronicBIndex = hadBindex;
                        		leptonicBIndex = lepBindex;
                        		jet1FromWIndex = jet1Index;
                        		jet2FromWIndex = jet2Index;
                        		chosen_Chi2Total = chi2;
                        		selectedNeutrino = currentSelectedNeutrino;
                        	}
                        	break;
                        }
                    }
                }
            }
	}
    std::sort(solutions.begin(), solutions.end(), compareSolutions);
    hadronicBJet = goodJets.at(hadronicBIndex);
    leptonicBJet = goodJets.at(leptonicBIndex);
    jet1FromW = goodJets.at(jet1FromWIndex);
    jet2FromW = goodJets.at(jet2FromWIndex);
    leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *electronFromW));
    leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *electronFromW));
    hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
    leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
    leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
    hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
    if (selectedNeutrino == 1)
    	ttbarResonance = ParticlePointer(new Particle(*leptonicTop1 + *hadronicTop));
    else
    	ttbarResonance = ParticlePointer(new Particle(*leptonicTop2 + *hadronicTop));
    doneReconstruction = true;
}


void TopPairCandidateTutorial::reconstructNeutrinos() {
    boost::array<double, 2> neutrinoPzs = computeNeutrinoPz();
    double energy1 = sqrt(met->et() * met->et() + neutrinoPzs.at(0) * neutrinoPzs.at(0));
    double energy2 = sqrt(met->et() * met->et() + neutrinoPzs.at(1) * neutrinoPzs.at(1));
    neutrino1 = ParticlePointer(new Particle(energy1, met->px(), met->py(), neutrinoPzs.at(0)));
    neutrino2 = ParticlePointer(new Particle(energy2, met->px(), met->py(), neutrinoPzs.at(1)));

    if (isnan(neutrino1->energy()) && isnan(neutrino2->energy()))
        throw ReconstructionException("No physical neutrino solution found");
    else if (isnan(neutrino1->energy()))
        neutrino1 = neutrino2;
    else if (isnan(neutrino2->energy()))
        neutrino2 = neutrino1;
}

const boost::array<double, 2> TopPairCandidateTutorial::computeNeutrinoPz() {
    if (electronFromW == 0)
        throw ReconstructionException("Could not reconstruct neutrinos: no isolated electrons found");
    if (met->energy() == 0)
        throw ReconstructionException("Could not reconstruct neutrinos: no MET found");
    boost::array<double, 2> neutrinoPzs;
    //    const ElectronPointer electron = goodIsolatedElectrons.front();

    double pz1(0), pz2(0);
    //    double M_W = 80.389;
    double M_e = 0.0005;
    double ee = electronFromW->energy();
    double pxe = electronFromW->px();
    double pye = electronFromW->py();
    double pze = electronFromW->pz();
    double pxnu = met->px();
    double pynu = met->py();

    double a = W_mass * W_mass - M_e * M_e + 2.0 * pxe * pxnu + 2.0 * pye * pynu;
    double A = 4.0 * (ee * ee - pze * pze);
    double B = -4.0 * a * pze;
    double C = 4.0 * ee * ee * (pxnu * pxnu + pynu * pynu) - a * a;

    double tmproot = B * B - 4.0 * A * C;
    if (tmproot < 0) {
        pz1 = pz2 = -B / (2 * A);
    } else {
        pz1 = (-B + TMath::Sqrt(tmproot)) / (2.0 * A);
        pz2 = (-B - TMath::Sqrt(tmproot)) / (2.0 * A);

    }
    neutrinoPzs[0] = pz1;
    neutrinoPzs[1] = pz2;
    return neutrinoPzs;
}


void TopPairCandidateTutorial::fillHypotheses() {
	TtbarHypothesisPointer hypothesis1(fillHypothesis(1));
	TtbarHypothesisPointer hypothesis2(fillHypothesis(2));
	solutions.push_back(hypothesis1);
	solutions.push_back(hypothesis2);

}

const TtbarHypothesisPointer TopPairCandidateTutorial::fillHypothesis(unsigned short int neutrinoSolution) {
	TtbarHypothesisPointer hypothesis(new TtbarHypothesis());
	hypothesis->electronFromW = electronFromW;
	hypothesis->leptonicBjet = leptonicBJet;
	hypothesis->hadronicBJet = hadronicBJet;
	hypothesis->jet1FromW = jet1FromW;
	hypothesis->jet2FromW = jet2FromW;
	hypothesis->hadronicW = hadronicW;
	hypothesis->hadronicTop = hadronicTop;
	hypothesis->hadronicChi2 = getHadronicChi2();
	if(neutrinoSolution == 1) {
		hypothesis->neutrinoFromW = neutrino1;
		hypothesis->leptonicW = leptonicW1;
		hypothesis->leptonicTop = leptonicTop1;
	}
	else {
		hypothesis->neutrinoFromW = neutrino2;
		hypothesis->leptonicW = leptonicW2;
		hypothesis->leptonicTop = leptonicTop2;
	}

	hypothesis->totalChi2 = getTotalChi2(neutrinoSolution);
	hypothesis->globalChi2 = getGlobalChi2(neutrinoSolution);
	hypothesis->leptonicChi2 = getLeptonicChi2(neutrinoSolution);
	ParticlePointer resonance(new Particle(*hypothesis->leptonicTop + *hypothesis->hadronicTop));
	hypothesis->resonance = resonance;
	return hypothesis;
}

void TopPairCandidateTutorial::selectNeutrinoSolution() {

    if (leptonicTop1->mass() < 0 && leptonicTop2->mass() < 0) {
        inspectReconstructedEvent();
        throw ReconstructionException("No valid neutrino solution found");
    } else if (leptonicTop1->mass() < 0 && leptonicTop2->mass() > 0) {
        currentSelectedNeutrino = 2;
    } else if (leptonicTop1->mass() > 0 && leptonicTop2->mass() < 0) {
        currentSelectedNeutrino = 1;
    } else {// both solutions give positive mass
        switch (usedNeutrinoSelection) {
        case NeutrinoSelectionCriterion::TopMassDifference:
            fabs(leptonicTop1->mass()-hadronicTop->mass()) < fabs(leptonicTop2->mass()-hadronicTop->mass()) ?
                      currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            break;
        case NeutrinoSelectionCriterion::chi2:
            getTotalChi2(1) < getTotalChi2(2) ? currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            break;

        case NeutrinoSelectionCriterion::pzClosestToLepton:
            fabs(neutrino1->pz() - electronFromW->pz()) < fabs(neutrino2->pz()
                    - electronFromW->pz()) ? currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            break;

        case NeutrinoSelectionCriterion::mostCentral:
            fabs(neutrino1->pz()) < fabs(neutrino2->pz()) ? currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            break;

        case NeutrinoSelectionCriterion::pzClosestToLeptonOrMostcentralIfAbove300:
            fabs(neutrino1->pz() - electronFromW->pz()) < fabs(neutrino2->pz()
                    - electronFromW->pz()) ? currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            if (fabs(neutrino1->pz()) > 300 || fabs(neutrino2->pz()) > 300)
                fabs(neutrino1->pz()) < fabs(neutrino2->pz()) ? currentSelectedNeutrino = 1 : currentSelectedNeutrino
                        = 2;
            break;

        case NeutrinoSelectionCriterion::largestValueOfCosine:
            TVector3 p3W, p3e;
            //TODO clean up
            p3W = leptonicW1->getFourVector().Vect();
            p3e = electronFromW->getFourVector().Vect();

            double sinthcm1 = 2. * (p3e.Perp(p3W)) / W_mass;
            p3W = leptonicW2->getFourVector().Vect();
            double sinthcm2 = 2. * (p3e.Perp(p3W)) / W_mass;

            double costhcm1 = TMath::Sqrt(1. - sinthcm1 * sinthcm1);
            double costhcm2 = TMath::Sqrt(1. - sinthcm2 * sinthcm2);
            costhcm1 > costhcm2 ? currentSelectedNeutrino = 1 : currentSelectedNeutrino = 2;
            break;

        }
    }

}

double TopPairCandidateTutorial::calculateTopMassDifference(unsigned short int neutrinoSolution) const {

  double LeptonicTop1MassDifference = fabs(leptonicTop1->mass()-hadronicTop->mass());
  double LeptonicTop2MassDifference = fabs(leptonicTop2->mass()-hadronicTop->mass());

  if (neutrinoSolution == 1)
    return LeptonicTop1MassDifference;
  else
    return LeptonicTop2MassDifference;

}

double TopPairCandidateTutorial::getTotalChi2() {
    double totalChi2(9999999);
    double firstTotalChi2 = getTotalChi2(1);
    double secondTotalChi2 = getTotalChi2(2);
    selectedNeutrino == 1 ? totalChi2 = firstTotalChi2 : totalChi2 = secondTotalChi2;
    return totalChi2;
}

double TopPairCandidateTutorial::getTotalChi2(unsigned short int neutrinoSolution) const {
    return getLeptonicChi2(neutrinoSolution) + getHadronicChi2() + getGlobalChi2(neutrinoSolution);
}

double TopPairCandidateTutorial::getLeptonicChi2(unsigned short int neutrinoSolution) const {
    double topMass(0);
    double angle = leptonicBJet->angle(electronFromW);
    if (neutrinoSolution == 1)
        topMass = leptonicTop1->mass();
    else
        topMass = leptonicTop2->mass();

    return getLeptonicChi2(topMass, angle);
}

// ********
// Reference values from MC used to check correctness of reconstruction
double const TopPairCandidateTutorial::matched_angle = 0.945666;
double const TopPairCandidateTutorial::matched_angle_sigma = 0.311091;
double const TopPairCandidateTutorial::matched_leptonic_top_mass = 169.0;
double const TopPairCandidateTutorial::matched_leptonic_top_mass_sigma = 16.3;
double const TopPairCandidateTutorial::matched_hadronic_W_mass = 83.;
double const TopPairCandidateTutorial::matched_hadronic_W_mass_sigma = 10.8995;
double const TopPairCandidateTutorial::matched_hadronic_top_mass = 174.7;
double const TopPairCandidateTutorial::matched_hadronic_top_mass_sigma = 14.6;
double const TopPairCandidateTutorial::matched_ptratio = 0.18552;
double const TopPairCandidateTutorial::matched_ptratio_sigma = 0.401973;
double const TopPairCandidateTutorial::matched_pt_ttbarSystem = 0.;
double const TopPairCandidateTutorial::matched_pt_ttbarSystem_sigma = 50.;
double const TopPairCandidateTutorial::matched_HTSystem = 1;
double const TopPairCandidateTutorial::matched_HTSystem_sigma = 0.1;
double const TopPairCandidateTutorial::W_mass = 80.389;


double TopPairCandidateTutorial::getLeptonicChi2(double topMass, double angle) const {
    double massDifference = TMath::Power(topMass - matched_leptonic_top_mass, 2);
    double massError = 2 * matched_leptonic_top_mass_sigma * matched_leptonic_top_mass_sigma;
    double massTerm = massDifference / massError;

    double angleDifference = TMath::Power(angle - matched_angle, 2);
    double angleError = 2 * matched_angle_sigma * matched_angle_sigma;
    double angleTerm = angleDifference / angleError;
    return 1 / sqrt(2) * (angleTerm + massTerm);
}

double TopPairCandidateTutorial::getHadronicChi2() const {
    double ptRatioDifference = TMath::Power(PtRatio() - matched_ptratio, 2);
    double ptRatioError = 2 * matched_ptratio_sigma * matched_ptratio_sigma;
    double ptRatioTerm = ptRatioDifference / ptRatioError;

    double WmassDifference = TMath::Power(hadronicW->mass() - matched_hadronic_W_mass, 2);
    double WmassError = 2 * matched_hadronic_W_mass_sigma * matched_hadronic_W_mass_sigma;
    double WmassTerm = WmassDifference / WmassError;

    double topMassDifference = TMath::Power(hadronicTop->mass() - matched_hadronic_top_mass, 2);
    double topMassError = 2 * matched_hadronic_top_mass_sigma * matched_hadronic_top_mass_sigma;
    double topMassTerm = topMassDifference / topMassError;
    return 1 / sqrt(3) * (topMassTerm + WmassTerm + ptRatioTerm);
    return 0;
}

double TopPairCandidateTutorial::PtRatio() const {
    return TMath::Log(hadronicTop->pt() / hadronicW->pt());
}

double TopPairCandidateTutorial::getGlobalChi2(unsigned short neutrinoSolution) const {
    double pttbar = PtTtbarSystem(neutrinoSolution);
    double pttbarDifference = TMath::Power(pttbar - matched_pt_ttbarSystem, 2);
    double pttbarError = (2 * matched_pt_ttbarSystem_sigma * matched_pt_ttbarSystem_sigma);
    double pttbarTerm = pttbarDifference / pttbarError;

    double htSystemDifference = TMath::Power(HTSystem() - matched_HTSystem, 2);
    double htSystemError = matched_HTSystem_sigma * matched_HTSystem_sigma * 2;
    double htSystemTerm = htSystemDifference / htSystemError;
    return 1 / sqrt(2) * (pttbarTerm + htSystemTerm);
}

double TopPairCandidateTutorial::PtTtbarSystem(unsigned short neutrinoSolution) const {
    ParticlePointer combined;
    if (neutrinoSolution == 1)
        combined = ParticlePointer(new Particle(*leptonicTop1 + *hadronicTop));
    else
        combined = ParticlePointer(new Particle(*leptonicTop2 + *hadronicTop));
    return combined->pt();
}

double TopPairCandidateTutorial::HT(unsigned short jetLimit) const {
    double HT(0);
    unsigned short limit = goodJets.size();
    if (limit > jetLimit + 1)
        limit = jetLimit + 1;

    for (unsigned short index = 0; index < limit; ++index)
        HT += goodJets.at(index)->pt();

    return HT;
}

double TopPairCandidateTutorial::HTSystem() const {
    return sumPt() / HT(8);
}

double TopPairCandidateTutorial::sumPt() const {
    return leptonicBJet->pt() + hadronicBJet->pt() + jet1FromW->pt() + jet2FromW->pt();
}


// **************************************************************
// ***** End ttbar candidate reconstruction section
// **************************************************************

// Analysis function that calls event selection and reconstruction functions
void IntroAnalysis::doTTBarAnalysis() {
    if (ttbarCandidate.passesNMinus1(TTbarEPlusJetsSelection::AtLeastFourGoodJets)) {
        histMan.H1D("numberOfJets")->Fill(ttbarCandidate.GoodJets().size());
    }

    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex)
            && ttbarCandidate.hasNoIsolatedMuon() && ttbarCandidate.isNotAZBosonEvent()
            && ttbarCandidate.hasAtLeastFourGoodJets()) {
        const ElectronCollection coll = ttbarCandidate.Electrons();
        for (unsigned int index = 0; index < coll.size(); ++index) {
            const ElectronPointer electron = coll.at(index);
            bool passesEt = electron->et() > 30;
            bool passesEta = fabs(electron->superClusterEta()) < 2.5
                    && electron->isInCrack() == false;
            bool passesID = electron->VBTF_W70_ElectronID();
            bool noConversion = electron->isFromConversion() == false;
            if (passesEt && passesEta && passesID && noConversion) {
                histMan.H1D("electronD0")->Fill(electron->d0(), weight);
            }
        }
    }
    // ******* Main event selection call -- See selection section above
    if (ttbarCandidate.passesFullTTbarEPlusJetSelection()) {
        histMan.H1D("numberOfBJets")->Fill(ttbarCandidate.GoodBJets().size(), weight);

        // ******** Now reconstruct candidates --
        // See recontruction section above
        if(Event::usePFIsolation) {
            ttbarCandidate.reconstructTTbar(ttbarCandidate.GoodPFIsolatedElectrons().front());

            // The following candidates use the ttbar selection.
            // See src/ToplikeCandidate.cpp for function definitions.
            tPrimeCandidate.recoTprimeUsingChi2(ttbarCandidate.GoodPFIsolatedElectrons().front());
            tPrimeCandidate.recoBestSingleTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
            loneTopsNoMassConstr.recoBestSingleTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
        } else {
            ttbarCandidate.reconstructTTbar(ttbarCandidate.GoodIsolatedElectrons().front());

            // The following candidates use the ttbar selection.
            // See src/ToplikeCandidate.cpp for function definitions.
            tPrimeCandidate.recoTprimeUsingChi2(ttbarCandidate.GoodIsolatedElectrons().front());
            tPrimeCandidate.recoBestSingleTop(ttbarCandidate.GoodIsolatedElectrons().front());
            loneTopsNoMassConstr.recoBestSingleTop(ttbarCandidate.GoodIsolatedElectrons().front());
			  }
        vector<TtbarHypothesisPointer> solutions = ttbarCandidate.Solutions();
        const ParticlePointer resonance = ttbarCandidate.getResonance();
        double mttbar = ttbarCandidate.mttbar();

        ParticlePointer leadingTop, nextToLeadingTop;

        if (ttbarCandidate.getHadronicTop()->pt() > ttbarCandidate.getLeptonicTop()->pt()) {
            leadingTop = ttbarCandidate.getHadronicTop();
            nextToLeadingTop = ttbarCandidate.getLeptonicTop();
        } else {
            leadingTop = ttbarCandidate.getLeptonicTop();
            nextToLeadingTop = ttbarCandidate.getHadronicTop();
        }
        double angleTops = leadingTop->angle(nextToLeadingTop);

// **************************************************************
// ***** Fill histograms
// **************************************************************

        histMan.H1D_BJetBinned("angleTops")->Fill(angleTops, weight);
        histMan.H2D_BJetBinned("angleTops_vs_mttbar")->Fill(mttbar, angleTops, weight);
        histMan.H1D_BJetBinned("mLeptonicTop")->Fill(ttbarCandidate.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mHadronicTop")->Fill(ttbarCandidate.getHadronicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mAllTop")->Fill(ttbarCandidate.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mAllTop")->Fill(ttbarCandidate.getHadronicTop()->mass(), weight);
        histMan.H1D_BJetBinned("chiHadronicTop")->Fill(solutions.at(0)->hadronicChi2,
        	weight);
        histMan.H1D_BJetBinned("chiLeptonicTop")->Fill(solutions.at(0)->leptonicChi2,
        	weight);
        histMan.H1D_BJetBinned("chiGlobal")->Fill(solutions.at(0)->globalChi2,
        	weight);
        histMan.H1D_BJetBinned("chiTotal")->Fill(solutions.at(0)->totalChi2,
        	weight);
        histMan.H1D_BJetBinned("mLeptonicTopLone")->Fill(tPrimeCandidate.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mHadronicTopLone")->Fill(tPrimeCandidate.getHadronicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mLepTopLoneNM")->Fill(loneTopsNoMassConstr.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mHadTopLoneNM")->Fill(loneTopsNoMassConstr.getHadronicTop()->mass(), weight);

        histMan.H1D_BJetBinned("MET")->Fill(ttbarCandidate.MET()->et(), weight);
        histMan.H2D_BJetBinned("METvsMttbar")->Fill(mttbar, ttbarCandidate.MET()->et(), weight);
        histMan.H1D_BJetBinned("HT")->Fill(ttbarCandidate.fullHT(), weight);
        histMan.H2D_BJetBinned("HTvsMttbar")->Fill(mttbar, ttbarCandidate.fullHT(), weight);
				histMan.H1D_BJetBinned("leadingJetMass")->Fill(ttbarCandidate.GoodJets().front()->mass(), weight);
        if (ttbarCandidate.isRealData() == false) {
						histMan.H1D_BJetBinned("leadingGenJetMass")->Fill(ttbarCandidate.GenJets().front()->mass(), weight);
						topTruth leptTop = tPrimeCandidate.getMCMatches(tPrimeCandidate.getLeptonicTop());
						topTruth hadTop = tPrimeCandidate.getMCMatches(tPrimeCandidate.getHadronicTop());
						histMan.H1D_BJetBinned("lepTopDeltaR")->Fill(leptTop.deltaR, weight);
						histMan.H1D_BJetBinned("hadTopDeltaR")->Fill(hadTop.deltaR, weight);
						histMan.H1D_BJetBinned("lepTopLep")->Fill(leptTop.leptonic, weight);
						histMan.H1D_BJetBinned("hadTopHad")->Fill(hadTop.hadronic, weight);
						histMan.H1D_BJetBinned("numHadTopMCMatches")->Fill(tPrimeCandidate.getNumMCMatchesHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopMCMatches")->Fill(tPrimeCandidate.getNumMCMatchesLTop(), weight);
						histMan.H1D_BJetBinned("numHadTopCorrectID")->Fill(tPrimeCandidate.getNumCorrectIDHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopCorrectID")->Fill(tPrimeCandidate.getNumCorrectIDLTop(), weight);
						topTruth leptTopNM = loneTopsNoMassConstr.getMCMatches(loneTopsNoMassConstr.getLeptonicTop());
						topTruth hadTopNM = loneTopsNoMassConstr.getMCMatches(loneTopsNoMassConstr.getHadronicTop());
						histMan.H1D_BJetBinned("lepTopDeltaRNM")->Fill(leptTopNM.deltaR, weight);
						histMan.H1D_BJetBinned("hadTopDeltaRNM")->Fill(hadTopNM.deltaR, weight);
						histMan.H1D_BJetBinned("numHadTopMCMatchesNM")->Fill(loneTopsNoMassConstr.getNumMCMatchesHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopMCMatchesNM")->Fill(loneTopsNoMassConstr.getNumMCMatchesLTop(), weight);
						histMan.H1D_BJetBinned("numHadTopCorrectIDNM")->Fill(loneTopsNoMassConstr.getNumCorrectIDHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopCorrectIDNM")->Fill(loneTopsNoMassConstr.getNumCorrectIDLTop(), weight);
				}

        if (Event::usePFIsolation)
            histMan.H1D_BJetBinned("mtW")->Fill(ttbarCandidate.transverseWmass(
                    ttbarCandidate.GoodPFIsolatedElectrons().front()), weight);
        else
            histMan.H1D_BJetBinned("mtW")->Fill(ttbarCandidate.transverseWmass(
                    ttbarCandidate.GoodIsolatedElectrons().front()), weight);
        histMan.H1D_BJetBinned("m3")->Fill(ttbarCandidate.M3(), weight);


        histMan.H1D_BJetBinned("mttbar")->Fill(mttbar, weight);
				/*
        histMan.H1D_BJetBinned("tPrimeMass")->Fill(tPrimeCandidate.tpmass(), weight);
        histMan.H1D_BJetBinned("tPrimeHT")->Fill(tPrimeCandidate.TPrimeHTSystem(), weight);
        histMan.H1D_BJetBinned("tPrimepT")->Fill(tPrimeCandidate.PtTPrimeSystem(), weight);
        histMan.H1D_BJetBinned("tPrime_pt")->Fill((tPrimeCandidate.getResonance())->pt(), weight);
        histMan.H1D_BJetBinned("tPrime_px")->Fill((tPrimeCandidate.getResonance())->px(), weight);
        histMan.H1D_BJetBinned("tPrime_py")->Fill((tPrimeCandidate.getResonance())->py(), weight);
        histMan.H1D_BJetBinned("tPrime_pz")->Fill((tPrimeCandidate.getResonance())->pz(), weight);
        histMan.H2D_BJetBinned("tPrime_pt_vs_tPrimeM")->Fill(tPrimeCandidate.tpmass(), (tPrimeCandidate.getResonance())->pt(), weight);
				*/

        histMan.H1D_BJetBinned("ttbar_pt")->Fill(resonance->pt(), weight);

        histMan.H1D_BJetBinned("ttbar_px")->Fill(resonance->px(), weight);
        histMan.H1D_BJetBinned("ttbar_py")->Fill(resonance->py(), weight);
        histMan.H1D_BJetBinned("ttbar_pz")->Fill(resonance->pz(), weight);

        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar")->Fill(mttbar, resonance->pt(), weight);


        histMan.H1D("electron_et")->Fill(ttbarCandidate.getElectronFromWDecay()->et(), weight);
        histMan.H1D_BJetBinned("neutrino_pz")->Fill(ttbarCandidate.getNeutrinoFromWDecay()->pz(), weight);

        histMan.H1D_BJetBinned("pt_leadingTop")->Fill(leadingTop->pt(), weight);
        histMan.H1D_BJetBinned("pt_NextToLeadingTop")->Fill(nextToLeadingTop->pt(), weight);
        histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar")->Fill(mttbar, leadingTop->pt(), weight);
        histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar")->Fill(mttbar, nextToLeadingTop->pt(), weight);

        if (ttbarCandidate.MET()->pt() > 20) {
            histMan.H1D_BJetBinned("angleTops_withMETCut")->Fill(angleTops, weight);
            histMan.H2D_BJetBinned("angleTops_vs_mttbar_withMETCut")->Fill(mttbar, angleTops, weight);
            histMan.H1D_BJetBinned("mttbar_withMETCut")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_withMETCut")->Fill(resonance->pt(), weight);
            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withMETCut")->Fill(mttbar, resonance->pt(), weight);

            histMan.H1D_BJetBinned("pt_leadingTop_withMETCut")->Fill(leadingTop->pt(), weight);
            histMan.H1D_BJetBinned("pt_NextToLeadingTop_withMETCut")->Fill(nextToLeadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETCut")->Fill(mttbar, leadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETCut")->Fill(mttbar, nextToLeadingTop->pt(), weight);

            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("angleTops_withMETAndAsymJets")->Fill(angleTops, weight);
                histMan.H2D_BJetBinned("angleTops_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, angleTops, weight);
                histMan.H1D_BJetBinned("mttbar_withMETAndAsymJets")->Fill(mttbar, weight);
                histMan.H1D_BJetBinned("ttbar_pt_withMETAndAsymJets")->Fill(resonance->pt(), weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, resonance->pt(), weight);

                histMan.H1D_BJetBinned("pt_leadingTop_withMETAndAsymJets")->Fill(leadingTop->pt(), weight);
                histMan.H1D_BJetBinned("pt_NextToLeadingTop_withMETAndAsymJets")->Fill(nextToLeadingTop->pt(), weight);
                histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, leadingTop->pt(), weight);
                histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, nextToLeadingTop->pt(), weight);
            }
        }

        if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
            histMan.H1D_BJetBinned("angleTops_withAsymJetsCut")->Fill(angleTops, weight);
            histMan.H2D_BJetBinned("angleTops_vs_mttbar_withAsymJetsCut")->Fill(mttbar, angleTops, weight);
            histMan.H1D_BJetBinned("mttbar_withAsymJetsCut")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_withAsymJetsCut")->Fill(resonance->pt(), weight);

            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withAsymJetsCut")->Fill(mttbar, resonance->pt(), weight);

            histMan.H1D_BJetBinned("pt_leadingTop_withAsymJetsCut")->Fill(leadingTop->pt(), weight);
            histMan.H1D_BJetBinned("pt_NextToLeadingTop_withAsymJetsCut")->Fill(nextToLeadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withAsymJetsCut")->Fill(mttbar, leadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut")->Fill(mttbar, nextToLeadingTop->pt(), weight);
        }

        for (unsigned int solutionIndex = 0; solutionIndex < solutions.size(); ++solutionIndex) {
            histMan.H1D_BJetBinned("mttbar_allSolutions")->Fill(solutions.at(solutionIndex)->resonance->mass(), weight);
            histMan.H1D_BJetBinned("ttbar_pt_allSolutions")->Fill(solutions.at(solutionIndex)->resonance->pt(), weight);
            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions")->Fill(
                    solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                    weight);
            if (solutionIndex == 1) {
                histMan.H1D_BJetBinned("mttbar_2ndSolution")->Fill(solutions.at(solutionIndex)->resonance->mass(),
                        weight);

                histMan.H1D_BJetBinned("ttbar_pt_2ndSolution")->Fill(solutions.at(solutionIndex)->resonance->pt(),
                        weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);

            }
            if (solutionIndex == 2) {
                histMan.H1D_BJetBinned("mttbar_3rdSolution")->Fill(solutions.at(solutionIndex)->resonance->mass(),
                        weight);

                histMan.H1D_BJetBinned("ttbar_pt_3rdSolution")->Fill(solutions.at(solutionIndex)->resonance->pt(),
                        weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);
            }

            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), weight);
                histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->pt(), weight);
                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);
                if (solutionIndex == 1) {
                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_2ndSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }
                if (solutionIndex == 2) {
                    histMan.H1D_BJetBinned("mttbar_3rdSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }

                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_BJetBinned("mttbar_allSolutions_withMETAndAsymJets")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withMETAndAsymJets")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETAndAsymJets")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }
            }

            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("mttbar_allSolutions_withAsymJetsCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), weight);
                histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withAsymJetsCut")->Fill(
                        solutions.at(solutionIndex)->resonance->pt(), weight);
                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withAsymJetsCut")->Fill(solutions.at(
                        solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);

                if (ttbarCandidate.MET()->et() > 20) {

                    if (solutionIndex == 1) {
                        histMan.H1D_BJetBinned("mttbar_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), weight);

                        histMan.H1D_BJetBinned("ttbar_pt_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->pt(), weight);
                        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    }

                    if (solutionIndex == 2) {
                        histMan.H1D_BJetBinned("mttbar_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), weight);

                        histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->pt(), weight);
                        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    }
                }

                if (solutionIndex == 1) {
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withAsymJetsCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);

                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                }

                if (solutionIndex == 2) {
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withAsymJetsCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);

                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                }

            }
        }
        if (ttbarCandidate.MET()->et() < 20) {
            histMan.H1D_BJetBinned("mttbar_QCDEnriched")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_QCDEnriched")->Fill(resonance->pt());
        }
        if (mttbar != mttbar) {//isnan
            ttbarCandidate.inspectReconstructedEvent();
        }
        if (ttbarCandidate.isRealData()) {
            cout << "run " << ttbarCandidate.runnumber() << ", event " << ttbarCandidate.eventnumber() << ", lumi "
                    << ttbarCandidate.lumiblock();
            cout << ", top pair invariant mass = " << mttbar << " GeV" << endl;
            interestingEvents.push_back(IntroOAT::InterestingEvent(ttbarCandidate, eventReader->getCurrentFile()));

            if (resonance->pt() > 100) {
                cout << "top pair pt = " << resonance->pt() << " GeV" << endl;
                ttbarCandidate.inspect();
                ttbarCandidate.inspectReconstructedEvent();
            }
        }

    }
}


// **************************************************************
// ***** Histogram creation function
// ***** Every histogram has to be defined here before it can be
// ***** filled in analysis function above
// **************************************************************

void IntroAnalysis::createHistograms() {
    histMan.setCurrentLumi(IntroAnalysis::luminosity);
    histMan.prepareForSeenDataTypes(eventReader->getSeenDatatypes());

    // addH1D_BJetBinned parameters are (name, name, number of bins, start bin,
    // end bin).  Creates set of histograms for 0, 1, 2, 3, and 4 jets and
    // 0, 1, 2, 3, and 4 or more jets.

    //histograms for Jet study
    histMan.addH1D_BJetBinned("AllJetMass", "AllJetMass", 500, 0, 500);
    histMan.addH1D_BJetBinned("AllGoodJetMass", "AllGoodJetMass", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastOneJets", "GoodJetMass_atLeastOneJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastTwoJets", "GoodJetMass_atLeastTwoJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastThreeJets", "GoodJetMass_atLeastThreeJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastFourJets", "GoodJetMass_atLeastFourJets", 500, 0, 500);
    //MC histograms
    histMan.addH1D("deltaRElectron", "delta R between truth and reco electron", 100, 0, 0.2);
    histMan.addH1D("deltaRLeptonicBjet", "delta R between truth and reco b-jet on leptonic side", 100, 0, 0.5);
    histMan.addH1D("deltaRHadronicBjet", "delta R between truth and reco b-jet on hadronic side", 100, 0, 0.5);
    histMan.addH1D("deltaRjet1fromW", "delta R between truth and reco jet1 from W decay", 100, 0, 0.5);
    histMan.addH1D("deltaRjet2fromW", "delta R between truth and reco jet2 from W decay", 100, 0, 0.5);

    histMan.addH1D("deltaRjet1", "delta R between quark from W and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet2", "delta R between antiquark from W and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet3", "delta R between b quark from top and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet4", "delta R between b quark from antitop and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet_sum", "summarized delta R between partons and genJets", 100, 0, 0.5);

    histMan.addH2D("deltaR_genJets_partons", "delta R between genJets from W as opposed to partons", 100, 0, 5, 100, 0, 5);

    histMan.addH1D("W_inv_mass_from_truth_partons", "W inv. mass from truth partons", 100, 0, 120);
    histMan.addH1D("W_inv_mass_from_genJets", "W inv. mass from genJets", 100, 0, 120);
    histMan.addH1D("top_leptonic_inv_mass_from_truth", "Leptonic top inv. mass from truth partons", 100, 100, 220);
    histMan.addH1D("top_hadronic_inv_mass_from_truth", "Haronic top inv. mass from truth partons", 100, 100, 220);
    histMan.addH1D("m3_mc", "M3 for truth event", 500, 0, 500);
    histMan.addH1D("m3_diff", "M3 difference between truth and reco", 500, 0, 500);

    //end of MC histograms

    histMan.addH1D("electron_et", "electron_et", 500, 0, 500);
    histMan.addH1D_JetBinned("MostPFIsolatedElectron_dPhiIn", "MostPFIsolatedElectron_dPhiIn", 50, 0, 0.1);
    histMan.addH1D_JetBinned("MostPFIsolatedElectron_dEtaIn", "MostPFIsolatedElectron_dEtaIn", 50, 0, 0.02 );
    histMan.addH1D_JetBinned("diElectronMass", "diElectronMass", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("tPrimeMass", "tPrimeMass", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrimepT", "tPrimepT", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrimeHT", "tPrimeHT", 600, 0, 2.0);

    histMan.addH1D_BJetBinned("mttbar_conversions", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_QCDEnriched", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withMETCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withMETCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withMETCut", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withMETAndAsymJets", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withMETAndAsymJets", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withMETAndAsymJets", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withAsymJetsCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withAsymJetsCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withAsymJetsCut", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_withMETCut", "mttbar_withMETCut", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_withMETAndAsymJets", "mttbar_withMETAndAsymJets", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_withAsymJetsCut", "mttbar_withAsymJetsCut", 5000, 0, 5000);
//
    histMan.addH1D_BJetBinned("mttbar_2ndSolution", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions", "mttbar_allSolutions", 5000, 0, 5000);
//
    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withMETCut", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withMETCut", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withMETCut", "mttbar_allSolutions", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withMETAndAsymJets", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withMETAndAsymJets", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withMETAndAsymJets", "mttbar_allSolutions", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withAsymJetsCut", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withAsymJetsCut", "mttbar_allSolutions", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mLeptonicTop", "mLeptonicTop", 500, 0, 500);
    histMan.addH1D_BJetBinned("mHadronicTop", "mHadronicTop", 500, 0, 500);
    histMan.addH1D_BJetBinned("mAllTop", "mAllTop", 500, 0, 500);
    histMan.addH1D_BJetBinned("m3", "m3", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mLeptonicTopLone", "mLeptonicTopLone", 600, 0, 600);
    histMan.addH1D_BJetBinned("mHadronicTopLone", "mHadronicTopLone", 600, 0, 600);
    histMan.addH1D_BJetBinned("mLepTopLoneNM", "mLepTopLoneNM", 600, 0, 600);
    histMan.addH1D_BJetBinned("mHadTopLoneNM", "mHadTopLoneNM", 600, 0, 600);
    histMan.addH1D_BJetBinned("mAllTop", "mAllTop", 600, 0, 600);
    histMan.addH1D_BJetBinned("chiHadronicTop", "chiHadronicTop", 400, 0, 200);
    histMan.addH1D_BJetBinned("chiLeptonicTop", "chiLeptonicTop", 400, 0, 200);
    histMan.addH1D_BJetBinned("chiGlobal", "chiGlobal", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("chiTotal", "chiTotal", 6000, 0, 6000);

    histMan.addH1D_BJetBinned("lepTopDeltaR", "lepTopDeltaR", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("hadTopDeltaR", "hadTopDeltaR", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("lepTopLep", "lepTopLep", 3, 0, 3);
    histMan.addH1D_BJetBinned("hadTopHad", "hadTopHad", 3, 0, 3);
    histMan.addH1D_BJetBinned("numHadTopMCMatches", "numHadTopMCMatches", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopMCMatches", "numLepTopMCMatches", 6, 0, 6);
    histMan.addH1D_BJetBinned("numHadTopCorrectID", "numHadTopCorrectID", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopCorrectID", "numLepTopCorrectID", 6, 0, 6);
    histMan.addH1D_BJetBinned("lepTopDeltaRNM", "lepTopDeltaRNM", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("hadTopDeltaRNM", "hadTopDeltaRNM", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("numHadTopMCMatchesNM", "numHadTopMCMatchesNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopMCMatchesNM", "numLepTopMCMatchesNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numHadTopCorrectIDNM", "numHadTopCorrectIDNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopCorrectIDNM", "numLepTopCorrectIDNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("tPrime_pt", "tPrime_pt", 6000, 0, 6000);
    histMan.addH2D_BJetBinned("tPrime_pt_vs_tPrimeM", "tPrime_pt_vs_tPrimeM", 600, 0, 6000, 600, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_px", "tPrime_px", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_py", "tPrime_py", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_pz", "tPrime_pz", 6000, 0, 6000);

    histMan.addH1D_BJetBinned("ttbar_pt", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_withMETCut", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_withMETAndAsymJets", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_withAsymJetsCut", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withMETCut", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withMETCut", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withMETCut", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withMETAndAsymJets", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withMETAndAsymJets", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withMETAndAsymJets", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withAsymJetsCut", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withAsymJetsCut", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution", "ttbar_pt_vs_mttbar_2ndSolution", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution", "ttbar_pt_vs_mttbar_3rdSolution", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions", "ttbar_pt_vs_mttbar_allSolutions", 500, 0, 5000, 500,
            0, 5000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500, 0,
            5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0,
            5000, 100, 0, 1000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
//
    histMan.addH1D_BJetBinned("ttbar_px", "ttbar_px", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_py", "ttbar_py", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pz", "ttbar_pz", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_QCDEnriched", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("HT", "HT", 5000, 0, 5000);
    histMan.addH2D_BJetBinned("HTvsMttbar", "HT vs mttbar", 500, 0, 5000, 500, 0, 5000);
    histMan.addH1D("numberOfJets", "numberOfJets", 10, 0, 10);
    histMan.addH1D("numberOfBJets", "numberOfBJets", 10, 0, 10);
    histMan.addH1D_BJetBinned("MET", "MET", 200, 0, 1000);
    histMan.addH2D_BJetBinned("METvsMttbar", "MET vs mttbar", 500, 0, 5000, 200, 0, 1000);
    histMan.addH1D_BJetBinned("leadingJetMass", "leadingJetMass", 200, 0, 200);
    histMan.addH1D_BJetBinned("leadingGenJetMass", "leadingGenJetMass", 200, 0, 200);
    histMan.addH1D_BJetBinned("mtW", "mtW", 600, 0, 600);
    histMan.addH1D("electronD0", "electronD0", 1000, 0, 0.2);
    histMan.addH1D_BJetBinned("neutrino_pz", "neutrino_pz", 1000, -500, 500);
    histMan.addH2D("ptRel_vs_DRmin", "ptRel_vs_DRmin", 100, 0, 1, 300, 0, 300);
    histMan.addH1D("ptRel_QCDenriched", "ptRel_QCDenriched", 300, 0, 300);
    histMan.addH1D("DRmin_QCDenriched", "DRmin_QCDenriched", 100, 0, 1);
    histMan.addH1D("ptRel_WZenriched", "ptRel_WZenriched", 300, 0, 300);
    histMan.addH1D("DRmin_WZenriched", "DRmin_WZenriched", 100, 0, 1);
    histMan.addH1D_JetBinned("QCDest_CombRelIso", "RelIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion", "RelIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_1btag", "RelIso (>=1 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion_1btag", "RelIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_2btag", "RelIso (>=2 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion_2btag", "RelIso control region", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithMETCut", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithMETCutAndAsymJetCuts", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithAsymJetCuts", "PFIso", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_1btag", "PFIso (>=1 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_1btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_2btag", "PFIso (>=2 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_2btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithMETCut", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithMETCutAndAsymJetCuts", "PFIso control region",
            1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithAsymJetCuts", "PFIso control region", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCut", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCutAndAsymJetCuts", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithAsymJetCuts", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_1btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_2btag", "PFIso control region", 1000, 0, 10);

    histMan.addH1D_BJetBinned("pt_leadingTop", "pt_leadingTop", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop", "pt_NextToLeadingTop", 1000, 0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar", "pt_leadingTop_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar", "pt_NextToLeadingTop_vs_mttbar", 500, 0, 5000, 100, 0,
            1000);
    histMan.addH1D_BJetBinned("pt_leadingTop_withMETCut", "pt_leadingTop_withMETCut", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withMETCut", "pt_NextToLeadingTop_withMETCut", 1000, 0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETCut", "pt_leadingTop_vs_mttbar_withMETCut", 500, 0,
            5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETCut", "pt_NextToLeadingTop_vs_mttbar_withMETCut",
            500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("pt_leadingTop_withMETAndAsymJets", "pt_leadingTop_withMETAndAsymJets", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withMETAndAsymJets", "pt_NextToLeadingTop_withMETAndAsymJets", 1000,
            0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETAndAsymJets",
            "pt_leadingTop_vs_mttbar_withMETAndAsymJets", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets",
            "pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets", 500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("pt_leadingTop_withAsymJetsCut", "pt_leadingTop_withMETAndAsymJets", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withAsymJetsCut", "pt_NextToLeadingTop_withMETAndAsymJets", 1000, 0,
            1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withAsymJetsCut", "pt_leadingTop_vs_mttbar_withAsymJetsCut",
            500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut",
            "pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut", 500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("angleTops", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withMETCut", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withMETAndAsymJets", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withAsymJetsCut", "angle between top quarks", 400, 0, 4);

    histMan.addH2D_BJetBinned("angleTops_vs_mttbar", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withMETCut", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withMETAndAsymJets", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withAsymJetsCut", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);

}

// **************************************************************
// ***** End of histogram creation function
// **************************************************************


// **************************************************************
// ***** Selection functions used for special-purpose histograms
// **************************************************************

bool TopPairCandidateTutorial::passesRelIsoSelection() const{
    bool passesFirst3 = passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex);
        bool passGoodElectrons = goodElectrons.size() > 0 && goodIsolatedElectrons.size() < 2;
        bool passesBothIsolationvetos = false;
        if (passGoodElectrons) {
            const ElectronPointer electron = MostIsolatedElectron();
            if (electron->isGood()) {
                if (useCustomConversionTagger) {
                    conversionTagger->calculateConversionVariables(electron, tracks, 3.8, 0.45);
                    passesBothIsolationvetos = electron->isFromConversion() == false && conversionTagger->isFromConversion(
                            0.02, 0.02) == false;
                }
                else{
                    passesBothIsolationvetos = electron->isFromConversion() == false && electron->isTaggedAsConversion(
                            0.02, 0.02) == false;
                }
            }

        }
        bool muonVeto = hasNoIsolatedMuon();
        bool Zveto = isNotAZBosonEvent();
        return passesFirst3 && passGoodElectrons && passesBothIsolationvetos && muonVeto && Zveto;
}

bool TopPairCandidateTutorial::passesRelIsoControlSelection() const{
    bool passesFirst3 = passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex);
       bool passGoodElectrons = allElectrons.size() > 0 && goodIsolatedElectrons.size() < 2;
       bool passesBothIsolationvetos = false;
       if (passGoodElectrons) {
           const ElectronPointer electron = MostIsolatedElectron();
           if (electron->isQCDElectron()) {
               if (useCustomConversionTagger) {
                   conversionTagger->calculateConversionVariables(electron, tracks, 3.8, 0.45);
                   passesBothIsolationvetos = electron->isFromConversion() == false && conversionTagger->isFromConversion(
                           0.02, 0.02) == false;
               } else {
                   passesBothIsolationvetos = electron->isFromConversion() == false && electron->isTaggedAsConversion(
                           0.02, 0.02) == false;
               }
           }

       }
       bool muonVeto = hasNoIsolatedMuon();
       bool Zveto = isNotAZBosonEvent();
       return passesFirst3 && passGoodElectrons && passesBothIsolationvetos && muonVeto && Zveto;
}

bool TopPairCandidateTutorial::passesPFIsoSelection() const{
    bool passesFirst3 = passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex);
        bool passGoodElectrons = goodElectrons.size() > 0 && goodPFIsolatedElectrons.size() < 2;
        bool passesBothIsolationvetos = false;
        if (passGoodElectrons) {
            const ElectronPointer electron = MostPFIsolatedElectron();
            if (electron->isGood()) {
                if (useCustomConversionTagger) {
                    conversionTagger->calculateConversionVariables(electron, tracks, 3.8, 0.45);
                    passesBothIsolationvetos = electron->isFromConversion() == false && conversionTagger->isFromConversion(
                            0.02, 0.02) == false;
                }
                else{
                    passesBothIsolationvetos = electron->isFromConversion() == false && electron->isTaggedAsConversion(
                            0.02, 0.02) == false;
                }
            }

        }
        bool muonVeto = hasNoIsolatedMuon();
        bool Zveto = isNotAZBosonEvent();
        return passesFirst3 && passGoodElectrons && passesBothIsolationvetos && muonVeto && Zveto;
}

bool TopPairCandidateTutorial::passesPFIsoControlSelection() const{
    bool passesFirst3 = passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex);
       bool passGoodElectrons = allElectrons.size() > 0 && goodPFIsolatedElectrons.size() < 2;
       bool passesBothIsolationvetos = false;
       if (passGoodElectrons) {
           const ElectronPointer electron = MostPFIsolatedElectron();
           if (electron->isQCDElectron()) {
               if (useCustomConversionTagger) {
                   conversionTagger->calculateConversionVariables(electron, tracks, 3.8, 0.45);
                   passesBothIsolationvetos = electron->isFromConversion() == false && conversionTagger->isFromConversion(
                           0.02, 0.02) == false;
               } else {
                   passesBothIsolationvetos = electron->isFromConversion() == false && electron->isTaggedAsConversion(
                           0.02, 0.02) == false;
               }
           }

       }
       bool muonVeto = hasNoIsolatedMuon();
       bool Zveto = isNotAZBosonEvent();
       return passesFirst3 && passGoodElectrons && passesBothIsolationvetos && muonVeto && Zveto;
}

bool TopPairCandidateTutorial::passesConversionSelection() const {
    bool passesFirst6 = passesSelectionStepUpTo(TTbarEPlusJetsSelection::Zveto);
    bool isConversion1 = isolatedElectronDoesNotComeFromConversion() == false;
    bool isConversion2 = isolatedElectronNotTaggedAsFromConversion() == false;
    bool atLeast4Jets = hasAtLeastFourGoodJets();
    return passesFirst6 && (isConversion1 || isConversion2) && atLeast4Jets;
}

bool TopPairCandidateTutorial::passesAntiIsolationSelection() const {
    //require at least one good electron and no isolated good electrons
    if (!(goodElectrons.size() > 0 && goodIsolatedElectrons.size() == 0))
            return false;

    bool passesFirst3 = passesSelectionStep(TTbarEPlusJetsSelection::GoodPrimaryvertex);


    bool muonVeto = passesSelectionStep(TTbarEPlusJetsSelection::LooseMuonVeto);
    bool zveto = passesSelectionStep(TTbarEPlusJetsSelection::Zveto);
    bool conversionVeto = (goodElectrons.front()->isFromConversion() || goodElectrons.front()->isTaggedAsConversion(
            0.2, 0.2)) == false;
    bool jets = hasAtLeastFourGoodJets();
    return passesFirst3 && muonVeto && zveto && conversionVeto && jets;
}

// **************************************************************
// ***** End of section of special selection functions
// **************************************************************


// **************************************************************
// ***** Additional analyses to the main analysis
// ***** These functions provide examples of how to add new types
// ***** of analyses to the main analysis
// **************************************************************

void IntroAnalysis::doDiElectronAnalysis() {
    ElectronCollection electrons = currentEvent.GoodElectrons();
    if (electrons.size() == 2) {
        ElectronPointer leadingElectron = electrons.front();
        ElectronPointer secondElectron = electrons.at(1);
        histMan.H1D_JetBinned("diElectronMass")->Fill(leadingElectron->invariantMass(secondElectron), weight);
    }
}


void IntroAnalysis::doJetAnalysis() {
    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex)) {
        for(unsigned short jetIndex = 0; jetIndex < ttbarCandidate.Jets().size(); ++jetIndex)
            histMan.H1D_BJetBinned("AllJetMass")->Fill(ttbarCandidate.Jets().at(jetIndex)->mass());

        for (unsigned short jetIndex = 0; jetIndex < ttbarCandidate.GoodJets().size(); ++jetIndex) {
            double jetMass = ttbarCandidate.Jets().at(jetIndex)->mass();
            histMan.H1D_BJetBinned("AllGoodJetMass")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastOneGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastOneJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastTwoGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastTwoJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastThreeGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastThreeJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastFourGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastFourJets")->Fill(jetMass, weight);
        }
    }
}


void IntroAnalysis::doNotePlots() {
    if (ttbarCandidate.GoodElectrons().size() >= 1 && ttbarCandidate.Jets().size() >= 2) {
        const ElectronCollection electrons = ttbarCandidate.GoodElectrons();
        ElectronCollection nonConversionElectrons;
        for (unsigned int index = 0; index < electrons.size(); ++index) {
            const ElectronPointer electron = electrons.at(index);
            if (electron->isFromConversion() == false && electron->isTaggedAsConversion(0.02,0.02) == false) {
//                ConversionTagger tagger = ConversionTagger();
//                tagger.calculateConversionVariables(electron, ttbarCandidate.Tracks(), 3.8, 0.45);
//                if (tagger.isFromConversion(0.02, 0.02) == false)
                    nonConversionElectrons.push_back(electron);
            }
        }
        if (nonConversionElectrons.size() == 1) {
            const ElectronPointer electron = nonConversionElectrons.front();
            JetCollection goodjets;
            for (unsigned index = 0; index < ttbarCandidate.Jets().size(); ++index) {
                if (ttbarCandidate.Jets().at(index)->isGood())
                    goodjets.push_back(ttbarCandidate.Jets().at(index));
            }
            if (goodjets.size() >= 2) {
                unsigned int closestID = electron->getClosestJetIndex(goodjets);
                float minDR = electron->deltaR(goodjets.at(closestID));
                float ptRel = electron->relativePtTo(goodjets.at(closestID));
                histMan.H2D("ptRel_vs_DRmin")->Fill(minDR, ptRel, weight);
                if (ttbarCandidate.MET()->et() < 20 && ttbarCandidate.transverseWmass(electron) < 35) {
                    histMan.H1D("DRmin_QCDenriched")->Fill(minDR, weight);
                    histMan.H1D("ptRel_QCDenriched")->Fill(ptRel, weight);
                } else if (ttbarCandidate.MET()->et() > 30 && ttbarCandidate.transverseWmass(electron) > 50) {
                    histMan.H1D("DRmin_WZenriched")->Fill(minDR, weight);
                    histMan.H1D("ptRel_WZenriched")->Fill(ptRel, weight);
                }
            }
        }

    }
}


void IntroAnalysis::doQCDStudy() {
    if (ttbarCandidate.passesRelIsoSelection()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_CombRelIso")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_CombRelIso_1btag")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_CombRelIso_2btag")->Fill(electron->relativeIsolation(), weight);
    }

    if (ttbarCandidate.passesPFIsoSelection() && NTupleEventReader::electronAlgorithm
            == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_PFIsolation")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_PFIsolation_1btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_PFIsolation_2btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.MET()->pt() > 20) {
            histMan.H1D_JetBinned("QCDest_PFIsolation_WithMETCut")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodJets().size() >= 2) {
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_JetBinned("QCDest_PFIsolation_WithMETCutAndAsymJetCuts")->Fill(electron->pfIsolation(),
                            weight);
                }
            }
        }

        if (ttbarCandidate.GoodJets().size() >= 2) {
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_JetBinned("QCDest_PFIsolation_WithAsymJetCuts")->Fill(electron->pfIsolation(), weight);

            }
        }
    }


    if (ttbarCandidate.passesRelIsoControlSelection()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion_1btag")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion_2btag")->Fill(electron->relativeIsolation(), weight);

        if (NTupleEventReader::electronAlgorithm == ElectronAlgorithm::ParticleFlow) {
            const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodBJets().size() >= 1)
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_1btag")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodBJets().size() >= 2)
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_2btag")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCut")->Fill(electron->pfIsolation(), weight);

                if (ttbarCandidate.GoodJets().size() >= 2) {
                    if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                        histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCutAndAsymJetCuts")->Fill(
                                electron->pfIsolation(), weight);
                    }
                }
            }

            if (ttbarCandidate.GoodJets().size() >= 2) {
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithAsymJetCuts")->Fill(electron->pfIsolation(), weight);

                }
            }
        }

    }

    if (ttbarCandidate.passesPFIsoControlSelection() && NTupleEventReader::electronAlgorithm
            == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion_1btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion_2btag")->Fill(electron->pfIsolation(), weight);
    }

    if (ttbarCandidate.passesRelIsoSelection() && ttbarCandidate.hasAtLeastFourGoodJets()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron(false);
        if (electron->isIsolated() == false && !isnan(electron->relativeIsolation()) && !isinf(
                electron->relativeIsolation())) {
            try {
                ttbarCandidate.reconstructTTbar(electron);
                const ParticlePointer resonance = ttbarCandidate.getResonance();
                histMan.H1D_BJetBinned("mttbar_controlRegion")->Fill(resonance->mass(), weight);
                if (ttbarCandidate.MET()->pt() > 20) {
                    histMan.H1D_BJetBinned("mttbar_controlRegion_withMETCut")->Fill(resonance->mass(), weight);
                    if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                        histMan.H1D_BJetBinned("mttbar_controlRegion_withMETAndAsymJets")->Fill(resonance->mass(),
                                weight);
                }

                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                    histMan.H1D_BJetBinned("mttbar_controlRegion_withAsymJetsCut")->Fill(resonance->mass(), weight);

            } catch (ReconstructionException& e) {
                cout << "Could not reconstruct event: " << e.what() << endl;
            }
        }
    }

    if (ttbarCandidate.passesConversionSelection()) {
        ElectronPointer electron;
        if (Event::usePFIsolation)
            electron = currentEvent.GoodPFIsolatedElectrons().front();
        else
            electron = currentEvent.GoodIsolatedElectrons().front();
        try {
            ttbarCandidate.reconstructTTbar(electron);
            const ParticlePointer resonance = ttbarCandidate.getResonance();
            histMan.H1D_BJetBinned("mttbar_conversions")->Fill(resonance->mass(), weight);
            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_conversions_withMETCut")->Fill(resonance->mass(), weight);
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                    histMan.H1D_BJetBinned("mttbar_conversions_withMETAndAsymJets")->Fill(resonance->mass(), weight);
            }
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                histMan.H1D_BJetBinned("mttbar_conversions_withAsymJetsCut")->Fill(resonance->mass(), weight);
        } catch (ReconstructionException& e) {
            cout << "Could not reconstruct event: " << e.what() << endl;
        }

    }

    if (ttbarCandidate.Electrons().size() > 0 && ttbarCandidate.GoodPFIsolatedElectrons().size() < 2
            && NTupleEventReader::electronAlgorithm == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron(true);
        histMan.H1D_JetBinned("MostPFIsolatedElectron_dPhiIn")->Fill(electron->dPhiIn(), weight);
        histMan.H1D_JetBinned("MostPFIsolatedElectron_dEtaIn")->Fill(electron->dEtaIn(), weight);
    }

    if (ttbarCandidate.passesAntiIsolationSelection()) {
        ElectronPointer electron = ttbarCandidate.GoodElectrons().front();
        try {
            ttbarCandidate.reconstructTTbar(electron);
            float mttbar = ttbarCandidate.mttbar();
            histMan.H1D_BJetBinned("mttbar_antiIsolated")->Fill(mttbar, weight);
            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_antiIsolated_withMETCut")->Fill(mttbar, weight);
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_BJetBinned("mttbar_antiIsolated_withMETAndAsymJets")->Fill(mttbar, weight);
                }
            }
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("mttbar_antiIsolated_withAsymJetsCut")->Fill(mttbar, weight);
            }

        } catch (ReconstructionException& e) {
            cout << "Could not reconstruct event: " << e.what() << endl;
        }
    }

}


void IntroAnalysis::doMCttbarReconstruction() {
	MCParticlePointer top, antitop, b_from_top, b_from_antitop, W_plus, W_minus, electron, neutrino, quark_from_W, antiquark_from_W;
	JetCollection genJets = currentEvent.GenJets();
	JetPointer topBjet, antitopBjet, jet1fromW, jet2fromW;
	TtbarHypothesis MCttbarEvent;
	bool ejets_event = false;
	bool leptonic_Wplus_found = false, leptonic_Wminus_found = false;
	bool hadronic_Wplus_found = false, hadronic_Wminus_found = false;
	bool fully_hadronic_event = false, fully_leptonic_event = false;
	bool non_electron_leptonic_channel = false;
	int index = 0;
	int top_index = -100, antitop_index = -100, W_plus_index = -100, W_minus_index = -100, electron_index = -100, neutrino_index = -100,
			b_from_top_index = -100, b_from_antitop_index = -100, quark_from_W_index = -100, antiquark_from_W_index = -100;

	// MC ttbar reconstruction
	for (MCParticleCollection::const_iterator mc_particle = currentEvent.GenParticles().begin(); mc_particle != currentEvent.GenParticles().end(); ++mc_particle, ++index) {

		if ((*mc_particle)->status() != 3) continue;
		//top quark
		if ((*mc_particle)->pdgId() == 6) {
			top = *mc_particle;
			top_index = index;
			continue;
		}

		//anti-top quark
		if ((*mc_particle)->pdgId() == -6) {
			antitop = *mc_particle;
			antitop_index = index;
			continue;
		}

		//W� bosons
		if (((*mc_particle)->pdgId() == 24) && ((*mc_particle)->motherIndex() == top_index)) {
			W_plus = *mc_particle;
			W_plus_index = index;
			continue;
		}

		if (((*mc_particle)->pdgId() == -24) && ((*mc_particle)->motherIndex() == antitop_index)) {
			W_minus = *mc_particle;
			W_minus_index = index;
			continue;
		}

		//b-quarks
		if (((*mc_particle)->pdgId() == 5) && ((*mc_particle)->motherIndex() == top_index)) {
			b_from_top = *mc_particle;
			b_from_top_index = index;
			continue;
		}
		if (((*mc_particle)->pdgId() == -5) && ((*mc_particle)->motherIndex() == antitop_index)) {
			b_from_antitop = *mc_particle;
			b_from_antitop_index = index;
			continue;
		}

		//W+ decay products
		if ((*mc_particle)->motherIndex()==W_plus_index) {
			if ((*mc_particle)->pdgId() == -11) {
				electron = *mc_particle;
				electron_index = index;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->pdgId() == 12) {
				neutrino = *mc_particle;
				neutrino_index = index;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->isLepton()) {
				non_electron_leptonic_channel = true;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->isQuark()  && ((*mc_particle)->pdgId()>0)) {
				quark_from_W = *mc_particle;
				quark_from_W_index = index;
				hadronic_Wplus_found = true;
			}

			else if ((*mc_particle)->isQuark() && ((*mc_particle)->pdgId()<0)) {
				antiquark_from_W = *mc_particle;
				antiquark_from_W_index = index;
				hadronic_Wplus_found = true;
			}

			else {
				cout << "Something went wrong: W+ has unusual decay products." << endl;
			}

		}

		//W- decay products
		if ((*mc_particle)->motherIndex()==W_minus_index) {
			if ((*mc_particle)->pdgId() == 11) {
				electron = *mc_particle;
				electron_index = index;
				leptonic_Wminus_found = true;
			}

			else if ((*mc_particle)->pdgId() == -12) {
				neutrino = *mc_particle;
				neutrino_index = index;
				leptonic_Wminus_found = true;
			}

			else if ((*mc_particle)->isLepton()) {
				leptonic_Wminus_found = true;
				non_electron_leptonic_channel = true;
			}

			else if ((*mc_particle)->isQuark()  && ((*mc_particle)->pdgId()>0)) {
				quark_from_W = *mc_particle;
				quark_from_W_index = index;
				hadronic_Wminus_found = true;
			}

			else if ((*mc_particle)->isQuark() && ((*mc_particle)->pdgId()<0)) {
				antiquark_from_W = *mc_particle;
				antiquark_from_W_index = index;
				hadronic_Wminus_found = true;
			}

			else {
				cout << "Something went wrong: W- has unusual decay products." << endl;
			}
		}
	}

	//classify the event
	if (((leptonic_Wplus_found) || (leptonic_Wminus_found)) && ((hadronic_Wplus_found) || (hadronic_Wminus_found))
			&& (!non_electron_leptonic_channel)) { ejets_event = true; }
	if (((leptonic_Wplus_found) || (leptonic_Wminus_found)) && (!hadronic_Wplus_found)
			&& (!hadronic_Wminus_found)) { fully_leptonic_event = true; }
	if (((hadronic_Wplus_found) || (hadronic_Wminus_found)) && (!leptonic_Wplus_found)
			&& (!leptonic_Wminus_found)) { fully_hadronic_event = true; }

	if (ejets_event) {
		//matching genJets and partons

		if (genJets.size()>0) {
			int closestJetQuarkFromWIndex = quark_from_W->getClosestJetIndex(genJets);
			float minDR_quarkW = quark_from_W->deltaR(genJets.at(closestJetQuarkFromWIndex));
			jet1fromW = genJets.at(closestJetQuarkFromWIndex);

			int closestJetAntiQuarkFromWIndex = antiquark_from_W->getClosestJetIndex(genJets);
			float minDR_antiquarkW = antiquark_from_W->deltaR(genJets.at(closestJetAntiQuarkFromWIndex));
			jet2fromW = genJets.at(closestJetAntiQuarkFromWIndex);

			int closestJetBfromTopIndex = b_from_top->getClosestJetIndex(genJets);
			float minDR_BfromTop = b_from_top->deltaR(genJets.at(closestJetBfromTopIndex));
			topBjet = genJets.at(closestJetBfromTopIndex);

			int closestJetBfromAntiTopIndex = b_from_antitop->getClosestJetIndex(genJets);
			float minDR_BfromAntiTop = b_from_antitop->deltaR(genJets.at(closestJetBfromAntiTopIndex));
			antitopBjet = genJets.at(closestJetBfromAntiTopIndex);

			//delta R between genJets and partons histograms
			histMan.H1D("deltaRjet1")->Fill(minDR_quarkW);
			histMan.H1D("deltaRjet2")->Fill(minDR_antiquarkW);
			histMan.H1D("deltaRjet3")->Fill(minDR_BfromTop);
			histMan.H1D("deltaRjet4")->Fill(minDR_BfromAntiTop);

			histMan.H1D("deltaRjet_sum")->Fill(minDR_quarkW);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_antiquarkW);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_BfromTop);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_BfromAntiTop);
		}

		if (leptonic_Wplus_found) {
			MCttbarEvent.leptonicTop = (ParticlePointer) top;
			MCttbarEvent.hadronicTop = (ParticlePointer) antitop;
			MCttbarEvent.leptonicW = (ParticlePointer) W_plus;
			MCttbarEvent.hadronicW = (ParticlePointer) W_minus;
			MCttbarEvent.leptonicBjet = topBjet;
			MCttbarEvent.hadronicBJet = antitopBjet;
			MCttbarEvent.jet1FromW = jet1fromW;
			MCttbarEvent.jet2FromW = jet2fromW;
			MCttbarEvent.neutrinoFromW = (ParticlePointer) neutrino;
			ElectronPointer e(new Electron(electron->energy(), electron->px(), electron->py(), electron->pz()));
			MCttbarEvent.electronFromW = e;
		}
		else if (hadronic_Wplus_found) {
			MCttbarEvent.leptonicTop = (ParticlePointer) antitop;
			MCttbarEvent.hadronicTop = (ParticlePointer) top;
			MCttbarEvent.leptonicW = (ParticlePointer) W_minus;
			MCttbarEvent.hadronicW = (ParticlePointer) W_plus;
			MCttbarEvent.leptonicBjet = antitopBjet;
			MCttbarEvent.hadronicBJet = topBjet;
			MCttbarEvent.jet1FromW = jet1fromW;
			MCttbarEvent.jet2FromW = jet2fromW;
			MCttbarEvent.neutrinoFromW = (ParticlePointer) neutrino;
			ElectronPointer e(new Electron(electron->energy(), electron->px(), electron->py(), electron->pz()));
			MCttbarEvent.electronFromW = e;
		}
		else cout << "ERROR: no hadronic or leptonic W's in semileptonic event (nonsense).\n";

		//comparing deltaR between genJets from W and closest partons
		histMan.H2D("deltaR_genJets_partons")->Fill(MCttbarEvent.jet1FromW->deltaR(MCttbarEvent.jet2FromW),quark_from_W->deltaR(antiquark_from_W));

		//invariant mass histograms
		histMan.H1D("W_inv_mass_from_truth_partons")->Fill(quark_from_W->invariantMass(antiquark_from_W));
		histMan.H1D("W_inv_mass_from_genJets")->Fill(MCttbarEvent.jet1FromW->invariantMass(MCttbarEvent.jet2FromW));
		histMan.H1D("top_leptonic_inv_mass_from_truth")->Fill(MCttbarEvent.leptonicW->invariantMass(MCttbarEvent.leptonicBjet));
		histMan.H1D("top_hadronic_inv_mass_from_truth")->Fill(MCttbarEvent.hadronicW->invariantMass(MCttbarEvent.hadronicBJet));

		histMan.H1D("m3_mc")->Fill(MCttbarEvent.M3());

		// comparing truth and reco objects
		if (ttbarCandidate.passesFullTTbarEPlusJetSelection()) {
			histMan.H1D("m3_diff")->Fill(fabs(MCttbarEvent.M3()-ttbarCandidate.M3()));

			histMan.H1D("deltaRElectron")->Fill(MCttbarEvent.electronFromW->deltaR(ttbarCandidate.getElectronFromWDecay()));
			histMan.H1D("deltaRLeptonicBjet")->Fill(MCttbarEvent.leptonicBjet->deltaR(ttbarCandidate.getLeptonicBJet()));
			histMan.H1D("deltaRHadronicBjet")->Fill(MCttbarEvent.hadronicBJet->deltaR(ttbarCandidate.getHadronicBJet()));
			histMan.H1D("deltaRjet1fromW")->Fill(MCttbarEvent.jet1FromW->deltaR(ttbarCandidate.getJet1FromHadronicW()));
			histMan.H1D("deltaRjet2fromW")->Fill(MCttbarEvent.jet2FromW->deltaR(ttbarCandidate.getJet2FromHadronicW()));
		}
	}
}


double TopPairCandidateTutorial::M3() const {
    double m3(0), max_pt(0);
    if (goodJets.size() >= 3) {
        for (unsigned int index1 = 0; index1 < goodJets.size() - 2; ++index1) {
            for (unsigned int index2 = index1 + 1; index2 < goodJets.size() - 1; ++index2) {
                for (unsigned int index3 = index2 + 1; index3 < goodJets.size(); ++index3) {
                    FourVector m3Vector(goodJets.at(index1)->getFourVector() + goodJets.at(index2)->getFourVector()
                            + goodJets.at(index3)->getFourVector());
                    double currentPt = m3Vector.Pt();
                    if (currentPt > max_pt) {
                        max_pt = currentPt;
                        m3 = m3Vector.M();
                    }
                }
            }
        }
    }

    return m3;
}


double TopPairCandidateTutorial::fullHT() const {
    double ht(met->pt());

    for (unsigned int index = 0; index < goodIsolatedElectrons.size(); ++index) {
        ht += goodIsolatedElectrons.at(index)->pt();
    }

    for (unsigned int index = 0; index < goodIsolatedMuons.size(); ++index) {
        ht += goodIsolatedMuons.at(index)->pt();
    }

    for (unsigned int index = 0; index < goodJets.size(); ++index) {
        ht += goodJets.at(index)->pt();
    }
    return ht;
}


double TopPairCandidateTutorial::transverseWmass(const ElectronPointer electron) const {
    double energySquared = pow(electron->et() + met->et(), 2);
    double momentumSquared = pow(electron->px() + met->px(), 2) + pow(electron->py() + met->py(), 2);
    double tMassSquared = energySquared - momentumSquared;

    if (tMassSquared > 0)
        return sqrt(tMassSquared);
    else
        return -1;
}


// **************************************************************
// ***** End of additional analyses section
// **************************************************************


// *** There are some additional utility functions called from this file.
// *** They can be found in bin/IntroAnalysisUtils.cpp.  They are just
// *** bookkeeping functions with no physics content.

