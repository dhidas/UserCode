/*
 * IntroAnalysis.h
 * Header file for declarations used by IntroAnalysis.cpp
 *
 */

#ifndef INTROANALYSIS_H_
#define INTROANALYSIS_H_
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include "../interface/Readers/NTupleEventReader.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "../interface/Event.h"
#include "../interface/TopPairCandidateTutorial.h"
// #include "../interface/ToplikeCandidate.h"
#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include "../interface/HistHelpers/HistogramManager.h"
#include "../interface/EventCounter.h"
#include "../interface/RecoObjects/Particle.h"
#include "Analysis.h"
#include "../interface/EventWeightProvider.h"

namespace IntroOAT {

struct InterestingEvent {
    InterestingEvent(unsigned long run, unsigned long event, std::string file) :
        candidate(), runNumber(run), eventNumber(event), fileName(file) {

    }

    InterestingEvent(BAT::TopPairCandidateTutorial cand, std::string file) :
        candidate(cand), runNumber(cand.runnumber()), eventNumber(cand.eventnumber()), fileName(file) {

    }
    ~InterestingEvent() {

    }
    BAT::TopPairCandidateTutorial candidate;
    unsigned long runNumber, eventNumber;
    std::string fileName;

    void print() {
        std::cout << "run " << candidate.runnumber() << ", event " << candidate.eventnumber() << " (Mttbar: "
                << candidate.mttbar() << ")" << std::endl;
        std::cout << "located in: " << fileName << std::endl << std::endl;
    }
};

} // end namespace IntroOAT

class IntroAnalysis {
private:
    boost::scoped_ptr<BAT::NTupleEventReader> eventReader;
    BAT::Event currentEvent;

    //************************************************
    // Reconstructed candidates must be declared here.
    //************************************************
    BAT::TopPairCandidateTutorial ttbarCandidate;
    BAT::ToplikeCandidate tPrimeCandidate;
    BAT::TopNoMassConstraint loneTopsNoMassConstr;

    BAT::HistogramManager histMan;
    cutarray cutflow;
    cutarray singleCuts;
    cutmap cutflowPerFile;
    cutmap singleCutsPerFile;
    std::vector<IntroOAT::InterestingEvent> interestingEvents, brokenEvents;
    std::map<unsigned long, std::vector<unsigned long> > eventCheck;
    BAT::EventWeightProvider weights;
    float weight;
    BAT::Counter cutflowPerSample;
public:
    static float luminosity;
    IntroAnalysis();
    virtual ~IntroAnalysis();
    void analyze();
    void addInputFile(const char * fileName);
    void setMaximalNumberOfEvents(long maxEvents);
    void setUsedNeutrinoSelectionForTopPairReconstruction(BAT::NeutrinoSelectionCriterion::value selection);
    void setUsedTTbarReconstructionCriterion(BAT::TTbarReconstructionCriterion::value selection);
    static void useJetAlgorithm(BAT::JetAlgorithm::value algo) {
        BAT::NTupleEventReader::jetAlgorithm = algo;
    }
    static void useElectronAlgorithm(BAT::ElectronAlgorithm::value algo) {
        BAT::NTupleEventReader::electronAlgorithm = algo;
    }
    static void useMETAlgorithm(BAT::METAlgorithm::value algo) {
        BAT::NTupleEventReader::metAlgorithm = algo;
    }
    static void useMuonAlgorithm(BAT::MuonAlgorithm::value algo){
        BAT::NTupleEventReader::muonAlgorithm = algo;
    }

    static void usePFIsolation(bool use){
        BAT::Event::usePFIsolation = use;
    }

    static void useCustomConversionTagger(bool use){
        BAT::TopPairCandidateTutorial::useCustomConversionTagger = use;
        //custom conversion tagger needs track information
        BAT::NTupleEventReader::loadTracks = use;
    }
private:
    void printNumberOfProccessedEventsEvery(unsigned long printEvery);
//    void doEcalSpikeAnalysis();
    void initiateEvent();
    void doDiElectronAnalysis();
    void doTTBarAnalysis();
    void doTTbarCutFlow();
    void doSynchExercise();
    void printInterestingEvents();
    void printSummary();
    void inspectEvents();
    void createHistograms();
    void doNotePlots();
    void doQCDStudy();
    void doMCttbarReconstruction();
    void checkForDuplicatedEvents();
    void checkForBrokenEvents();
    void doJetAnalysis();
};

#endif /* INTROANALYSIS_H_ */
