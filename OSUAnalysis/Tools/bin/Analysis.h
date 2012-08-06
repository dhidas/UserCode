/*
 * Analysis.h
 *
 *  Created on: 12 Jul 2010
 *      Author: kreczko
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_
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
#include "../interface/TopPairEventCandidate.h"
#include "../interface/ToplikeCandidate.h"
#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include "../interface/HistHelpers/HistogramManager.h"
#include "../interface/EventCounter.h"
#include "../interface/RecoObjects/Particle.h"
#include "../interface/EventWeightProvider.h"
#include "TFile.h"
#include "TTree.h"

struct InterestingEvent {
    InterestingEvent(unsigned long run, unsigned long event, std::string file) :
        candidate(), runNumber(run), eventNumber(event), fileName(file) {

    }

    InterestingEvent(BAT::TopPairEventCandidate cand, std::string file) :
        candidate(cand), runNumber(cand.runnumber()), eventNumber(cand.eventnumber()), fileName(file) {

    }
    ~InterestingEvent() {

    }
    BAT::TopPairEventCandidate candidate;
    unsigned long runNumber, eventNumber;
    std::string fileName;

    void print() {
        std::cout << "run " << candidate.runnumber() << ", event " << candidate.eventnumber() << " (Mttbar: "
                << candidate.mttbar() << ")" << std::endl;
        std::cout << "located in: " << fileName << std::endl << std::endl;
    }
};
typedef boost::array<unsigned long, BAT::TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS> cutarray;
typedef boost::unordered_map<std::string, cutarray> cutmap;

class Analysis {
private:
    boost::scoped_ptr<BAT::NTupleEventReader> eventReader;
    BAT::Event currentEvent;
    BAT::TopPairEventCandidate ttbarCandidate;
    BAT::TopPlusXCandidates tPlusXCandidates;
    BAT::TwoNonResTops twoTopsNonRes;
    BAT::TwoNonResTops goodTopBadTop;
    BAT::TopNoMassConstraint loneTopsNoMassConstr;
    BAT::LeptoTopNoMassConstraint leptoTopNoMassConstr;
public:
    BAT::HistogramManager histMan;
private:
    cutarray cutflow;
    cutarray singleCuts;
    cutmap cutflowPerFile;
    cutmap singleCutsPerFile;
    cutarray mucutflow;
    cutarray musingleCuts;
    cutmap mucutflowPerFile;
    cutmap musingleCutsPerFile;
    std::vector<InterestingEvent> interestingEvents, brokenEvents;
    std::map<unsigned long, std::vector<unsigned long> > eventCheck;
    BAT::EventWeightProvider weights;
    double weight;
    BAT::Counter cutflowPerSample;
    BAT::Counter mucutflowPerSample;

    unsigned int wpTopPos, wpTopNeg, specTopPos, specTopNeg;

    TFile *outputFile;
    TTree *microTuple;

    // general event information: 
    int    type, eventNumber, runNumber, numberOfJets, numberOfBJets, eSEL, muSEL, leptoCharge;
    double ST, MET;
    // leading non-b-jet:
    int    leadingJetPdgId, leadingJetIndGen;
    double leadingJetPtGen, leadingJetEtaGen, leadingJetPhiGen;
    double leadingJetPtRec, leadingJetEtaRec, leadingJetPhiRec;

    double JetPx[100];
    double JetPy[100];
    double JetPz[100];
    double JetE[100];
    int    JetBTag[100];

    // highest pT leftover jet (gen-level: radiation for the background and d-jet in signals)
    int    freeJetPdgId,    freeJetIndGen;
    double freeJetPtGen,    freeJetEtaGen,    freeJetPhiGen;
    double freeJetPtRec,    freeJetEtaRec,    freeJetPhiRec;
    // a candidate b-jet from the leptonic top:
    double bJetTlPtGen, bJetTlEtaGen, bJetTlPhiGen;
    double bJetTlPtRec, bJetTlEtaRec, bJetTlPhiRec;
    // a candidate b-jet from the hadronic top:
    double bJetThPtGen, bJetThEtaGen, bJetThPhiGen;
    double bJetThPtRec, bJetThEtaRec, bJetThPhiRec;
    // first and second (ordered in pT) jets from the hadronic W
    int    jet1WhPdgId, jet2WhPdgId,  jet1WhIndGen, jet2WhIndGen;
    double jet1WhPtGen, jet1WhEtaGen, jet1WhPhiGen;
    double jet1WhPtRec, jet1WhEtaRec, jet1WhPhiRec;
    double jet2WhPtGen, jet2WhEtaGen, jet2WhPhiGen;
    double jet2WhPtRec, jet2WhEtaRec, jet2WhPhiRec;
    // neutrino and missing ET
    double metPtGen,   metPhiGen,   metEtaGen;
    double metPtRec,   metPhiRec,   metEtaRec;
    // lepton
    int    leptonPdgId;
    double leptonPtGen, leptonEtaGen, leptonPhiGen;
    double leptonPtRec, leptonEtaRec, leptonPhiRec;
    // composite objects 
    //  leptonic and hadronic W and top
    int    tlPdgId,    thPdgId,     wlPdgId,   whPdgId;
    double wlPtGen,    wlEtaGen,    wlPhiGen,  wlMgen;
    double wlPtRec,    wlEtaRec,    wlPhiRec,  wlMrec;
    double tlPtGen,    tlEtaGen,    tlPhiGen,  tlMgen;
    double tlPtRec,    tlEtaRec,    tlPhiRec,  tlMrec;
    double whPtGen,    whEtaGen,    whPhiGen,  whMgen;
    double whPtRec,    whEtaRec,    whPhiRec,  whMrec;
    double thPtGen,    thEtaGen,    thPhiGen,  thMgen;
    double thPtRec,    thEtaRec,    thPhiRec,  thMrec;
    //  W' candidate
    int    wpPdgId;
    double wpPtGen,    wpEtaGen,    wpPhiGen;
    double tPosJetMassGen, tPosJetPtGen, tPosJetEtaGen, tPosJetPhiGen;
    double tNegJetMassGen, tNegJetPtGen, tNegJetEtaGen, tNegJetPhiGen;
    double tPosJetMassRec, tPosJetPtRec, tPosJetEtaRec, tPosJetPhiRec;
    double tNegJetMassRec, tNegJetPtRec, tNegJetEtaRec, tNegJetPhiRec;

    // rec-gen matching for the jets used in the top reco 
    double dRtlBjet, dRthBjet, dRwhJet1, dRwhJet2;
    // matching of gen jets to any jets around
    //double matchRtlJet, matchRthJet, matchRwhJet1, matchRwhJet2;

public:
    static float luminosity;
    Analysis();
    virtual ~Analysis();
    void analyze(int const Section = -1, TString const ListName = "");
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
        BAT::TopPairEventCandidate::useCustomConversionTagger = use;
        //custom conversion tagger needs track information
        BAT::NTupleEventReader::loadTracks = use;
    }
private:
    void printNumberOfProccessedEventsEvery(unsigned long printEvery);
//    void doEcalSpikeAnalysis();
    bool initiateEvent();
    void doDiElectronAnalysis();
    void doTTBarAnalysis();
		void fillMTplots(const double mt, const  float charge,
			const std::string histPrefix, const TString &leptSuffix,
			const unsigned int step);
		void doTransverseMassPlots();
    void doTTbarCutFlow();
		void doToplikeCutFlow(cutarray &singleCuts, cutmap &singleCutsPerFile,
		    cutarray &cutflow, cutmap &cutflowPerFile, BAT::Counter &cutflowPerSample,
					const BAT::ToplikeSelectionSteps::Step cutList[], unsigned int cutListSiz);
    void doSynchExercise();
    void printInterestingEvents();
    void printSummary();
    void inspectEvents();
    void createHistograms();
    void initMicroNtuple();
    void doMicroNtuple();
    void doNotePlots();
    void doQCDStudy();
    void doMCttbarReconstruction();
    void checkForDuplicatedEvents();
    void checkForBrokenEvents();
    void doJetAnalysis();
//    void doPileUpStudy();
};

#endif /* ANALYSIS_H_ */
