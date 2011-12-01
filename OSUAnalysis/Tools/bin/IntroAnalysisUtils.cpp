
#include <cmath>
#include <math.h>
// ***** IntroAnalysisUtils.cpp
// *****  Utility functions for supporting the main analysis code
// *****  There is no physics content in this file, only bookkeeping functions.

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>

#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"

#include "IntroAnalysis.h"
#include "../interface/EventCounter.h"
#include "../interface/Printers/EventTablePrinter.h"


using namespace ROOT;
using namespace std;
using namespace BAT;


float IntroAnalysis::luminosity = 36.145;  // Initialization
// Reset luminosity to correct value in "main".


void IntroAnalysis::printNumberOfProccessedEventsEvery(unsigned long printEvery) {
    unsigned long eventIndex = eventReader->getNumberOfProccessedEvents();
    if (eventIndex % printEvery == 0 || eventIndex == 1) {
        cout << "Analysing event no " << eventIndex << ", sample: " << DataType::names[currentEvent.getDataType()]
                << endl;
        cout << "File: " << eventReader->getCurrentFile() << endl;
    }

}


void IntroAnalysis::inspectEvents() {
    std::vector<IntroOAT::InterestingEvent> eventsToInspect;

    for (unsigned int index = 0; index < eventsToInspect.size(); ++index) {
        if ((ttbarCandidate.runnumber() == eventsToInspect.at(index).runNumber && ttbarCandidate.eventnumber()
                == eventsToInspect.at(index).eventNumber)) {
            cout << "file: " << eventReader->getCurrentFile() << endl;
            ttbarCandidate.inspect();
        }
    }

}

void IntroAnalysis::doSynchExercise() {
    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::ConversionFinder)) {
        cout << ttbarCandidate.runnumber() << ":" << ttbarCandidate.eventnumber() << ":" << endl;//electron->et() << endl;
        if (ttbarCandidate.eventnumber() == 450622) {
            ttbarCandidate.inspect();
        }
    }
}


// Count number of events that pass each cut
void IntroAnalysis::doTTbarCutFlow() {
    for (unsigned int cut = 0; cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
        if (ttbarCandidate.passesSelectionStep((TTbarEPlusJetsSelection::Step) cut)) {
            ++singleCuts[cut];
            singleCutsPerFile[eventReader->getCurrentFile()][cut]++;
        }

        if (ttbarCandidate.passesSelectionStepUpTo((TTbarEPlusJetsSelection::Step) cut)) {
            cutflow[cut] += 1;
            cutflowPerFile[eventReader->getCurrentFile()][cut]++;
            unsigned int njet = ttbarCandidate.GoodJets().size();
            if (njet >= JetBin::NUMBER_OF_JET_BINS)
                njet = JetBin::NUMBER_OF_JET_BINS - 1;
            cutflowPerSample.increase(ttbarCandidate.getDataType(), cut, njet, weight);
        }
    }
}


void IntroAnalysis::printInterestingEvents() {
    cout << "Interesting events:" << endl;
    for (unsigned int index = 0; index < interestingEvents.size(); ++index) {
        interestingEvents.at(index).print();
    }
}

void IntroAnalysis::printSummary() {
    EventTablePrinter::printCutFlowLatexTable(cutflowPerSample);
    EventTablePrinter::printUnweightedCutFlowLatexTable(cutflowPerSample);

    cout << "total number of processed events: " << eventReader->getNumberOfProccessedEvents() << endl;
    cout << endl;
    for (unsigned int cut = 0; cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
        cout << "Selection step '" << TTbarEPlusJetsSelection::StringSteps[cut] << "'" << endl;
        cout << "passed events (single cut): " << singleCuts.at(cut) << endl;
        if (cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS)
            cout << "passed events (up to this cut):" << cutflow.at(cut) << endl;
        else
            cout << "passed events (full selection):" << cutflow.at(cut) << endl;
        cout << endl;
    }

    cout << "number of events without electrons: " << brokenEvents.size() << endl;
}

IntroAnalysis::IntroAnalysis() :
    eventReader(new NTupleEventReader()),
//    eventFilter(Filter::makeTopPairEPlusJetsFilter()),
    currentEvent(),
    ttbarCandidate(),
    histMan(),
    cutflow(),
    singleCuts(),
    cutflowPerFile(),
    singleCutsPerFile(),
    interestingEvents(),
    brokenEvents(),
    eventCheck(),
    weights(IntroAnalysis::luminosity/*current lumi*/),
    weight(0),
    cutflowPerSample(DataType::NUMBER_OF_DATA_TYPES, TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS,
                    JetBin::NUMBER_OF_JET_BINS) {
    //    outputfile->SetCompressionLevel(7);
    for (unsigned int cut = 0; cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
        cutflow[cut] = 0;
        singleCuts[cut] = 0;
    }
}

IntroAnalysis::~IntroAnalysis() {
    histMan.writeToDisk();
}

void IntroAnalysis::addInputFile(const char* fileName) {
    eventReader->addInputFile(fileName);
}

void IntroAnalysis::setMaximalNumberOfEvents(long maxEvents) {
    if (maxEvents > 0) {
        eventReader->setMaximumNumberOfEvents(maxEvents);
    }
}

void IntroAnalysis::setUsedNeutrinoSelectionForTopPairReconstruction(NeutrinoSelectionCriterion::value selection) {
    TopPairCandidateTutorial::usedNeutrinoSelection = selection;
}

void IntroAnalysis::setUsedTTbarReconstructionCriterion(TTbarReconstructionCriterion::value selection) {
	TopPairCandidateTutorial::usedTTbarReconstruction = selection;
}

void IntroAnalysis::checkForDuplicatedEvents(){
    map<unsigned long, std::vector<unsigned long> >::const_iterator iter;
    std::vector<pair<unsigned long, unsigned long> > duplicateEvents;

    for(iter = eventCheck.begin(); iter != eventCheck.end(); ++iter){
        std::vector<unsigned long> events = (*iter).second;
        std::sort(events.begin(), events.end());
        for(unsigned long ev = 0; ev < events.size() -1; ++ev){
            if(events.at(ev) == events.at(ev +1)){
                duplicateEvents.push_back(make_pair((*iter).first, events.at(ev)));
            }
        }
    }

    if (duplicateEvents.size() > 0){
        cout << "found duplicate events" << endl;
        for(unsigned long ev = 0; ev < duplicateEvents.size() -1; ++ev){
            cout << "run: " << duplicateEvents.at(ev).first << " event: " << duplicateEvents.at(ev).second << endl;
        }
    }
}

void IntroAnalysis::checkForBrokenEvents(){
    if(ttbarCandidate.Electrons().size() == 0){
        brokenEvents.push_back(IntroOAT::InterestingEvent(ttbarCandidate, eventReader->getCurrentFile()));
    }

    if(ttbarCandidate.eventnumber() == 1019245){
        cout << "broken event" << endl;
        ttbarCandidate.inspect();
    }
}




