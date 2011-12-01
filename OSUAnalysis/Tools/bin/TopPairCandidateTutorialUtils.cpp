// Utility members of the TopPairCandidateTutorial class.
// This file contains no physics details.

#include "../interface/TopPairCandidateTutorial.h"
#include <iostream>
#include <iomanip>


using namespace std;

namespace BAT {

NeutrinoSelectionCriterion::value TopPairCandidateTutorial::usedNeutrinoSelection; 
TTbarReconstructionCriterion::value TopPairCandidateTutorial::usedTTbarReconstruction;


TopPairCandidateTutorial::TopPairCandidateTutorial() :
    Event(),
    electronFromW(),
    leptonicBJet(),
    hadronicBJet(),
    jet1FromW(),
    jet2FromW(),
    neutrino1(),
    neutrino2(),
    leptonicW1(),
    leptonicW2(),
    hadronicW(),
    leptonicTop1(),
    leptonicTop2(),
    hadronicTop(),
    selectedNeutrino(0),
    currentSelectedNeutrino(0),
    hadronicBIndex(0),
    leptonicBIndex(0),
    jet1FromWIndex(0),
    jet2FromWIndex(0),
    doneReconstruction(false),
    conversionTagger(new ConversionTagger()),
    doneConversionTagging(false),
    solutions(),
    compareSolutions(){
}

TopPairCandidateTutorial::TopPairCandidateTutorial(const Event& event) :
    Event(event),
    electronFromW(),
    leptonicBJet(),
    hadronicBJet(),
    jet1FromW(),
    jet2FromW(),
    neutrino1(),
    neutrino2(),
    leptonicW1(),
    leptonicW2(),
    hadronicW(),
    leptonicTop1(),
    leptonicTop2(),
    hadronicTop(),
    selectedNeutrino(0),
    currentSelectedNeutrino(0),
    hadronicBIndex(0),
    leptonicBIndex(0),
    jet1FromWIndex(0),
    jet2FromWIndex(0),
    doneReconstruction(false),
    conversionTagger(new ConversionTagger()),
    doneConversionTagging(false) {

}

TopPairCandidateTutorial::~TopPairCandidateTutorial() {
}



bool TopPairCandidateTutorial::hasOneGoodPrimaryVertex() const {
    return PrimaryVertex()->isGood();
}

bool TopPairCandidateTutorial::hasOnlyOneGoodIsolatedElectron() const {
    if(Event::usePFIsolation)
        return goodPFIsolatedElectrons.size() == 1;
    else
        return goodIsolatedElectrons.size() == 1;
}

bool TopPairCandidateTutorial::isolatedElectronDoesNotComeFromConversion() const {
    bool passConversion = false;
    if (Event::usePFIsolation) {
        if (goodPFIsolatedElectrons.size() > 0)
            passConversion = goodPFIsolatedElectrons.front()->isFromConversion() == false;
    } else {
        if (goodIsolatedElectrons.size() > 0)
            passConversion = goodIsolatedElectrons.front()->isFromConversion() == false;
    }

    return passConversion;
}

bool TopPairCandidateTutorial::isolatedElectronNotTaggedAsFromConversion() const {
    bool passConversion = false;
    ElectronPointer electron;
    if (Event::usePFIsolation) {
        if (goodPFIsolatedElectrons.size() > 0)
            electron = goodPFIsolatedElectrons.front();
    } else {
        if (goodIsolatedElectrons.size() > 0)
            electron = goodIsolatedElectrons.front();
    }
    if (electron != 0) {
        if (useCustomConversionTagger) {
            conversionTagger->calculateConversionVariables(electron, tracks, 3.8, 0.45);
            passConversion = conversionTagger->isFromConversion(0.02, 0.02) == false;
        } else {
            passConversion = electron->isTaggedAsConversion(0.02, 0.02) == false;
        }
    }

    return passConversion;
}

bool TopPairCandidateTutorial::hasNoIsolatedMuon() const {
    return goodIsolatedMuons.size() == 0;
}

bool TopPairCandidateTutorial::hasAtLeastOneGoodJet() const {
    return goodJets.size() >= 1;
}

bool TopPairCandidateTutorial::hasAtLeastTwoGoodJets() const {
    return goodJets.size() >= 2;
}

bool TopPairCandidateTutorial::hasAtLeastThreeGoodJets() const {
    return goodJets.size() >= 3;
}

bool TopPairCandidateTutorial::hasAtLeastFourGoodJets() const {
    return goodJets.size() >= 4;
}


void TopPairCandidateTutorial::throwExpeptionIfNotReconstructed() const {
    if (doneReconstruction == false)
        throw ReconstructionException("Can't access reconstructed particles before reconstruction.");
}

const ElectronPointer TopPairCandidateTutorial::getElectronFromWDecay() const {
    return electronFromW;
}

const ParticlePointer TopPairCandidateTutorial::getNeutrinoFromWDecay() const {
    throwExpeptionIfNotReconstructed();
    if (selectedNeutrino == 1)
        return neutrino1;
    else
        return neutrino2;
}

const JetPointer TopPairCandidateTutorial::getHadronicBJet() const {
    throwExpeptionIfNotReconstructed();
    return hadronicBJet;
}

const JetPointer TopPairCandidateTutorial::getLeptonicBJet() const {
    throwExpeptionIfNotReconstructed();
    return leptonicBJet;
}


const ParticlePointer TopPairCandidateTutorial::getLeptonicW() const
{
    throwExpeptionIfNotReconstructed();
    if (selectedNeutrino == 1)
        return leptonicW1;
    else
        return leptonicW2;
}


const JetPointer TopPairCandidateTutorial::getJet1FromHadronicW() const {
    throwExpeptionIfNotReconstructed();
    return jet1FromW;
}

const JetPointer TopPairCandidateTutorial::getJet2FromHadronicW() const {
    throwExpeptionIfNotReconstructed();
    return jet2FromW;
}

const ParticlePointer TopPairCandidateTutorial::getLeptonicTop() const {
    throwExpeptionIfNotReconstructed();
    if (selectedNeutrino == 1)
        return leptonicTop1;
    else
        return leptonicTop2;
}

const ParticlePointer TopPairCandidateTutorial::getHadronicTop() const {
    throwExpeptionIfNotReconstructed();
    return hadronicTop;
}

const ParticlePointer TopPairCandidateTutorial::getResonance() const {
    throwExpeptionIfNotReconstructed();
    return ttbarResonance;
}



double TopPairCandidateTutorial::mttbar() const {
    return getResonance()->mass();
}

const std::vector<TtbarHypothesisPointer>& TopPairCandidateTutorial::Solutions() const{
	return solutions;
}

void TopPairCandidateTutorial::inspectReconstructedEvent() const {
    cout << "run " << runNumber << ", event " << eventNumber << endl;
    cout << "leptonic b jet, goodJet index " << leptonicBIndex << endl;
    EventContentPrinter::printJet(leptonicBJet);

    cout << "electron from W" << endl;
    EventContentPrinter::printElectron(electronFromW);

    cout << "MET" << endl;
    EventContentPrinter::printParticle(met);
    cout << endl;

    cout << "reconstructed neutrino 1(selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(neutrino1);
    cout << endl;

    cout << "reconstructed neutrino 2(selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(neutrino2);
    cout << endl;

    cout << "leptonic W 1 (selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(leptonicW1);
    cout << endl;

    cout << "leptonic W 2 (selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(leptonicW2);
    cout << endl;

    cout << "leptonic top 1 (selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(leptonicTop1);
    cout << endl;

    cout << "leptonic top 2 (selected: " << selectedNeutrino << ")" << endl;
    EventContentPrinter::printParticle(leptonicTop2);
    cout << endl;

    cout << "hadronic b jet, goodJet index " << hadronicBIndex << endl;
    EventContentPrinter::printJet(hadronicBJet);

    cout << "jet1 from W, goodJet index " << jet1FromWIndex << endl;
    EventContentPrinter::printJet(jet1FromW);

    cout << "jet 2 from W, goodJet index " << jet2FromWIndex << endl;
    EventContentPrinter::printJet(jet2FromW);

    cout << "hadronic W" << endl;
    EventContentPrinter::printParticle(hadronicW);
    cout << endl;

    cout << "hadronic top" << endl;
    EventContentPrinter::printParticle(hadronicTop);
    cout << endl;
}


}
