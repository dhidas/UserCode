/*
 * TopPairCandidateTutorial.h
 *
 */

#ifndef TOPPAIRCANDIDATETUTORIAL_H_
#define TOPPAIRCANDIDATETUTORIAL_H_

#include <boost/shared_ptr.hpp>
#include "Event.h"
#include "RecoObjects/Particle.h"
#include "RecoObjects/Electron.h"
#include "RecoObjects/Jet.h"
#include "Taggers/ConversionTagger.h"
#include "TtbarHypothesis.h"
#include "TopPairEventCandidate.h"


namespace BAT {


class TopPairCandidateTutorial: public Event {

protected:
    static double const matched_angle;// = 0.945666;
    static double const matched_angle_sigma;// = 0.311091;
    static double const matched_leptonic_top_mass;// = 178.377;
    static double const matched_leptonic_top_mass_sigma;// = 31.050;
    static double const matched_hadronic_W_mass;// = 89.9153;
    static double const matched_hadronic_W_mass_sigma;// = 13.8711;
    static double const matched_hadronic_top_mass;// = 182.191;
    static double const matched_hadronic_top_mass_sigma;// = 22.1484;
    static double const matched_ptratio;// = 0.18552;
    static double const matched_ptratio_sigma;// = 0.401973;
    static double const matched_pt_ttbarSystem;// = 0.0760939;
    static double const matched_pt_ttbarSystem_sigma;// = 0.0700391;
    static double const matched_HTSystem;// = 1;
    static double const matched_HTSystem_sigma;// = 0.1;
    static double const W_mass;

    ElectronPointer electronFromW;
    JetPointer leptonicBJet, hadronicBJet, jet1FromW, jet2FromW;
    ParticlePointer neutrino1, neutrino2, leptonicW1, leptonicW2, hadronicW, leptonicTop1, leptonicTop2, hadronicTop,
            ttbarResonance;
    unsigned short selectedNeutrino, currentSelectedNeutrino, hadronicBIndex, leptonicBIndex, jet1FromWIndex,
            jet2FromWIndex;
    bool doneReconstruction;
    ConversionTaggerPointer conversionTagger;
    bool doneConversionTagging;
    std::vector<TtbarHypothesisPointer> solutions;
    compare_totalChi2 compareSolutions;
public:
    static NeutrinoSelectionCriterion::value usedNeutrinoSelection;
    static TTbarReconstructionCriterion::value usedTTbarReconstruction;
    TopPairCandidateTutorial();
    TopPairCandidateTutorial(const Event& event);
    virtual ~TopPairCandidateTutorial();

    virtual const JetPointer getLeptonicBJet() const;
    virtual const JetPointer getHadronicBJet() const;
    virtual const JetPointer getJet1FromHadronicW() const;
    virtual const JetPointer getJet2FromHadronicW() const;
    virtual const ElectronPointer getElectronFromWDecay() const;
    virtual const ParticlePointer getNeutrinoFromWDecay() const;
    virtual const ParticlePointer getLeptonicW() const;
    // virtual const ParticlePointer getHadronicW() const;
    virtual const ParticlePointer getLeptonicTop() const;
    virtual const ParticlePointer getHadronicTop() const;
    virtual const ParticlePointer getResonance() const;

    virtual bool passesSelectionStep(TTbarEPlusJetsSelection::Step step) const;
    virtual bool passesSelectionStepUpTo(TTbarEPlusJetsSelection::Step upToStep) const;
    virtual bool passesNMinus1(TTbarEPlusJetsSelection::Step omittedStep) const;
    virtual bool passesRelIsoSelection() const;
    virtual bool passesRelIsoControlSelection() const;
    virtual bool passesPFIsoSelection() const;
    virtual bool passesPFIsoControlSelection() const;
    virtual bool passesAntiIsolationSelection() const;

    virtual bool passesConversionSelection() const;

    virtual bool passesScrapingFilter() const;
    virtual bool passesHighLevelTrigger() const;
    virtual bool hasOneGoodPrimaryVertex() const;
    virtual bool hasOnlyOneGoodIsolatedElectron() const;
    virtual bool isolatedElectronDoesNotComeFromConversion() const;
    virtual bool isolatedElectronNotTaggedAsFromConversion() const;
    virtual bool hasNoIsolatedMuon() const;

    virtual bool hasAtLeastOneGoodJet() const;
    virtual bool hasAtLeastTwoGoodJets() const;
    virtual bool hasAtLeastThreeGoodJets() const;
    virtual bool hasAtLeastFourGoodJets() const;

    virtual bool isNotAZBosonEvent() const;
    virtual bool passesFullTTbarEPlusJetSelection() const;
    // virtual bool hasIsolatedElectronInBarrelRegion() const;

		virtual void reconstructTTbar(ElectronPointer electron);
    // virtual void reconstructUsingSubjets();
    // virtual void reconstructUsingMCTruth();

		virtual double calculateTopMassDifference(unsigned short neutrinoSolution) const;
    virtual double getLeptonicChi2(unsigned short neutrinoSolution) const;
    virtual double getLeptonicChi2(double topMass, double angle) const;
    virtual double getHadronicChi2() const;
    virtual double getGlobalChi2(unsigned short neutrinoSolution) const;
    virtual double getTotalChi2(unsigned short neutrinoSolution) const;
    virtual double getTotalChi2();

    virtual double M3() const;
    virtual double mttbar() const;
    virtual double sumPt() const;
    virtual double HT(unsigned short jetLimit) const;
    virtual double fullHT() const;
    virtual double transverseWmass(const ElectronPointer electron) const;
    virtual double HTSystem() const;
    virtual double PtRatio() const;
    virtual double PtTtbarSystem(unsigned short neutrinoSolution) const;
    virtual const boost::array<double, 2> computeNeutrinoPz();
    virtual void reconstructNeutrinos();
    // virtual bool hasNextJetCombination();
    // virtual void selectNextJetCombination();
    virtual void inspectReconstructedEvent() const;
    // virtual unsigned int NJet() const;
    virtual const std::vector<TtbarHypothesisPointer>& Solutions() const;

protected:
    void throwExpeptionIfNotReconstructed() const;
    virtual void selectNeutrinoSolution();
    virtual void fillHypotheses();
    virtual const TtbarHypothesisPointer fillHypothesis(unsigned short int neutrinoSolution);

};

}

#endif /* TOPPAIRCANDIDATETUTORIAL_H_ */
