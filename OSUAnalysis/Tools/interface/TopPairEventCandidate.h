/*
 * TopPairEventCandidate.h
 *
 *  Created on: 9 Jul 2010
 *      Author: kreczko
 */

#ifndef TOPPAIREVENTCANDIDATE_H_
#define TOPPAIREVENTCANDIDATE_H_

#include <boost/shared_ptr.hpp>

#include "Event.h"
#include "ReconstructionRules/BasicRule.h"
#include "RecoObjects/Particle.h"
#include "RecoObjects/Lepton.h"
#include "RecoObjects/Jet.h"
#include "Taggers/ConversionTagger.h"
#include "TtbarHypothesis.h"

namespace BAT {

namespace NeutrinoSelectionCriterion {
enum value {
    pzClosestToLepton,
    mostCentral,
    pzClosestToLeptonOrMostcentralIfAbove300,
    largestValueOfCosine,
    chi2,
    TopMassDifference
};
}

namespace TTbarReconstructionCriterion {
enum value {
	TopMassDifference,
	chi2
};
}

namespace TTbarEPlusJetsSelection {
enum Step {
    FilterOutScraping,
    HighLevelTrigger,
    GoodPrimaryvertex,
    OneIsolatedElectron,
    LooseMuonVeto,
    Zveto,
    ConversionRejection,
    ConversionFinder,
    AtLeastOneGoodJets,
    AtLeastTwoGoodJets,
    AtLeastThreeGoodJets,
    AtLeastFourGoodJets,
    MissingTransverseEnergy,
    AsymmetricJetCuts,
    AtLeastOneBtag,
    AtLeastTwoBtags,
    NUMBER_OF_SELECTION_STEPS
};

const std::string StringSteps[NUMBER_OF_SELECTION_STEPS] = {
        "Scraping Filter",
        "High Level Trigger",
        "good primary vertex",
        "exactly one isolated electron",
        "loose muon veto",
        "Z veto",
        "Conversion veto (missing hits)",
        "Conversion finder (partner track)",
        ">= 1 jets",
        ">= 2 jets",
        ">= 3 jets",
		">= 4 jets",
        "MET > 20GeV",
        "jet pT: 1st > 70 GeV, 2nd > 50 GeV",
        ">=1 SSV b-tag",
        ">=2 SSV b-tag"

};
}

struct ReconstructionException: public std::exception {
    TString msg;
    ReconstructionException(TString message) :
        msg(message) {
    }
    ~ReconstructionException() throw () {
    }

    const char* what() const throw () {
        return msg;
    }
};

class TopPairEventCandidate: public Event {

protected:
		static const double pdgWmass;
		static const double wMassWidth;
		static const double pdgTopMass;
		static const double topMassWidth;

		// The following are obsolete and are only retained in
		// case the old ttbar code might be useful during
		// testing.
    static const double W_mass;
    static const double matched_leptonic_top_mass;
    static const double matched_leptonic_top_mass_sigma;
    static const double matched_hadronic_W_mass;
    static const double matched_hadronic_W_mass_sigma;
    static const double matched_hadronic_top_mass;
    static const double matched_hadronic_top_mass_sigma;
    static const double matched_angle;
    static const double matched_angle_sigma;
    static const double matched_ptratio;
    static const double matched_ptratio_sigma;
    static const double matched_pt_ttbarSystem;
    static const double matched_pt_ttbarSystem_sigma;
    static const double matched_HTSystem;
    static const double matched_HTSystem_sigma;

    // ElectronPointer electronFromW;
    LeptonPointer leptonFromW;
    JetPointer leptonicBJet, hadronicBJet, jet1FromW, jet2FromW;
    ParticlePointer neutrino1, neutrino2, leptonicW1, leptonicW2, hadronicW, leptonicTop1, leptonicTop2, hadronicTop,
            ttbarResonance;
public:
    unsigned short selectedNeutrino, currentSelectedNeutrino, hadronicBIndex, leptonicBIndex, jet1FromWIndex,
            jet2FromWIndex;
    virtual bool isReconstructed(void) const { return doneReconstruction; }

protected:
    bool doneReconstruction;
    ConversionTaggerPointer conversionTagger;
    bool doneConversionTagging;
    std::vector<TtbarHypothesisPointer> solutions;
    compare_totalChi2 compareSolutions;
public:
    static NeutrinoSelectionCriterion::value usedNeutrinoSelection;
    static TTbarReconstructionCriterion::value usedTTbarReconstruction;
    TopPairEventCandidate();
    TopPairEventCandidate(const Event& event);
    virtual ~TopPairEventCandidate();

    virtual const JetPointer getLeptonicBJet() const;
    virtual const JetPointer getHadronicBJet() const;
    virtual const JetPointer getJet1FromHadronicW() const;
    virtual const JetPointer getJet2FromHadronicW() const;
    virtual const LeptonPointer getLeptonFromWDecay() const;
    virtual const ParticlePointer getNeutrinoFromWDecay() const;
    virtual const ParticlePointer getLeptonicW() const;
    virtual const ParticlePointer getHadronicW() const;
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

    virtual bool passesMETCut() const;
    virtual bool passesAsymmetricJetCuts() const;
    virtual bool hasAtLeastOneBtag() const;
    virtual bool hasAtLeastTwoBtags() const;

		virtual void reconstructTTbar(LeptonPointer lepton);
    // virtual void reconstructUsingSubjets();
    // virtual void reconstructUsingMCTruth();
		void reconstruct(const Rule rule, const LeptonPointer lepton);

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
    virtual double transverseWmass(const LeptonPointer lepton) const;
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
    void throwExpeptionIfNotReconstructed(TString location = "") const;
    virtual void selectNeutrinoSolution();
    virtual void fillHypotheses();
    virtual const TtbarHypothesisPointer fillHypothesis(unsigned short int neutrinoSolution);

};

}

#endif /* TOPPAIREVENTCANDIDATE_H_ */
