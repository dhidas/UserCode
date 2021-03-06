/*
 * TTbarHypothesis.h
 *
 *  Created on: Dec 4, 2010
 *      Author: lkreczko
 */

#ifndef TTBARHYPOTHESIS_H_
#define TTBARHYPOTHESIS_H_
#include "RecoObjects/Electron.h"
#include "RecoObjects/Jet.h"
#include "RecoObjects/MET.h"
#include "RecoObjects/Particle.h"

namespace BAT {
class TtbarHypothesis {
public:
	TtbarHypothesis();
	virtual ~TtbarHypothesis();
	double totalChi2, leptonicChi2, hadronicChi2, globalChi2, disc;
	ParticlePointer hadronicTop, leptonicTop, leptonicW, hadronicW, resonance, neutrinoFromW;
	JetPointer leptonicBjet, hadronicBJet, jet1FromW, jet2FromW;
	LeptonPointer leptonFromW;
	ElectronPointer electronFromW;	// For compatibility with older code
	METPointer met;

	double M3() const;

	bool operator==(const TtbarHypothesis& hyp) const;
	bool operator<(const TtbarHypothesis& hyp) const;

};

typedef boost::shared_ptr<TtbarHypothesis> TtbarHypothesisPointer;

struct compare_totalChi2 {
    bool operator ()(TtbarHypothesisPointer lhs, TtbarHypothesisPointer rhs) {
        return lhs->totalChi2 < rhs->totalChi2;
    }

    bool operator ()(TtbarHypothesis lhs, TtbarHypothesis rhs) {
        return lhs.totalChi2 < rhs.totalChi2;
    }
};

struct compare_disc {
    bool operator ()(TtbarHypothesisPointer lhs, TtbarHypothesisPointer rhs) {
        return lhs->disc <= rhs->disc;
    }

    bool operator ()(TtbarHypothesis lhs, TtbarHypothesis rhs) {
        return lhs.disc <= rhs.disc;
    }
};
} // namespace BAT

#endif /* TTBARHYPOTHESIS_H_ */
