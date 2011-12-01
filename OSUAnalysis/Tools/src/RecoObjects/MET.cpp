/*
 * MET.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#include "../../interface/RecoObjects/MET.h"

namespace BAT {


MET::MET() :
    Particle(),
    usedAlgorithm(METAlgorithm::Calo) {
}

MET::MET(float ex, float ey) :
    Particle(sqrt(ex * ex + ey * ey), ex, ey, 0), usedAlgorithm(METAlgorithm::Calo) {

}

MET::~MET() {
}

bool MET::isGood() const {
    return et() > 20;
}

void MET::setUsedAlgorithm(METAlgorithm::value algo) {
    usedAlgorithm = algo;
}

METAlgorithm::value MET::getUsedAlgorithm() const {
    return usedAlgorithm;
}


void MET::adjForUnc(const JetCollection &jets)
{
	if (Jet::correctDirection == JetCorrDirection::NONE ||
		jets.size() <= 0)
		return;
	FourVector new4vec = getFourVector();
	for (unsigned int index = 0; index < jets.size(); ++index)
		new4vec += jets.at(index)->DiffVec;
	new4vec.SetPz(0.0);  // Make sure no longitudinal component
	new4vec.SetE(new4vec.Pt()); // Make sure it stays massless
	setFourVector(new4vec);
}



}
