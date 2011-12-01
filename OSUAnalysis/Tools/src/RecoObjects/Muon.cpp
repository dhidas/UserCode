/*
 * Muon.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#include "../../interface/RecoObjects/Muon.h"

namespace BAT {

Muon::Muon() :
	Lepton(), is_Global(false), ecal_Isolation(0.), hcal_Isolation(0), tracker_Isolation(0) {

}

Muon::Muon(float energy, float px, float py, float pz) :
	Lepton(energy, px, py, pz), is_Global(false), ecal_Isolation(0.), hcal_Isolation(0), tracker_Isolation(0) {

}

Muon::~Muon() {
}

bool Muon::isGlobal() const {
	return is_Global;
}

void Muon::makeGlobal(bool global) {
	is_Global = global;
}

float Muon::ecalIsolation() const {
	return ecal_Isolation;
}

void Muon::setEcalIsolation(float isolation) {
	ecal_Isolation = isolation;
}

float Muon::hcalIsolation() const {
	return hcal_Isolation;
}

void Muon::setHcalIsolation(float isolation) {
	hcal_Isolation = isolation;
}

float Muon::trackerIsolation() const {
	return tracker_Isolation;
}

void Muon::setTrackerIsolation(float isolation) {
	tracker_Isolation = isolation;
}

float Muon::relativeIsolation() const {
	return (ecal_Isolation + hcal_Isolation + tracker_Isolation) / pt();
}

bool Muon::isIsolated() const{
	return pfIsolation() < 0.1;
}

bool Muon::isGood() const{
	bool passesPt = pt() > 35;
	bool passesEta = fabs(eta()) < 2.1;

	bool passesD0 = fabs(d0()) < 0.02; //cm

	bool passesDistanceToPV = fabs(zDistanceToPrimaryVertex) < 1.0;
	return passesPt && passesEta && passesD0 && passesDistanceToPV && is_Global;
}


bool Muon::isLoose() const {
    bool passesPt = pt() > 10;
    bool passesEta = fabs(eta()) < 2.5;
		bool looseisolated = pfIsolation() < 0.2;
    return passesPt && passesEta && looseisolated && is_Global;

}

}
