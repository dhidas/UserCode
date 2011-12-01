/*
 * Lepton.cpp
 *
 */

#include "../../interface/RecoObjects/Lepton.h"

namespace BAT {

static const float initialBigValue = 123456789;

Lepton::Lepton() :
	Particle(),
	xyDistanceToPrimaryVertex(initialBigValue),
	zDistanceToPrimaryVertex(initialBigValue),
	PFGamma_Isolation(initialBigValue),
	PFChargedHadron_Isolation(initialBigValue),
	PFNeutralHadron_Isolation(initialBigValue)
{
}

Lepton::Lepton(float energy, float px, float py, float pz) :
	Particle(energy, px, py, pz),
	xyDistanceToPrimaryVertex(initialBigValue),
	zDistanceToPrimaryVertex(initialBigValue),
	PFGamma_Isolation(initialBigValue),
	PFChargedHadron_Isolation(initialBigValue),
	PFNeutralHadron_Isolation(initialBigValue)
{
}

Lepton::~Lepton() {
}

void Lepton::setXyDistanceToPrimaryVertex (float dist) {
    xyDistanceToPrimaryVertex = dist;
}

void Lepton::setZDistanceToPrimaryVertex(float dist) {
    zDistanceToPrimaryVertex = dist;
}

void Lepton::setPFGammaIsolation(float pfGammaIso) {
    PFGamma_Isolation = pfGammaIso;
}

void Lepton::setPFChargedHadronIsolation(float chargedHadronIso) {
    PFChargedHadron_Isolation = chargedHadronIso;
}

void Lepton::setPFNeutralHadronIsolation(float neutralHadronIso) {
    PFNeutralHadron_Isolation = neutralHadronIso;
}

float Lepton::PFGammaIsolation() const {
    return PFGamma_Isolation;
}

float Lepton::PFChargedHadronIsolation() const {
    return PFChargedHadron_Isolation;
}

float Lepton::PFNeutralHadronIsolation() const {
    return PFNeutralHadron_Isolation;
}

float Lepton::XyDistanceToPrimaryVertex () const {
    return xyDistanceToPrimaryVertex ;
}

float Lepton::ZDistanceToPrimaryVertex() const {
    return zDistanceToPrimaryVertex;
}

float Lepton::pfIsolation() const {
    return (PFGamma_Isolation + PFChargedHadron_Isolation + PFNeutralHadron_Isolation) / pt();
		// Previously was divided by et, but I think pt is correct.
}

}
