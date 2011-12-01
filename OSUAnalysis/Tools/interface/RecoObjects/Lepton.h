/*
 * Lepton.h
 *
 *  Shared features of electrons & muons
 */

#ifndef LEPTON_H_
#define LEPTON_H_
#include "Particle.h"
#include <vector>
#include <string>
#include <boost/array.hpp>

namespace BAT {

class Lepton: public Particle {
public:
	Lepton();
	Lepton(float energy, float px, float py, float pz);
	virtual ~Lepton();
	float XyDistanceToPrimaryVertex() const;
	float ZDistanceToPrimaryVertex() const;
	float PFGammaIsolation() const;
	float PFChargedHadronIsolation() const;
	float PFNeutralHadronIsolation() const;
	float pfIsolation() const;

	void setXyDistanceToPrimaryVertex (float dist);
	void setZDistanceToPrimaryVertex(float dist);
	void setPFGammaIsolation(float pfGammaIso);
	void setPFChargedHadronIsolation(float chargedHadronIso);
	void setPFNeutralHadronIsolation(float neutralHadronIso);

protected:
	float xyDistanceToPrimaryVertex, zDistanceToPrimaryVertex;
	float PFGamma_Isolation, PFChargedHadron_Isolation, PFNeutralHadron_Isolation;
};

typedef boost::shared_ptr<Lepton> LeptonPointer;

}

#endif /* LEPTON_H_ */
