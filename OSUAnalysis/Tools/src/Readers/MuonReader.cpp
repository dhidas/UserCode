/*
 * MuonReader.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#include "../../interface/Readers/MuonReader.h"

namespace BAT {

MuonReader::MuonReader() :
//	numberOfMuonsReader(),
	energyReader(),
	pxReader(),
	pyReader(),
	pzReader(),
	chargeReader(),
	ecalIsolationReader(),
	hcalIsolationReader(),
	trackerIsolationReader(),
	vertex_dist_xy(),
	vertex_dist_z(),
	PFGammaIsolationReader(),
	PFChargedHadronIsolationReader(),
	PFNeutralHadronIsolationReader(),
	d0_PV_Reader(),
	isGlobalReader() {

}

MuonReader::MuonReader(TChainPointer input, MuonAlgorithm::value algo) :
	energyReader(input, MuonAlgorithm::prefixes.at(algo) + ".Energy"),
	pxReader(input, MuonAlgorithm::prefixes.at(algo) + ".Px"),
	pyReader(input, MuonAlgorithm::prefixes.at(algo) + ".Py"),
	pzReader(input, MuonAlgorithm::prefixes.at(algo) + ".Pz"),
	chargeReader(input, MuonAlgorithm::prefixes.at(algo) + ".Charge"),
	ecalIsolationReader(input, MuonAlgorithm::prefixes.at(algo) + ".EcalIso03"),
	hcalIsolationReader(input, MuonAlgorithm::prefixes.at(algo) + ".HcalIso03"),
	trackerIsolationReader(input, MuonAlgorithm::prefixes.at(algo) + ".TrkIso03"),
	vertex_dist_xy(input, MuonAlgorithm::prefixes.at(algo) + ".VtxDistXY"),
	vertex_dist_z(input, MuonAlgorithm::prefixes.at(algo) + ".VtxDistZ"),
	PFGammaIsolationReader(input,MuonAlgorithm::prefixes.at(algo) + ".PFGammaIso"),
	PFChargedHadronIsolationReader(input,MuonAlgorithm::prefixes.at(algo) + ".PfChargedHadronIso"),
	PFNeutralHadronIsolationReader(input,MuonAlgorithm::prefixes.at(algo) + ".PfNeutralHadronIso"),
	d0_PV_Reader(input, MuonAlgorithm::prefixes.at(algo) + ".dB"),
	isGlobalReader(input, MuonAlgorithm::prefixes.at(algo) + ".isGoodGlobalMuon") {

}

void MuonReader::initialise() {
    energyReader.initialise();
    pxReader.initialise();
    pyReader.initialise();
    pzReader.initialise();
    chargeReader.initialise();

    ecalIsolationReader.initialise();
    hcalIsolationReader.initialise();
    trackerIsolationReader.initialise();

    isGlobalReader.initialise();
    d0_PV_Reader.initialise();
		PFGammaIsolationReader.initialise();
		PFChargedHadronIsolationReader.initialise();
		PFNeutralHadronIsolationReader.initialise();
    vertex_dist_xy.initialise();
    vertex_dist_z.initialise();
}

const MuonCollection& MuonReader::getMuons() {
    if (muons.empty() == false)
        muons.clear();
    readMuons();
    return muons;
}

void MuonReader::readMuons() {
    for (unsigned int index = 0; index < energyReader.size(); index++) {
        float energy = energyReader.getVariableAt(index);
        float px = pxReader.getVariableAt(index);
        float py = pyReader.getVariableAt(index);
        float pz = pzReader.getVariableAt(index);
        MuonPointer muon(new Muon(energy, px, py, pz));
        muon->setCharge(chargeReader.getIntVariableAt(index));
        muon->setEcalIsolation(ecalIsolationReader.getVariableAt(index));
        muon->setHcalIsolation(hcalIsolationReader.getVariableAt(index));
        muon->setTrackerIsolation(trackerIsolationReader.getVariableAt(index));
        muon->makeGlobal(isGlobalReader.getBoolVariableAt(index));
        muon->setD0(d0_PV_Reader.getVariableAt(index));
        muon->setXyDistanceToPrimaryVertex(vertex_dist_xy.getVariableAt(index));
        muon->setZDistanceToPrimaryVertex(vertex_dist_z.getVariableAt(index));
				muon->setPFGammaIsolation(PFGammaIsolationReader.getVariableAt(index));
				muon->setPFChargedHadronIsolation(PFChargedHadronIsolationReader.getVariableAt(index));
				muon->setPFNeutralHadronIsolation(PFNeutralHadronIsolationReader.getVariableAt(index));
        muons.push_back(muon);
    }
}
MuonReader::~MuonReader() {
}

}
