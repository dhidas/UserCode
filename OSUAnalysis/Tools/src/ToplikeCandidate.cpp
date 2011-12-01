#include <iostream>
#include <iomanip>

#include "../interface/ToplikeCandidate.h"


using namespace std;
using namespace BAT;

namespace BAT {

static const double pi = 3.1415927;
static const int DQ_ID = 1, BQ_ID = 5, TOP_ID = 6, W_ID = 24, ELECTRON_ID = 11;
static const int TAU_ID = 15, NU_E_ID = 12, MU_ID = 13, NU_MU_ID = 14;
static const int WPRIME_ID = 4000023;
static const unsigned int MCLIMIT = 90;  // Number of MC particles to check
static const float DELTA_R_LIMIT = 0.5; // Biggest deltaR to consider in range
static const double pt_tprime_sys_mean = 48.49;
static const double pt_tprime_sys_sigma = 27.87;

static map <ToplikeSelectionSteps::Step, TTbarEPlusJetsSelection::Step> selSteps;


ToplikeCandidate::ToplikeCandidate() :
	jet3FromWIndex(0),
  jet4FromWIndex(0),
	doneReconstructiontop(false)
{
}
		
	
ToplikeCandidate::ToplikeCandidate(const Event& event) :
	TopPairEventCandidate(event),
	jet3FromWIndex(0),
  jet4FromWIndex(0),
	doneReconstructiontop(false)
{
	static bool mapset = false;
	if (mapset == false) {
		for(unsigned int ind = 0; ind < sizeof(selStepArr)/sizeof(selPair); ++ind) {
			selSteps[selStepArr[ind].key] = selStepArr[ind].value;
		}
		mapset= true;
	}
}


ToplikeCandidate::~ToplikeCandidate()
{
}


bool ToplikeCandidate::hasAtLeastFiveGoodJets() const {
    return goodJets.size() >= 5;
}

/*
const ParticlePointer ToplikeCandidate::getResonance() const {
    throwExpeptionIfNotReconstructed("ToplikeCandidate::getResonance");
    return tPrime;
}
*/


double ToplikeCandidate::tpmass() const {
		/*
		ParticlePointer tmp = ParticlePointer(new Particle(*leptonicTop1 + *hadronicW));
		cout << "tmp px = " << tmp->px() << " tp py = " << tmp->py() << " tp pz = " << tmp->pz() << endl;
		cout << "energy = " << tmp->energy() << endl;
    return tmp->mass();
		cout << "tp px = " << tPrime->px() << " tp py = " << tPrime->py() << " tp pz = " << tPrime->pz() << endl;
		cout << "energy = " << tPrime->energy() << endl;
		if (tPrime->px() <= 0)
			tPrime->setMass(0.0);
		*/
    return tPrime->mass();
}


double ToplikeCandidate::sumPt() const {
    return (jet1FromW->pt() + jet2FromW->pt() + hadronicBJet->pt());
}


double ToplikeCandidate::getGlobalChi2() const {
    double pttbar = PtTtbarSystem();
    double pttbarTerm = TMath::Power(pttbar - pt_tprime_sys_mean, 2) / (2 *  pt_tprime_sys_sigma * pt_tprime_sys_sigma);

    double htSystemDifference = TMath::Power(HTSystem() - matched_HTSystem, 2);
    double htSystemError = matched_HTSystem_sigma * matched_HTSystem_sigma * 2;
    double htSystemTerm = htSystemDifference / htSystemError;
    return 1 / sqrt(2) * (pttbarTerm + htSystemTerm);
}


double ToplikeCandidate::PtTtbarSystem() const {
    ParticlePointer combined;
		combined = ParticlePointer(new Particle(*hadronicWtPrime + *hadronicTop));
    return combined->pt() / HTSystem();
}


double ToplikeCandidate::TPrimeHTSystem() const {
	return htSystem;
}


double ToplikeCandidate::PtTPrimeSystem() const {
	return ptTprimeSystem;
}


double ToplikeCandidate::PtTtbarSystem(unsigned short neutrinoSolution) const {
    ParticlePointer combined;
    if (neutrinoSolution == 1)
        combined = ParticlePointer(new Particle(*leptonicTop1 + *hadronicW));
    else
        combined = ParticlePointer(new Particle(*leptonicTop2 + *hadronicW));
    return combined->pt() / HTSystem();
}


double ToplikeCandidate::getWChi2() const {
	double WmassDifference = TMath::Power(hadronicWtPrime->mass() - pdgWmass, 2);
	double WmassError = 2 * wMassWidth * wMassWidth;
	double WmassTerm = WmassDifference / WmassError;
	return (WmassTerm);
}


double ToplikeCandidate::getLoneHadChi2() const
{
	return (TopPairEventCandidate::getHadronicChi2());
}


/*
double ToplikeCandidate::getTotalChi2(unsigned short int neutrinoSolution) const {
    return getLeptonicChi2(neutrinoSolution) + getHadronicChi2() + getGlobalChi2(neutrinoSolution);
}
*/


void ToplikeCandidate::recoTprimeUsingChi2(ElectronPointer electron) {
  // cout << "Starting tp reco " << endl;
	if (goodJets.size() < 3)
		throw ReconstructionException("Not enough jets available to reconstruct Tprime using Chi2 method.");
	double chosen_Chi2Total(9999999.);
	leptonFromW = electron;
  selectedNeutrino = 0;
  currentSelectedNeutrino = 0;
  reconstructNeutrinos();

  // cout << "Starting jet loop" << endl;
	for (unsigned short hadBindex = 0; hadBindex < goodJets.size(); ++hadBindex) {
		for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
			if (jet1Index == hadBindex)
				continue;
			for (unsigned short jet2Index = 0; jet2Index < goodJets.size(); ++jet2Index) {
				if (jet2Index == jet1Index || jet2Index == hadBindex)
					continue;
				hadronicBJet = goodJets.at(hadBindex);
				leptonicBJet = hadronicBJet; // hack
				jet1FromW = goodJets.at(jet1Index);
				jet2FromW = goodJets.at(jet2Index);
        leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
        leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
        leptonicTop1 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW1));
        leptonicTop2 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW2));
				hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
				// fillHypotheses();
				selectNeutrinoSolution();
				double chi2 = TopPairEventCandidate::getTotalChi2(currentSelectedNeutrino);
				if (chi2 < chosen_Chi2Total) {
						hadronicBIndex = hadBindex;
						jet1FromWIndex = jet1Index;
						jet2FromWIndex = jet2Index;
						chosen_Chi2Total = chi2;
            selectedNeutrino = currentSelectedNeutrino;
				}
			}
		}
	}
	// cout << "Done with loop" << endl;
	hadronicBJet = goodJets.at(hadronicBIndex);
	leptonicBJet = hadronicBJet; // hack
	jet1FromW = goodJets.at(jet1FromWIndex);
	jet2FromW = goodJets.at(jet2FromWIndex);
	hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
  leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
  leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
  leptonicTop1 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW1));
  leptonicTop2 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW2));
  ptTprimeSystem = PtTtbarSystem(selectedNeutrino);
  htSystem = HTSystem();
  if (selectedNeutrino == 1)
		tPrime = ParticlePointer(new Particle(*leptonicTop1 + *hadronicW));
	else 
		tPrime = ParticlePointer(new Particle(*leptonicTop2 + *hadronicW));
  // cout << "tPrime mass = " << tPrime->mass() << endl;
	doneReconstruction = true;
}


// Omit lepton-jet angle from criteria.
double ToplikeCandidate::getLeptonicChi2(double topMass, double unused)
	const
{
    double massDifference = TMath::Power(topMass - pdgTopMass, 2);
    double massError = 2 * topMassWidth * topMassWidth;
    double massTerm = massDifference / massError;
    return (massTerm);
}


double ToplikeCandidate::getHadronicChi2() const {
    double WmassDifference = TMath::Power(hadronicW->mass() - pdgWmass, 2);
    double WmassError = 2 * wMassWidth * wMassWidth;
    double WmassTerm = WmassDifference / WmassError;

    double topMassDifference = TMath::Power(hadronicTop->mass() - pdgTopMass, 2);
    double topMassError = 2 * topMassWidth * topMassWidth;
    double topMassTerm = topMassDifference / topMassError;
		return ((1.0 / sqrt(2.0)) * (topMassTerm + WmassTerm));
}


void ToplikeCandidate::recoBestSingleTop(LeptonPointer lepton) {
	if (goodJets.size() < 3)
		throw ReconstructionException("Not enough jets available to reconstruct top using Chi2 method.");
	double bestLepChi2(9999999.), bestHadChi2(9999999.);
	leptonFromW = lepton;
  selectedNeutrino = 0;
  currentSelectedNeutrino = 0;
  reconstructNeutrinos();

	for (unsigned short hadBindex = 0; hadBindex < goodJets.size(); ++hadBindex) {
		for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
			if (hadBindex == jet1Index)
				continue;
			for (unsigned short jet2Index = 0; jet2Index < goodJets.size(); ++jet2Index) {
				if (jet2Index == jet1Index || jet2Index == hadBindex)
					continue;
				hadronicBJet = goodJets.at(hadBindex);
				leptonicBJet = hadronicBJet; // hack
				jet1FromW = goodJets.at(jet1Index);
				jet2FromW = goodJets.at(jet2Index);
				leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
				leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
				leptonicTop1 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW1));
				leptonicTop2 = ParticlePointer(new Particle(*hadronicBJet + *leptonicW2));
				hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
				hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
				// fillHypotheses();
				selectNeutrinoSolution();
				double ltopchi2 = TopPairEventCandidate::getLeptonicChi2(currentSelectedNeutrino);
				double htopchi2 = getHadronicChi2();
				if (ltopchi2 < bestLepChi2) {
					selectedNeutrino = currentSelectedNeutrino;
					leptonicBIndex = hadBindex;
					bestLepChi2 = ltopchi2;
				}
				if (htopchi2 < bestHadChi2) {
					hadronicBIndex = hadBindex;
					jet1FromWIndex = jet1Index;
					jet2FromWIndex = jet2Index;
					bestHadChi2 = htopchi2;
				}
			}
		}
	}
	hadronicBJet = goodJets.at(hadronicBIndex);
	leptonicBJet = goodJets.at(leptonicBIndex);
	jet1FromW = goodJets.at(jet1FromWIndex);
	jet2FromW = goodJets.at(jet2FromWIndex);
	hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
  leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
  leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
  hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
	leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
	leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
  // cout << "Found lone tops " << endl;
 	doneReconstructiontop = true;
 	doneReconstruction = true;
}



static bool cmpMcEntries(const mcObj &firstArg, const mcObj &secondArg)
{
	return (firstArg.deltaR < secondArg.deltaR);
}


mcList ToplikeCandidate::getMCList(const Particle *const recoObject) const
{
	const MCParticleCollection &mcColl = genParticles;
	mcList mcObjs;
	if (mcColl.size() > 0) {
		for (unsigned int i = 0; i < mcColl.size(); ++i) {
			mcObj mcEntry;
			mcEntry.deltaR = recoObject->deltaR(mcColl.at(i));
			mcEntry.ptr = mcColl.at(i);
			mcEntry.gpInd = i;
			mcObjs.push_back(mcEntry);
		}
		sort(mcObjs.begin(), mcObjs.end(), cmpMcEntries);
	}
	return (mcObjs);
}


mcList ToplikeCandidate::getMCListForJets(const Jet *const recoObject)
	const
{
	const JetCollection &mcColl = genJets;
	mcList mcObjs;
	if (mcColl.size() > 0) {
		for (unsigned int i = 0; i < mcColl.size(); ++i) {
			mcObj mcEntry;
			mcEntry.deltaR = recoObject->deltaR(mcColl.at(i));
			mcEntry.gpInd = i;
			mcObjs.push_back(mcEntry);
		}
		sort(mcObjs.begin(), mcObjs.end(), cmpMcEntries);
	}
	return (mcObjs);
}


int ToplikeCandidate::missRecoJetsForGenJets() const
{
	const JetCollection &mcColl = genJets;
	int numGenJets = mcColl.size();
	int matchJets = 0;
	for (int i = 0; i < numGenJets; ++i) {
		bool match = false;
		for (unsigned int recoInd = 0;
				recoInd < GoodJets().size() && match == false; ++recoInd)
			if (GoodJets().at(recoInd)->deltaR(mcColl.at(i)) < 0.3) {
					++matchJets;
					match = true;
			}
	}
	return (numGenJets - matchJets);
}



topTruth ToplikeCandidate::getMCMatches(const ParticlePointer &recoObject)
const
{
	const float nullDeltaR = 3;
	bool deltaRSet = false;
	topTruth result = {nullDeltaR, false, false};
	mcList mcObjs = getMCList(recoObject.get());
	if (mcObjs.size() <= 0)
		return (result);
	/*
	if (mcObjs.begin()->first > 0.9)  // Nothing close enough
		return ("");
	TString idList;
  idList += mcObjs.at(0).first;
	idList += " ";
	for (unsigned int i = 0; i < mcObjs.size() && i < 12; ++i) {
		idList += PdgStrs::pdgStr(mcObjs.at(i).second);
		idList += " ";
	}
	*/
	float topDeltaR = nullDeltaR;
	bool wFound = false, lepFound = false;
	for (unsigned int i = 0; i < mcObjs.size() && i < MCLIMIT &&
		(deltaRSet == false || lepFound == false || wFound == false); ++i) {
		if (deltaRSet == false && fabs(mcObjs.at(i).ptr->pdgId()) == TOP_ID) {
			topDeltaR = mcObjs.at(i).deltaR;
			deltaRSet = true;
			if (topDeltaR > nullDeltaR)
				topDeltaR = nullDeltaR;
		}
		if (fabs(mcObjs.at(i).ptr->pdgId()) == ELECTRON_ID || 
				fabs(mcObjs.at(i).ptr->pdgId()) == MU_ID)
			lepFound = true;
		if (fabs(mcObjs.at(i).ptr->pdgId()) == W_ID)
			wFound = true;
	}
	if (deltaRSet == false)
		cout << " No top found in GPs! Increase MCLIMIT.\n" ;
	result.deltaR = topDeltaR;
	result.leptonic = wFound && lepFound;
	result.hadronic = wFound;
	return (result);
}


static int findMCParticle(const mcList &mcObjs, int pdgId, int antiId = 0,
	unsigned int startInd = 0)
{
	if (antiId == 0)
		antiId = pdgId;
	for (unsigned int i = startInd; i < mcObjs.size(); ++i) {
		if (pdgId == 0 || mcObjs.at(i).ptr->pdgId() == pdgId ||
				mcObjs.at(i).ptr->pdgId() == antiId) {
			if (mcObjs.at(i).deltaR < DELTA_R_LIMIT)
				return (i);
			else return (-1);
		}
	}
	return (-1);
}


// Scans up the decay chain when mother & daughter are the same.

bool ToplikeCandidate::chkIfChild(int gpInd, int motherGpInd) const
{
	if (motherGpInd < 0)
		return (false);
	const MCParticleCollection &mcColl = genParticles;
	if (mcColl.size() <= 0)
		return (false);
	int mother = mcColl.at(gpInd)->motherIndex();
	if (mother == motherGpInd)
		return (true);
	while (mother > 0 &&
		mcColl.at(gpInd)->pdgId() == mcColl.at(mother)->pdgId())
	{
		gpInd = mother;
		mother = mcColl.at(gpInd)->motherIndex();
		if (mother == motherGpInd)
			return (true);
	}
	int gMother = mcColl.at(motherGpInd)->motherIndex();

	// W' sample b's have W' as mother rather than top.
	if (gMother > 0 &&
			fabs(mcColl.at(gMother)->pdgId()) == WPRIME_ID && mother == gMother)
		return (true);
	while (gMother > 0 &&
		mcColl.at(motherGpInd)->pdgId() == mcColl.at(gMother)->pdgId())
	{
		motherGpInd = gMother;
		gMother = mcColl.at(motherGpInd)->motherIndex();
		if (mother == motherGpInd)
			return (true);
	}
	/*
	if (mcColl.at(mother)->motherIndex() == mcColl.at(motherGpInd)->motherIndex())
	{
		cout << "gm match: ind cpdg mi mpdg gm pdg omi ompdg " << gpInd << " ";
		cout << mcColl.at(gpInd)->pdgId() << " " << mother << " ";
		cout << mcColl.at(mother)->pdgId() << " ";
		int gm = mcColl.at(mother)->motherIndex();
		cout << gm << " " << mcColl.at(gm)->pdgId() << " ";
		cout << motherGpInd << " " << mcColl.at(motherGpInd)->pdgId() << endl;
		return (true);
	}
	*/
	return (false);
}


int ToplikeCandidate::chkMCMatch(const Particle *const recoObject, 
	int motherInd, const int pdgId, int antiId) const
{
	mcList truList = getMCList(recoObject);
	if (truList.size() > 0) {
		int index = -1;
		while ((index = findMCParticle(truList, pdgId, antiId, index + 1)) >= 0)
			if (motherInd == -1 || chkIfChild(truList.at(index).gpInd, motherInd))
				return (truList.at(index).gpInd);
	}
	return (-1);
}


// Look for a non-top quark whose mother is specified by motherInd and which
// doesn't match badInd.  If motherInd == -1, no mother is checked for.

int ToplikeCandidate::chkForQuark(const JetPointer recoObject, 
	const int motherInd, const int badInd) const
{
	mcList truList = getMCList(recoObject.get());
	if (truList.size() > 0) {
		for (int i = 0; i < (int) truList.size() &&
				truList.at(i).deltaR < DELTA_R_LIMIT; ++i) {
			if (truList.at(i).ptr->isQuark() && fabs(truList.at(i).ptr->pdgId()) !=
					TOP_ID && (motherInd == -1 || chkIfChild(truList.at(i).gpInd, motherInd)) && truList.at(i).gpInd != badInd)
				return (truList.at(i).gpInd);
				// Checking for q-qbar pair not necessary since both daughters of mother
		}
	}
	return (-1);
}


// Checks if ind1 is a b quark, other two are daughters of a W, and
// b & W are siblings.
bool ToplikeCandidate::chkBW(const int ind1, const int ind2, const int ind3)
	const
{
	const MCParticleCollection &mcColl = genParticles;
	if (mcColl.size() <= 0)
		return (false);
	if (fabs(mcColl.at(ind1)->pdgId()) != BQ_ID)
		return (false);
	int jetMother = mcColl.at(ind2)->motherIndex();
	if (jetMother > 0 && fabs(mcColl.at(jetMother)->pdgId()) == W_ID &&
			mcColl.at(ind3)->motherIndex() == jetMother && 
			mcColl.at(ind1)->motherIndex() == mcColl.at(jetMother)->motherIndex())
		return (true);
	return (false);
}


int ToplikeCandidate::goodHTopMatching() const
{
	int numFound = 0, mostFound = 0;
	mcList bTruList = getMCList(hadronicBJet.get());
	mcList j1TruList = getMCList(jet1FromW.get());
	mcList j2TruList = getMCList(jet2FromW.get());
	int bInd = -1;
	while (bTruList.size() > 0 &&
			(bInd = findMCParticle(bTruList, 0, 0, bInd + 1)) >= 0)
	{
		if (bTruList.at(bInd).ptr->isQuark())
		{
			numFound = 1;
			int bGPInd = bTruList.at(bInd).gpInd;
			int j1Ind = -1;
			while (j1TruList.size() > 0 &&
					(j1Ind = findMCParticle(j1TruList, 0, 0, j1Ind + 1)) >= 0)
			{
				int j1GPInd = j1TruList.at(j1Ind).gpInd;
				if (j1GPInd != bGPInd && j1TruList.at(j1Ind).ptr->isQuark())
				{
					numFound = 2;
					int j2Ind = -1;
					while (j2TruList.size() > 0 &&
							(j2Ind = findMCParticle(j2TruList, 0, 0, j2Ind + 1)) >= 0)
					{
						int j2GPInd = j2TruList.at(j2Ind).gpInd;
						if (j2GPInd != bGPInd && j2GPInd != j1GPInd &&
							j2TruList.at(j2Ind).ptr->isQuark())
						{
							numFound = 3;
							if (chkBW(bGPInd, j1GPInd, j2GPInd) ||
									chkBW(j1GPInd, bGPInd, j2GPInd) ||
									chkBW(j2GPInd, j1GPInd, bGPInd))
								return (4);
						}
					}
					if (numFound > mostFound)
						mostFound = numFound;
				}
			}
		}
	}
	if (numFound > mostFound)
		mostFound = numFound;
	return (mostFound);
}



int ToplikeCandidate::strictHTopMatching(const mcObj &topTru) const
{
	int numFound = 1; 	// Start with finding the top.
	short antiMult = 1;
	if (topTru.ptr->pdgId() < 0)
		antiMult = -1;
	// cout << " hti " << topTru.gpInd << " ";
	if (chkMCMatch(hadronicBJet.get(), topTru.gpInd, BQ_ID * antiMult) >= 0)
		++numFound;
	int wInd = chkMCMatch(hadronicW.get(), topTru.gpInd, W_ID * antiMult);
	if (wInd >= 0)
		++numFound;
	else return (numFound);

	// Try to find 2 quarks from W.  If both fail, trying again will still fail
	// but is harmless.  If one fails, try again to see if can get both if we
	// let the other jet grab the 1st quark.  If that doesn't work, one quark
	// will still be found.
	bool pairFnd = false, last1 = false;
	int takenInd2 = -1;
	int takenInd = chkForQuark(jet1FromW, wInd);
	if (takenInd >= 0 && chkForQuark(jet2FromW, wInd, takenInd) >= 0) {
		numFound += 2;
		pairFnd = true;
	} else {
		takenInd2 = chkForQuark(jet2FromW, wInd);
		if (takenInd2 >= 0)
			++numFound;
		if (chkForQuark(jet1FromW, wInd, takenInd2) >= 0) {
			++numFound;
			last1 = true;
		}
	}
	/*
	if (numFound == 5 || numFound == 4) {
		cout << "num pair wInd takenInd taken2 last1 " << numFound << " ";
		cout << pairFnd << " " << wInd << " ";
		cout << takenInd << " " << takenInd2 << " " << last1 << endl;
	}
	*/
	return (numFound);
}


int ToplikeCandidate::strictLTopMatching(const mcObj &topTru) const
{
	int numFound = 1;
	short antiMult = 1;
	if (topTru.ptr->pdgId() < 0)
		antiMult = -1;
	if (chkMCMatch(getLeptonicBJet().get(), topTru.gpInd, BQ_ID * antiMult) >= 0)
		++numFound;
	int wInd = chkMCMatch(getLeptonicW().get(), topTru.gpInd, W_ID * antiMult);
	if (wInd >= 0)
		++numFound;
	else return (numFound);
	// Lepton is opposite from top to preserve lepton number.
	if (chkMCMatch(getLeptonFromWDecay().get(), wInd, -ELECTRON_ID * antiMult)
			>= 0 || 
			chkMCMatch(getLeptonFromWDecay().get(), wInd, -MU_ID * antiMult) >= 0)
		++numFound;
	if (chkMCMatch(getNeutrinoFromWDecay().get(), wInd, NU_E_ID * antiMult) >= 0
			|| chkMCMatch(getNeutrinoFromWDecay().get(), wInd, NU_MU_ID * antiMult)
			>= 0)
		++numFound;
	return (numFound);
}


int ToplikeCandidate::goodLTopMatching() const
{
	const MCParticleCollection &mcColl = genParticles;
	int numParticles = mcColl.size();
	if (numParticles <= 0)
		return (0);
	int numFound = 0, mostFound = 0;
	mcList lepTruList = getMCList(getLeptonFromWDecay().get());
	mcList bTruList = getMCList(getLeptonicBJet().get());
	for (int lepTypes = 0; lepTypes < 2; ++lepTypes) {
		int lepInd = -1;
		int lepPdgId = ELECTRON_ID;
		if (lepTypes == 1)
			lepPdgId = MU_ID;
		while (lepTruList.size() > 0 && (lepInd =
				findMCParticle(lepTruList, lepPdgId, -lepPdgId, lepInd + 1)) >= 0)
		{
			numFound = 1;
			int bInd = -1;
			int lepGPInd = lepTruList.at(lepInd).gpInd;
			if (lepGPInd >= 0 && lepGPInd < numParticles) {
				int lepMothInd = mcColl.at(lepGPInd)->motherIndex();
				if (lepMothInd >= 0 && lepMothInd < numParticles && 
						fabs(mcColl.at(lepMothInd)->pdgId()) == W_ID) {
					numFound = 2;
					while (bTruList.size() > 0 &&
							(bInd = findMCParticle(bTruList, BQ_ID, -BQ_ID, bInd + 1)) >= 0)
					{
						numFound = 3;
						int bGPInd = bTruList.at(bInd).gpInd;
						if ((bGPInd >= 0 && bGPInd < numParticles) &&
								(chkIfChild(bGPInd, mcColl.at(lepMothInd)->motherIndex()) ||
								chkIfChild(lepMothInd, mcColl.at(bGPInd)->motherIndex())))
							return(4);
					}
				}
			}
			if (numFound > mostFound)
				mostFound = numFound;
		}
	}
	return (mostFound);
}



/*
int ToplikeCandidate::getNumMCPermutsForTop(const ParticlePointer &recoTop, 
	bool hadronic) const
{
	mcList topTru = getMCList(recoTop.get());
	if (topTru.size() <= 0)
		return (0);
	int highestNumFound = 0;
	int topInd = -1;
	while ((topInd = findMCParticle(topTru, TOP_ID, -TOP_ID, topInd + 1)) >= 0) {
		int numFound = 0;
		if (hadronic)
			 numFound = hTopDauMatching(topTru.at(topInd));
		else numFound = lTopDauMatching(topTru.at(topInd));
		if (numFound > highestNumFound)
			highestNumFound = numFound;
	}
	return (highestNumFound);
}
*/


int ToplikeCandidate::getNumMCMatchesForTop(const ParticlePointer &recoTop, 
	bool hadronic)
	const
{
	mcList topTru = getMCList(recoTop.get());
	if (topTru.size() <= 0)
		return (0);
	int highestNumFound = 0;
	int topInd = -1;
	while ((topInd = findMCParticle(topTru, TOP_ID, -TOP_ID, topInd + 1)) >= 0) {
		int numFound = 0;
		if (hadronic)
			 numFound = strictHTopMatching(topTru.at(topInd));
		else numFound = strictLTopMatching(topTru.at(topInd));
		if (numFound > highestNumFound)
			highestNumFound = numFound;
	}
	return (highestNumFound);
}


int ToplikeCandidate::getNumMCMatchesHTop() const
{
	return (getNumMCMatchesForTop(getHadronicTop(), true));
}

/*
int ToplikeCandidate::getNumDauMatchesHTop() const
{
	return (getNumMCPermutsForTop(getHadronicTop(), true));
}


int ToplikeCandidate::getNumDauMatchesLTop() const
{
	return (getNumMCPermutsForTop(getLeptonicTop(), false));
}
*/


int ToplikeCandidate::getNumCorrectIDHTop() const
{
	mcList topTru = getMCList(getHadronicTop().get());
	int numFound = 0;
	bool topFnd = false, bFnd = false, jet1Fnd = false, jet2Fnd = false;
	if (topTru.size() > 0 && findMCParticle(topTru, TOP_ID, -TOP_ID) >= 0) {
		topFnd = true;
		++numFound;
	}
	
	// Take either particle or anti-particle.
	if (chkMCMatch(hadronicBJet.get(), -1, BQ_ID, -BQ_ID) >= 0) {
		++numFound;
		bFnd = true;
	}
	int wInd = chkMCMatch(hadronicW.get(), -1, W_ID, -W_ID);
	if (wInd >= 0)
		++numFound;

	// Look for any quark near a jet.  Even allow the same quark for both jets.
	if (chkForQuark(jet1FromW, -1) >= 0) {
			++numFound;
			jet1Fnd = true;
	}
	if (chkForQuark(jet2FromW, -1) >= 0) {
		jet2Fnd = true;
		++numFound;
	}
	/*
	if (numFound == 4 || numFound == 3) {
		cout << "had num t b w j1 j2 " << numFound << " " << topFnd << " " << bFnd << " " << wInd;
		cout << " " << jet1Fnd << " " << jet2Fnd << endl;
	}
	*/
	return (numFound);
}


int ToplikeCandidate::getNumMCMatchesLTop() const
{
	return (getNumMCMatchesForTop(getLeptonicTop(), false));
}


int ToplikeCandidate::getNumCorrectIDLTop() const
{
	mcList topTru = getMCList(getLeptonicTop().get());
	int numFound = 0;
	if (topTru.size() > 0 && findMCParticle(topTru, TOP_ID, -TOP_ID) >= 0)
		++numFound;
	if (chkMCMatch(getLeptonicBJet().get(), -1, BQ_ID, -BQ_ID) >= 0)
		++numFound;
	if (chkMCMatch(getLeptonicW().get(), -1, W_ID, -W_ID) >= 0)
		++numFound;
	if (chkMCMatch(getLeptonFromWDecay().get(), -1, -ELECTRON_ID,
			ELECTRON_ID) >= 0 ||
			chkMCMatch(getLeptonFromWDecay().get(), -1, -MU_ID,	MU_ID) >= 0)
		++numFound;
	if (chkMCMatch(getNeutrinoFromWDecay().get(), -1, NU_E_ID, -NU_E_ID) >= 0 ||
			chkMCMatch(getNeutrinoFromWDecay().get(), -1, NU_MU_ID, -NU_MU_ID) >= 0)
		++numFound;
	return (numFound);
}


int ToplikeCandidate::getTausInJets() const
{
	int numTaus = 0;
	for (unsigned int index = 0; index < GoodJets().size(); ++index) {
		mcList truList = getMCList(GoodJets().at(index).get());
		int tauIndex = -1;
		if (truList.size() > 0 &&
				(tauIndex = findMCParticle(truList, TAU_ID, -TAU_ID)) >= 0)
			++numTaus;
		}
	// cout << " Found taus: " << numTaus;
	return (numTaus);
}


int ToplikeCandidate::numTaus() const
{
	const MCParticleCollection &mcColl = genParticles;
	int numTaus = 0;
	if (mcColl.size() > 0) {
		for (unsigned int i = 0; i < mcColl.size(); ++i) {
			if (fabs(mcColl.at(i)->pdgId()) == TAU_ID)
				++numTaus;
		}
	}
	return (numTaus);
}


int ToplikeCandidate::missGenJetsForRecoJets() const
{
	int matchJets = 0, numJets = GoodJets().size();
	for (int index = 0; index < numJets; ++index) {
		mcList mcObjs = getMCListForJets(GoodJets().at(index).get());
		if (mcObjs.size() > 0 && mcObjs.at(0).deltaR < 0.3)
			++matchJets;
	}
	return (numJets - matchJets);
}



TopPlusXCandidates::TopPlusXCandidates() :
	dJetFromWpIndex(-1)
{
}
		
	
TopPlusXCandidates::TopPlusXCandidates(const Event& event) :
	ToplikeCandidate(event),
	dJetFromWpIndex(-1)
{
}


TopPlusXCandidates::~TopPlusXCandidates()
{
}


const JetPointer TopPlusXCandidates::nulljet;


bool TopPlusXCandidates::passesMETCut() const {
    return met->pt() > 20.;
}


bool TopPlusXCandidates::passesHTCut() const {
    return fullHT() > 700.;
}


bool TopPlusXCandidates::hasAtLeastOneGoodJet() const {
    return goodJets.size() >= 1 && goodJets.front()->pt() > 180.0;
}


bool TopPlusXCandidates::hasAtLeastTwoGoodJets() const {
    return goodJets.size() >= 2 && goodJets.at(1)->pt() > 90.0;
}


bool TopPlusXCandidates::hasOnlyOneGoodIsolatedMuon() const {
		return goodIsolatedMuons.size() == 1;
}


bool TopPlusXCandidates::hasNoGoodIsolatedMuon() const {
		return goodIsolatedMuons.size() <= 0;
}


bool TopPlusXCandidates::hasNoGoodIsolatedElectron() const {
    if(Event::usePFIsolation)
        return goodPFIsolatedElectrons.size() <= 0;
    else
        return goodIsolatedElectrons.size() <= 0;
}


bool TopPlusXCandidates::hasNoIsolatedElectron() const {
    return looseElectrons.size() <= 0;
}

bool TopPlusXCandidates::passesExtraElectronCut() const{
		ElectronPointer electron;
    if(Event::usePFIsolation) {
        if (goodPFIsolatedElectrons.size() <= 0)
					return (false);
				electron = GoodPFIsolatedElectrons().front();
		} else if (goodIsolatedElectrons.size() > 0) {
			electron = GoodIsolatedElectrons().front();
		} else return (false);
    return (electron->pt() > 45.);
}


bool TopPlusXCandidates::passesSelectionStep(enum
		ToplikeSelectionSteps::Step step) const {
	switch (step) {
    case ToplikeSelectionSteps::AtLeastFiveGoodJets:
			return (hasAtLeastFiveGoodJets());
    case ToplikeSelectionSteps::FortyfiveGeVElectron:
			return (passesExtraElectronCut());
    case ToplikeSelectionSteps::OneIsolatedMuon:
			return (hasOnlyOneGoodIsolatedMuon());
    case ToplikeSelectionSteps::LooseElectronVeto:
			return (hasNoIsolatedElectron());
    case ToplikeSelectionSteps::TightElectronVeto:
			return (hasNoGoodIsolatedElectron());
    case ToplikeSelectionSteps::TightMuonVeto:
			return (hasNoGoodIsolatedMuon());
    case ToplikeSelectionSteps::HT700GeV:
			return (passesHTCut());
		default:
			return (TopPairEventCandidate::passesSelectionStep(selSteps[step]));
	}
}
	

bool TopPlusXCandidates::passesSelectionStepUpTo(unsigned int step,
		const ToplikeSelectionSteps::Step steps[]) const {
	if (step <= 0)
		return (passesSelectionStep(steps[0]));
	else {
		return (passesSelectionStep(steps[step]) &&
			passesSelectionStepUpTo(step - 1, steps));
	}
}


bool TopPlusXCandidates::passesFullEPlusJetSelection() const {
	return (passesSelectionStepUpTo(toplikeElectronSelSize - 1, toplikeElectronSelection));
}


bool TopPlusXCandidates::passesFullMuPlusJetSelection() const {
	return (passesSelectionStepUpTo(toplikeMuonSelSize - 1, toplikeMuonSelection));
}


bool TopPlusXCandidates::passesFullEPlusJetSelectionKK() const {
	return (passesSelectionStepUpTo(toplikeElectronSelSize - 1, toplikeElectronSelection));
}


double TopPlusXCandidates::jetsHTNoTop(unsigned int numJets) const {
    double ht = 0.0;
		JetCollection topJets;

    for (unsigned int index = 0; index < goodJets.size(); ++index) {
        if (goodJets.at(index) != leptonicBJet) {
						if (topJets.size() <= 0 || topJets.back()->pt() > goodJets.at(index)->pt())
							topJets.push_back(goodJets.at(index));
						else {
							unsigned int srchind = 0;
							while (topJets.at(srchind)->pt() > goodJets.at(index)->pt())
									++srchind;
							topJets.insert(topJets.begin() + srchind, goodJets.at(index));
						}
				}
    }
    for (unsigned int index = 0; index < topJets.size() && index < numJets; ++index) {
				ht += topJets.at(index)->pt();
    }
    return ht;
}


double TopPlusXCandidates::mass2jets(unsigned int jet1, unsigned int jet2) const
{
		int hijets[3] = {-1, -1, -1};
		unsigned int hijindex = 0;

    for (unsigned int index = 0; index < goodJets.size(); ++index) {
			if (goodJets.at(index) != leptonicBJet && hijets[hijindex] == -1) {
				hijets[hijindex] = index;
				if (hijindex < 2)
					hijindex++;
			}
		}
		if (hijets[jet1] > -1 && hijets[jet2] > -1) {
				ParticlePointer candidate =
					ParticlePointer(new Particle(*goodJets.at(hijets[jet1]) +
					*goodJets.at(hijets[jet2])));
				return (candidate->mass());
		}
		return (0.0);
}


double TopPlusXCandidates::mass4jets() const
{
		unsigned int cnt = 0;

		ParticlePointer candidate;
    for (unsigned int index = 0; index < goodJets.size(); ++index) {
			if (goodJets.at(index) != leptonicBJet) {
				if (cnt == 0) {
					candidate = ParticlePointer(new Particle(*goodJets.at(index)));
				} else if (cnt < 3) {
					candidate = ParticlePointer(new Particle(*goodJets.at(index) + *candidate));
				}
				++cnt;
			}
		}
		if (cnt > 0)
			return (candidate->mass());
		return (0.0);
}


JetPointer TopPlusXCandidates::leadnontopjet() const {
		double hiJetPt = 0;
		unsigned int hiInd = 0;

    for (unsigned int index = 0; index < goodJets.size(); ++index) {
        if (goodJets.at(index) != leptonicBJet && goodJets.at(index)->pt() > hiJetPt) {
          hiJetPt = goodJets.at(index)->pt();
					hiInd = index;
			}
    }
		return (goodJets.at(hiInd));
}


// ST is the sum of jet and lepton pT minus the pT of the leading lepton and jet.
double TopPlusXCandidates::ST() const {
    double st = fullHT() - met->pt();

		double hiLepPt = 0;
    for (unsigned int index = 0; index < goodPFIsolatedElectrons.size(); ++index) {
         if (goodPFIsolatedElectrons.at(index)->pt() > hiLepPt)
          hiLepPt = goodPFIsolatedElectrons.at(index)->pt();
    }
    for (unsigned int index = 0; index < goodIsolatedMuons.size(); ++index) {
         if (goodIsolatedMuons.at(index)->pt() > hiLepPt)
          hiLepPt = goodIsolatedMuons.at(index)->pt();
    }
		double hiJetPt = 0;
    for (unsigned int index = 0; index < goodJets.size(); ++index) {
         if (goodJets.at(index)->pt() > hiJetPt)
          hiJetPt = goodJets.at(index)->pt();
    }
		st -= hiLepPt + hiJetPt;
    return st;
}

double TopPlusXCandidates::getTopMetAngle() const
{
    return (getLeptonicTop()->getFourVector().Theta() - met->getFourVector().Theta());
}


void TopPlusXCandidates::recoBestLeptoTop(LeptonPointer lepton,
		const JetPointer hadTopJet, const JetPointer wJet1, const JetPointer wJet2)
{
	if (goodBJets.size() < 1)
		throw ReconstructionException("Not enough jets available to reconstruct top using Chi2 method.");
	double bestLepChi2(9999999.);
	// JetCollection *topJets = &goodBJets;
	JetCollection *topJets = &goodJets;
	if (topJets->size() <= 1 && hadTopJet != 0)
		topJets = &goodJets;
	leptonFromW = lepton;
  selectedNeutrino = 0;
  currentSelectedNeutrino = 0;
  reconstructNeutrinos();

	for (unsigned short jet1Index = 0; jet1Index < topJets->size(); ++jet1Index) {
		if (hadTopJet == topJets->at(jet1Index) ||
				dJetFromWp == topJets->at(jet1Index) ||
				wJet1 == topJets->at(jet1Index) ||
				wJet2 == topJets->at(jet1Index))
			continue;
		leptonicBJet = topJets->at(jet1Index);
		leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
		leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
		leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
		leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
		selectNeutrinoSolution();
		double ltopchi2 = TopPairEventCandidate::getLeptonicChi2(currentSelectedNeutrino);
		if (ltopchi2 < bestLepChi2) {
			selectedNeutrino = currentSelectedNeutrino;
			leptonicBIndex = jet1Index;
			bestLepChi2 = ltopchi2;
		}
	}
	leptonicBJet = topJets->at(leptonicBIndex);
  leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
  leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
	leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
	leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
 	doneReconstructiontop = true;
 	doneReconstruction = true;
}


void TopPlusXCandidates::findDJet(const JetPointer topJet)
{
	dJetFromWpIndex = -1;
	// dJetFromWp = 0;  // unnecessary
	unsigned short cnt = 0;
	do {
		for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index)
		{
			if ((goodJets.at(jet1Index)->isBJet() == false || cnt >= 1) &&
					goodJets.at(jet1Index) != topJet &&
					jet1FromW != goodJets.at(jet1Index) &&
					jet2FromW != goodJets.at(jet1Index))
			{
				dJetFromWpIndex = jet1Index;
				dJetFromWp = goodJets.at(jet1Index);
				break;
			}
		}
		++cnt;
	} while (dJetFromWpIndex == -1 && cnt < 2);
	if (dJetFromWpIndex == -1)
		throw ReconstructionException("No d jet found -- all jets used for tops.");
}


void TopPlusXCandidates::recoHadronicTop(LeptonPointer lepton,
		const JetPointer lepTopJet)
{
	if (goodJets.size() < 4)
		throw ReconstructionException("Not enough jets available to reconstruct top using Chi2 method.");
	double bestHadChi2(9999999.);
	// JetCollection *topJets = &goodBJets;
	JetCollection *topJets = &goodJets;
	if (topJets->size() <= 1 && lepTopJet != 0)
		topJets = &goodJets;

	for (unsigned short hadBindex = 0; hadBindex < topJets->size(); ++hadBindex) {
		if (lepTopJet == topJets->at(hadBindex) ||
				dJetFromWp == topJets->at(hadBindex))
			continue;
		for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
			if (lepTopJet == goodJets.at(jet1Index) || 
					goodJets.at(jet1Index) == topJets->at(hadBindex) ||
					jet1Index == dJetFromWpIndex)
				continue;
			for (unsigned short jet2Index = 0; jet2Index < goodJets.size();
				++jet2Index) {
				if (jet2Index == jet1Index ||
						goodJets.at(jet2Index) == topJets->at(hadBindex) ||
						jet2Index == dJetFromWpIndex ||
						lepTopJet == goodJets.at(jet2Index))
					continue;
				hadronicBJet = topJets->at(hadBindex);
				jet1FromW = goodJets.at(jet1Index);
				jet2FromW = goodJets.at(jet2Index);
				hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
				hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
				double htopchi2 = getHadronicChi2();
				if (htopchi2 < bestHadChi2) {
					hadronicBIndex = hadBindex;
					jet1FromWIndex = jet1Index;
					jet2FromWIndex = jet2Index;
					bestHadChi2 = htopchi2;
				}
			}
		}
	}
	hadronicBJet = topJets->at(hadronicBIndex);
	jet1FromW = goodJets.at(jet1FromWIndex);
	jet2FromW = goodJets.at(jet2FromWIndex);
	hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
  hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
 	doneReconstructiontop = true;
 	doneReconstruction = true;
}


TwoNonResTops::TwoNonResTops() :
	wpGenParticlesSet(false)
{
}
		
	
TwoNonResTops::TwoNonResTops(const Event& event) :
	TopPlusXCandidates(event),
	wpGenParticlesSet(false)
{
}


TwoNonResTops::~TwoNonResTops()
{
}

double TwoNonResTops::getGlobalChi2() const {	// No global chi2 for this class
	return (0.0);
}


void TwoNonResTops::recoNonResTops(LeptonPointer lepton) {
    if (goodJets.size() < 4)
      throw ReconstructionException("Not enough jets available to reconstruct two tops");
		// const JetCollection &topJets = goodBJets;
		const JetCollection &topJets = goodJets;
    leptonFromW = lepton;
    selectedNeutrino = 0;
    currentSelectedNeutrino = 0;
    reconstructNeutrinos();
    double chosen_TopMassDifference(9999999.);
    double chosen_Chi2Total(9999999.);

		for (unsigned short hadBindex = 0; hadBindex < topJets.size(); ++hadBindex) {
					for (unsigned short lepBindex = 0; lepBindex < topJets.size(); ++lepBindex) {
							if (lepBindex == hadBindex)
									continue;
							for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
									if (goodJets.at(jet1Index) == topJets.at(lepBindex) ||
											goodJets.at(jet1Index) == topJets.at(hadBindex))
										continue;
									for (unsigned short jet2Index = 0; jet2Index < goodJets.size(); ++jet2Index) {
											if (jet2Index == jet1Index ||
													goodJets.at(jet2Index) == topJets.at(lepBindex) ||
													goodJets.at(jet2Index) == topJets.at(hadBindex))
												continue;
											hadronicBJet = topJets.at(hadBindex);
											leptonicBJet = topJets.at(lepBindex);
											jet1FromW = goodJets.at(jet1Index);
											jet2FromW = goodJets.at(jet2Index);

											leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
											leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
											hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
											leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
											leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
											hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
											fillHypotheses();
											selectNeutrinoSolution();
											double TopMassDifference = calculateTopMassDifference(currentSelectedNeutrino);
											double chi2 = getTotalChi2(currentSelectedNeutrino);
											switch (usedTTbarReconstruction) {
											case TTbarReconstructionCriterion::TopMassDifference:
												if (TopMassDifference < chosen_TopMassDifference) {
													hadronicBIndex = hadBindex;
													leptonicBIndex = lepBindex;
													jet1FromWIndex = jet1Index;
													jet2FromWIndex = jet2Index;
													chosen_TopMassDifference = TopMassDifference;
													selectedNeutrino = currentSelectedNeutrino;
												}
												break;
											case TTbarReconstructionCriterion::chi2:
												if (chi2 < chosen_Chi2Total) {
													hadronicBIndex = hadBindex;
													leptonicBIndex = lepBindex;
													jet1FromWIndex = jet1Index;
													jet2FromWIndex = jet2Index;
													chosen_Chi2Total = chi2;
													selectedNeutrino = currentSelectedNeutrino;
												}
												break;
											}
									}
							}
					}
	}
	std::sort(solutions.begin(), solutions.end(), compareSolutions);
	hadronicBJet = topJets.at(hadronicBIndex);
	leptonicBJet = topJets.at(leptonicBIndex);
	jet1FromW = goodJets.at(jet1FromWIndex);
	jet2FromW = goodJets.at(jet2FromWIndex);
	leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
	leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
	hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
	leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
	leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
	hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
	if (selectedNeutrino == 1)	// Resonance may be real of fake -- can exclude ttbar?
		ttbarResonance = ParticlePointer(new Particle(*leptonicTop1 + *hadronicTop));
	else
		ttbarResonance = ParticlePointer(new Particle(*leptonicTop2 + *hadronicTop));
	doneReconstruction = true;
}


double TwoNonResTops::getTotalChi2DJ(unsigned short neutrinoSolution) const
{
	double chi2 = TopPairEventCandidate::getTotalChi2(neutrinoSolution);
	double pt = dJetFromWp->pt();
	if (pt < jet1FromW->pt())
		chi2 += (0.5/sqrt(2.0));
	if (pt < jet2FromW->pt())
		chi2 += (0.5/sqrt(2.0));
	double deltaphi = dJetFromWp->deltaPhi(hadronicBJet);
	double deltaphi2 = dJetFromWp->deltaPhi(leptonicBJet);
	if (deltaphi2 > deltaphi)
		deltaphi = deltaphi2;
	double dphichi = deltaphi - pi;
	chi2 += (1.0/sqrt(2.0)) * (dphichi * dphichi)/(pi * pi/2.0);
	if (dJetFromWp->isBJet())
		chi2 += (0.5/sqrt(2.0));
	return (chi2);
}


void TwoNonResTops::recoNonResTopsWdjet(LeptonPointer lepton) {
	if (goodJets.size() < 5)
		throw ReconstructionException("Not enough jets available to reconstruct two tops");
	const JetCollection &topJets = goodJets;
	leptonFromW = lepton;
	selectedNeutrino = 0;
	currentSelectedNeutrino = 0;
	reconstructNeutrinos();
	double chosen_TopMassDifference(9999999.);
	double chosen_Chi2Total(9999999.);

	for (unsigned short hadBindex = 0; hadBindex < topJets.size(); ++hadBindex) {
		for (unsigned short lepBindex = 0; lepBindex < topJets.size(); ++lepBindex) {
			if (lepBindex == hadBindex)
					continue;
			for (unsigned short jet1Index = 0; jet1Index < goodJets.size(); ++jet1Index) {
				if (goodJets.at(jet1Index) == topJets.at(lepBindex) ||
						goodJets.at(jet1Index) == topJets.at(hadBindex))
					continue;
				for (unsigned short jet2Index = 0; jet2Index < goodJets.size(); ++jet2Index) {
					if (jet2Index == jet1Index ||
							goodJets.at(jet2Index) == topJets.at(lepBindex) ||
							goodJets.at(jet2Index) == topJets.at(hadBindex))
						continue;
					for (unsigned short djetIndex = 0; djetIndex < goodJets.size(); ++djetIndex) {
						if (djetIndex == jet1Index || jet2Index == djetIndex ||
								goodJets.at(djetIndex) == topJets.at(lepBindex) ||
								goodJets.at(djetIndex) == topJets.at(hadBindex))
							continue;
						hadronicBJet = topJets.at(hadBindex);
						leptonicBJet = topJets.at(lepBindex);
						jet1FromW = goodJets.at(jet1Index);
						jet2FromW = goodJets.at(jet2Index);
						dJetFromWp = goodJets.at(djetIndex);

						leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
						leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
						hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
						leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
						leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
						hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
						fillHypotheses();
						selectNeutrinoSolution();
						double TopMassDifference = calculateTopMassDifference(currentSelectedNeutrino);
						double chi2 = getTotalChi2DJ(currentSelectedNeutrino);
						switch (usedTTbarReconstruction) {
						case TTbarReconstructionCriterion::TopMassDifference:
							if (TopMassDifference < chosen_TopMassDifference) {
								hadronicBIndex = hadBindex;
								leptonicBIndex = lepBindex;
								jet1FromWIndex = jet1Index;
								jet2FromWIndex = jet2Index;
								chosen_TopMassDifference = TopMassDifference;
								selectedNeutrino = currentSelectedNeutrino;
							}
							break;
						case TTbarReconstructionCriterion::chi2:
							if (chi2 < chosen_Chi2Total) {
								hadronicBIndex = hadBindex;
								leptonicBIndex = lepBindex;
								jet1FromWIndex = jet1Index;
								jet2FromWIndex = jet2Index;
								dJetFromWpIndex = djetIndex;
								chosen_Chi2Total = chi2;
								selectedNeutrino = currentSelectedNeutrino;
							}
							break;
						}
					}
				}
			}
		}
	}
	std::sort(solutions.begin(), solutions.end(), compareSolutions);
	hadronicBJet = topJets.at(hadronicBIndex);
	leptonicBJet = topJets.at(leptonicBIndex);
	jet1FromW = goodJets.at(jet1FromWIndex);
	jet2FromW = goodJets.at(jet2FromWIndex);
	dJetFromWp = goodJets.at(dJetFromWpIndex);
	leptonicW1 = ParticlePointer(new Particle(*neutrino1 + *leptonFromW));
	leptonicW2 = ParticlePointer(new Particle(*neutrino2 + *leptonFromW));
	hadronicW = ParticlePointer(new Particle(*jet1FromW + *jet2FromW));
	leptonicTop1 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW1));
	leptonicTop2 = ParticlePointer(new Particle(*leptonicBJet + *leptonicW2));
	hadronicTop = ParticlePointer(new Particle(*hadronicBJet + *hadronicW));
	if (selectedNeutrino == 1)	// Resonance may be real of fake -- can exclude ttbar?
		ttbarResonance = ParticlePointer(new Particle(*leptonicTop1 + *hadronicTop));
	else
		ttbarResonance = ParticlePointer(new Particle(*leptonicTop2 + *hadronicTop));
	doneReconstruction = true;
}


const JetPointer TwoNonResTops::getLtNonTopJet() const {
    throwExpeptionIfNotReconstructed("TwoNonResTops::getLtNonTopJet");
    if (goodJets.size() < 5)
      throw ReconstructionException("Not enough jets available to reconstruct Wprime");
		int ltJetInd = -1, bJetInd = -1;
    unsigned int index = 0;
    while (index < goodJets.size() && ltJetInd < 0) {
			JetPointer currJet = goodJets.at(index);
			if (currJet->isBJet() == false && currJet != leptonicBJet &&
					currJet != hadronicBJet && currJet != jet1FromW && currJet != jet2FromW) {
				ltJetInd = index;
			} else if (bJetInd < 0 && currJet != leptonicBJet && currJet != hadronicBJet &&
					currJet != jet1FromW && currJet != jet2FromW) {
				bJetInd = index;
			}
			++index; 
		}
		if (ltJetInd < 0 && bJetInd >= 0)	// Take a B jet if no others found.
			ltJetInd = bJetInd;
		if (ltJetInd >= 0)
    	return (goodJets.at(ltJetInd));
		throw ReconstructionException(TString("getLtNonTopJet:Fifth jet missing"));
}


const ParticlePointer TwoNonResTops::getWprime() const {
    throwExpeptionIfNotReconstructed("TwoNonResTops::getWprime");
		ParticlePointer ltJet = getLtNonTopJet();
		ParticlePointer hiPtTop = getLeptonicTop();

		// Generator study indicates lower Pt top is usually spectator
		if (hiPtTop->pt() < hadronicTop->pt())
			hiPtTop = hadronicTop;
		return (ParticlePointer(new Particle(*hiPtTop + *ltJet)));
}


const ParticlePointer TwoNonResTops::getNegWprime() const {
    throwExpeptionIfNotReconstructed("TwoNonResTops::getNegWprime");
		ParticlePointer ltJet = getLtNonTopJet();
		ParticlePointer negTop = hadronicTop;
		if (leptonFromW->charge() < 0)
			negTop = getLeptonicTop();

		return (ParticlePointer(new Particle(*negTop + *ltJet)));
		// return (ParticlePointer(new Particle(*negTop + *dJetFromWp)));
}


// Deliberately reco the incorrect W' to see if ttbar still looks the same
const ParticlePointer TwoNonResTops::getNegWprimeBad() const {
    throwExpeptionIfNotReconstructed("TwoNonResTops::getNegWprimeBad");
		ParticlePointer ltJet = getLtNonTopJet();
		ParticlePointer negTop = hadronicTop;
		if (leptonFromW->charge() > 0)
			negTop = getLeptonicTop(); // Use wrong top

		return (ParticlePointer(new Particle(*negTop + *ltJet)));
		// return (ParticlePointer(new Particle(*negTop + *dJetFromWp)));
}


const ParticlePointer TwoNonResTops::getBtBWprime() const {
    throwExpeptionIfNotReconstructed("TwoNonResTops::getBtBWprim");
		ParticlePointer ltJet = getLtNonTopJet();
		ParticlePointer hiAnglTop = getLeptonicTop();
		if (getHadTopDAngle() > hiAnglTop->angle(ltJet))
			hiAnglTop = hadronicTop;
		return (ParticlePointer(new Particle(*hiAnglTop + *ltJet)));
}


double TwoNonResTops::getHadTopDAngle() const
{
 ParticlePointer jet = getLtNonTopJet();
 return (hadronicTop->angle(jet));
}


void TwoNonResTops::recoWpFrom1Top(LeptonPointer lepton) {
	if (goodJets.size() < 5)
		throw
			ReconstructionException("Not enough jets available to reconstruct W'");
	ParticlePointer negTop, posTop;
	if (lepton->charge() < 0) {
		// findDJet();
		recoBestLeptoTop(lepton);
		// findDJet(leptonicBJet);
		recoHadronicTop(lepton);
		// recoHadronicTop(lepton, leptonicBJet);
		negTop = getLeptonicTop();
		posTop = hadronicTop;
	} else {
		// findDJet();
		recoHadronicTop(lepton);
		// findDJet(hadronicBJet);
		recoBestLeptoTop(lepton);
		// recoBestLeptoTop(lepton, hadronicBJet, jet1FromW, jet2FromW);
		negTop = hadronicTop;
		posTop = getLeptonicTop();
	}
	JetPointer ltJet = getLtNonTopJet();
	wpFrom1Top = ParticlePointer(new Particle(*negTop + *ltJet));
	wpBadFrom1Top = ParticlePointer(new Particle(*posTop + *ltJet));
	topFromWp = negTop;
	dJetFromWp = ltJet;
	// wpFrom1Top = ParticlePointer(new Particle(*negTop + *dJetFromWp));
	// wpBadFrom1Top = ParticlePointer(new Particle(*posTop + *dJetFromWp));
}


void TwoNonResTops::getTrueWpChain()
{
	const MCParticleCollection &mcColl = genParticles;
	const int nullInd = -99;
	int wpInd = nullInd, wpTopInd = nullInd, specTopInd = nullInd, dQInd = nullInd;
	int wpTopBInd = nullInd, specTopBInd = nullInd;
	if (mcColl.size() > 0) {
		for (unsigned int i = 0; i < mcColl.size(); ++i) {
			if (fabs(mcColl.at(i)->pdgId()) == WPRIME_ID) {
				if (wpInd == nullInd)
					wpInd = i;
				else if (mcColl.at(i)->motherIndex() != wpInd)
					cout << "Stray W': " << i << " mother " << mcColl.at(i)->motherIndex() << endl;
			}
		}
		for (unsigned int i = 0; i < mcColl.size(); ++i) {
			if (fabs(mcColl.at(i)->pdgId()) == TOP_ID) {
				if (mcColl.at(i)->motherIndex() == wpInd) {
					if (wpTopInd == nullInd)
						wpTopInd = i;
					else cout << "Additional W' top: " << i << endl;
				} else if (specTopInd == nullInd)
						specTopInd = i;
					else cout << "Additional spec top: " << i << endl;
			} else if (fabs(mcColl.at(i)->pdgId()) == DQ_ID) {
				if (mcColl.at(i)->motherIndex() == wpInd) {
					if (dQInd == nullInd)
						dQInd = i;
					else cout << "Additional d daughter: " << i << endl;
				}
			}
		}
		if (wpTopInd >= 0 && specTopInd >= 0) {
			int wptprevPdgId = 0;
			int sptprevPdgId = 0;
			for (unsigned int i = 0; i < mcColl.size(); ++i) {
				int pdgId = fabs(mcColl.at(i)->pdgId());
				if (mcColl.at(i)->motherIndex() == wpTopInd ||
						mcColl.at(i)->motherIndex() == wpInd) {
					if (mcColl.at(i)->isQuark() && pdgId % 2 == 1) {
						if (wpTopBInd == nullInd || pdgId > wptprevPdgId) {
							wpTopBInd = i;
							wptprevPdgId = pdgId;
						} else if (pdgId == BQ_ID)
							cout << "Additional W' top b quark: " << i << " mother " << mcColl.at(i)->motherIndex() << endl;
					}
				} else if (mcColl.at(i)->motherIndex() == specTopInd ||
						mcColl.at(i)->motherIndex() == mcColl.at(specTopInd)->motherIndex()) {
					if (mcColl.at(i)->isQuark() && pdgId % 2 == 1) {
						if (specTopBInd == nullInd || pdgId > sptprevPdgId) {
							specTopBInd = i;
							sptprevPdgId = pdgId;
						} else if (pdgId == BQ_ID)
							cout << "Additional spec top b quark: " << i << " mother " << mcColl.at(i)->motherIndex() << endl;
					}
				} else if (pdgId == BQ_ID && fabs(mcColl.at(mcColl.at(i)->motherIndex())->pdgId()) != BQ_ID) {
					// cout << "Additional b quark: " << i << " mother " << mcColl.at(i)->motherIndex();
					// cout << " ID " << mcColl.at(mcColl.at(i)->motherIndex())->pdgId() << endl;
				}
			}
		}
		if (wpInd >= 0 && wpTopInd >= 0 && specTopInd >= 0 && dQInd >= 0 &&
				wpTopBInd >= 0 && specTopBInd >= 0) {
			wpTru = mcColl.at(wpInd);
			wpTopTru = mcColl.at(wpTopInd);
			specTopTru = mcColl.at(specTopInd);
			wpTopBTru = mcColl.at(wpTopBInd);
			specTopBTru = mcColl.at(specTopBInd);
			wpDDauTru = mcColl.at(dQInd);
			wpGenParticlesSet = true;
		} else {
			// cout << "#11 : " << mcColl.at(11)->pdgId() << " mother " << mcColl.at(11)->motherIndex();
			// cout << " ID " << mcColl.at(mcColl.at(11)->motherIndex())->pdgId() << endl;
			cout << "wp wpTop specTop dq wpTopB specTopB " << wpInd << " " << wpTopInd << " " << specTopInd << " " << dQInd;
			cout << " " << wpTopBInd << " " << specTopBInd << endl;
		}
	}
}


MCParticlePointer TwoNonResTops::getWpTru() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (wpTru);
}

MCParticlePointer TwoNonResTops::getWpTopTru() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (wpTopTru);
}

MCParticlePointer TwoNonResTops::getSpecTop() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (specTopTru);
}

MCParticlePointer TwoNonResTops::getWpTopBTru() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (wpTopBTru);
}

MCParticlePointer TwoNonResTops::getSpecTopB() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (specTopBTru);
}

MCParticlePointer TwoNonResTops::getWpDDauTru() const
{
	if (wpGenParticlesSet == false)
		throw ReconstructionException("Need to run TwoNonResTops::getTrueWpChain first");
	return (wpDDauTru);
}


TopNoMassConstraint::TopNoMassConstraint()
{
}
		
	
TopNoMassConstraint::TopNoMassConstraint(const Event& event) :
	ToplikeCandidate(event)
{
}


TopNoMassConstraint::~TopNoMassConstraint()
{
}


double TopNoMassConstraint::getLoneHadChi2() const
{
    double ptRatioDifference = TMath::Power(PtRatio() - matched_ptratio, 2);
    double ptRatioError = 2 * matched_ptratio_sigma * matched_ptratio_sigma;
    double ptRatioTerm = ptRatioDifference / ptRatioError;

    double WmassDifference = TMath::Power(hadronicW->mass() - pdgWmass, 2);
    double WmassError = 2 * wMassWidth * wMassWidth;
    double WmassTerm = WmassDifference / WmassError;
    return 1.0 / sqrt(2.0) * (WmassTerm + ptRatioTerm);
}

double TopNoMassConstraint::getLepBAngle() const
{
 return (leptonicBJet->angle(leptonFromW));
}


LeptoTopNoMassConstraint::LeptoTopNoMassConstraint()
{
}
		
	
LeptoTopNoMassConstraint::LeptoTopNoMassConstraint(const Event& event) :
	TopNoMassConstraint(event)
{
}


LeptoTopNoMassConstraint::~LeptoTopNoMassConstraint()
{
}


/*
Tried for neutrino selection, but getGlobalChi2 still uses hadronic top quantities.
double LeptoTopNoMassConstraint::getTotalChi2(unsigned short int neutrinoSolution) const {
    return TopPairEventCandidate::getLeptonicChi2(neutrinoSolution) + TopPairEventCandidate::getGlobalChi2(neutrinoSolution);
}
*/

void LeptoTopNoMassConstraint::inspectReconstructedEvent() const {
    cout << "run " << runNumber << ", event " << eventNumber << endl;
    cout << "leptonic b jet, goodJet index " << leptonicBIndex << endl;
    EventContentPrinter::printJet(leptonicBJet);

    cout << "lepton from W" << endl;
    // EventContentPrinter::printElectron(electronFromW);

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

    cout << endl;
}


}
