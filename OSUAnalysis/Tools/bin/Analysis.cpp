/*
 * Analysis.cpp
 *
 *  Created on: 12 Jul 2010
 *      Author: kreczko
 */

#include "Analysis.h"
#include "TROOT.h"
#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>
#include "../interface/EventCounter.h"
#include <cmath>
#include <math.h>
#include "../interface/Printers/EventTablePrinter.h"
#include "TTree.h"

using namespace BAT;
using namespace std;
float Analysis::luminosity = 36.145;

void Analysis::analyze(int const Section, TString const ListName) {
  if (Section == -1) {
    outputFile = new TFile("micro.root","RECREATE");
  } else {
    if (ListName == "") {
      outputFile = new TFile(TString::Format("micro_%i.root", Section),"RECREATE");
    } else {
      outputFile = new TFile(ListName + TString::Format("_micro_%i.root", Section),"RECREATE");
    }
  }
    microTuple = new TTree("micro","micro");
    createHistograms();
    initMicroNtuple();
    cout << "detected samples:" << endl;
    for (unsigned int sample = 0; sample < DataType::NUMBER_OF_DATA_TYPES; ++sample) {
        if (eventReader->getSeenDatatypes()[sample])
            cout << DataType::names[sample] << endl;
    }
    while (eventReader->hasNextEvent() ) {
        if( !initiateEvent() ) continue;
        printNumberOfProccessedEventsEvery(100000);
        inspectEvents();
        //        doEcalSpikeAnalysis();
        //                doSynchExercise();
//        if (currentEvent.GoodElectrons().size() > 0) {
//        }
//        doPileUpStudy();
				// cout << "Cutflow step\n";
        doToplikeCutFlow(singleCuts, singleCutsPerFile,
				    cutflow, cutflowPerFile, cutflowPerSample,
						toplikeElectronSelection, toplikeElectronSelSize);
        doToplikeCutFlow(musingleCuts, musingleCutsPerFile,
				    mucutflow, mucutflowPerFile, mucutflowPerSample,
						toplikeMuonSelection, toplikeMuonSelSize);
        // doTTbarCutFlow();
       // doDiElectronAnalysis();
				// cout << "analysis step\n";
				doTransverseMassPlots();
        doTTBarAnalysis();
        doMicroNtuple();
        // doNotePlots();
        // doQCDStudy();
        // doJetAnalysis();
        // if (currentEvent.getDataType() == DataType::ttbar)
        	// doMCttbarReconstruction();

//        if(currentEvent.getDataType() == DataType::DATA)
//            eventCheck[currentEvent.runnumber()].push_back(currentEvent.eventnumber());
//        checkForBrokenEvents();
    }
//    checkForDuplicatedEvents();
    printInterestingEvents();
    printSummary();
}

void Analysis::printNumberOfProccessedEventsEvery(unsigned long printEvery) {
    unsigned long eventIndex = eventReader->getNumberOfProccessedEvents();
    if (eventIndex % printEvery == 0 || eventIndex == 1) {
        cout << "Analysing event no " << eventIndex << ", sample: " << DataType::names[currentEvent.getDataType()]
                << endl;
        cout << "File: " << eventReader->getCurrentFile() << endl;
    } else {
    	/*
    	static TString prevFile;
    	const char *currFile = eventReader->getCurrentFile();
    	if (prevFile != currFile) {
        cout << "File: " << currFile << endl;
        prevFile = currFile;
    	}
    	*/
    }

}

extern bool isGood(int,int); 
extern bool isGoodMay10ReRecoV3(int,int);
extern bool isGoodReReco5AugV3(int,int);

bool Analysis::initiateEvent() {
    currentEvent = eventReader->getNextEvent();
    string filename( eventReader->getCurrentFile() );

    if( currentEvent.isRealData() && (
		       !isGood(currentEvent.runnumber(), currentEvent.lumiblock())
      	      /* (filename.find("PromptReco") != string::npos && !isGood             (currentEvent.runnumber(), currentEvent.lumiblock())) ||
 (filename.find("May10ReReco")!= string::npos && !isGoodMay10ReRecoV3(currentEvent.runnumber(), currentEvent.lumiblock())) ||
 (filename.find("ReReco5Aug") != string::npos && !isGoodReReco5AugV3 (currentEvent.runnumber(), currentEvent.lumiblock())) ||
 (filename.find("PromptReco")==string::npos && filename.find("May10ReReco")==string::npos && filename.find("ReReco5Aug")==string::npos )*/ // Unknown
                                     )
    ) return false;

    ttbarCandidate = TopPairEventCandidate(currentEvent);
    weight = weights.getWeight(currentEvent.getDataType());
    if(!currentEvent.isRealData() //&&
//			currentEvent.getDataType() != DataType::WprimeTToTTD_M600 &&
//			currentEvent.getDataType() != DataType::WprimeTToTTD_M1000
    ) {
			// W' samples lack pile-up.
        weight *= weights.reweightPileUp(currentEvent.numberOfGeneratedPileUpVertices());
	PileUpWeight= weights.reweightPileUp(currentEvent.numberOfGeneratedPileUpVertices()); 
	//	std::cout<<PileUpWeight<<" "<<weight<<std::endl;
    }
    else if (currentEvent.isRealData())  PileUpWeight=1;
    tPlusXCandidates = ToplikeCandidate(currentEvent);
    twoTopsNonRes = TwoNonResTops(currentEvent);
    goodTopBadTop = TwoNonResTops(currentEvent);
    // loneTopsNoMassConstr = TopNoMassConstraint(currentEvent);
    // leptoTopNoMassConstr =  LeptoTopNoMassConstraint(currentEvent);
    histMan.setCurrentDataType(ttbarCandidate.getDataType());
    histMan.setCurrentJetBin(currentEvent.GoodJets().size());
    histMan.setCurrentBJetBin(currentEvent.GoodBJets().size());
    return true;
}

void Analysis::inspectEvents() {
    std::vector<InterestingEvent> eventsToInspect;

    for (unsigned int index = 0; index < eventsToInspect.size(); ++index) {
        if ((ttbarCandidate.runnumber() == eventsToInspect.at(index).runNumber && ttbarCandidate.eventnumber()
                == eventsToInspect.at(index).eventNumber)) {
            cout << "file: " << eventReader->getCurrentFile() << endl;
            ttbarCandidate.inspect();
        }
    }

}

void Analysis::doSynchExercise() {
    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::ConversionFinder)) {
        cout << ttbarCandidate.runnumber() << ":" << ttbarCandidate.eventnumber() << ":" << endl;//electron->et() << endl;
        if (ttbarCandidate.eventnumber() == 450622) {
            ttbarCandidate.inspect();
        }
    }
}


void Analysis::doToplikeCutFlow(cutarray &singleCuts, cutmap &singleCutsPerFile,
		cutarray &cutflow, cutmap &cutflowPerFile, BAT::Counter &cutflowPerSample,
		const ToplikeSelectionSteps::Step cutList[], unsigned int cutListSiz)
{
	for (unsigned int cut = 0; cut < cutListSiz; ++cut) {
		if (tPlusXCandidates.passesSelectionStep(cutList[cut])) {
				++singleCuts[cut];
				singleCutsPerFile[eventReader->getCurrentFile()][cut]++;
		}
		if (tPlusXCandidates.passesSelectionStepUpTo(cut, cutList)) {
				cutflow[cut] += 1;
				cutflowPerFile[eventReader->getCurrentFile()][cut]++;
				unsigned int njet = tPlusXCandidates.GoodJets().size();
				if (njet >= JetBin::NUMBER_OF_JET_BINS)
						njet = JetBin::NUMBER_OF_JET_BINS - 1;
				cutflowPerSample.increase(tPlusXCandidates.getDataType(), cut, njet, weight);
		}
	}
}


/*
void Analysis::doTTbarCutFlow() {
    for (unsigned int cut = 0; cut < TTbarEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
        if (ttbarCandidate.passesSelectionStep((TTbarEPlusJetsSelection::Step) cut)) {
            ++singleCuts[cut];
            singleCutsPerFile[eventReader->getCurrentFile()][cut]++;
        }

        if (ttbarCandidate.passesSelectionStepUpTo((TTbarEPlusJetsSelection::Step) cut)) {
            cutflow[cut] += 1;
            cutflowPerFile[eventReader->getCurrentFile()][cut]++;
            unsigned int njet = ttbarCandidate.GoodJets().size();
            if (njet >= JetBin::NUMBER_OF_JET_BINS)
                njet = JetBin::NUMBER_OF_JET_BINS - 1;
            cutflowPerSample.increase(ttbarCandidate.getDataType(), cut, njet, weight);
        }
    }
}
*/

void Analysis::doJetAnalysis() {
    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex)) {
        for(unsigned short jetIndex = 0; jetIndex < ttbarCandidate.Jets().size(); ++jetIndex)
            histMan.H1D_BJetBinned("AllJetMass")->Fill(ttbarCandidate.Jets().at(jetIndex)->mass());

        for (unsigned short jetIndex = 0; jetIndex < ttbarCandidate.GoodJets().size(); ++jetIndex) {
            double jetMass = ttbarCandidate.Jets().at(jetIndex)->mass();
            histMan.H1D_BJetBinned("AllGoodJetMass")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastOneGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastOneJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastTwoGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastTwoJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastThreeGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastThreeJets")->Fill(jetMass, weight);

            if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::AtLeastFourGoodJets))
                histMan.H1D_BJetBinned("GoodJetMass_atLeastFourJets")->Fill(jetMass, weight);
        }
    }
}

void Analysis::doDiElectronAnalysis() {
    ElectronCollection electrons = currentEvent.GoodElectrons();
    if (electrons.size() == 2) {
        ElectronPointer leadingElectron = electrons.front();
        ElectronPointer secondElectron = electrons.at(1);
        histMan.H1D_JetBinned("diElectronMass")->Fill(leadingElectron->invariantMass(secondElectron), weight);
    }
}

void Analysis::fillMTplots(const double mt, const  float charge,
	const std::string histPrefix, const TString &leptSuffix, const unsigned int step)
{
	std::string bothname = histPrefix, plusname = histPrefix, minusname = histPrefix;
	bothname += "bothchg";
	plusname += "plus";
	minusname += "minus";
	std::string bothlepname = histPrefix, pluslepname = histPrefix,
		minuslepname = histPrefix;
	bothlepname += "bothchglep";
	pluslepname += "pluslep";
	minuslepname += "minuslep";

	bothname += leptSuffix;
	plusname += leptSuffix;
	minusname += leptSuffix;
	TString suffix = "4j";
	if (step == 1)
		suffix = "5j";
	bothname += suffix;
	plusname += suffix;
	minusname += suffix;
	bothlepname += suffix;
	pluslepname += suffix;
	minuslepname += suffix;
	histMan.H1D_BJetBinned(bothname)->Fill(mt, weight);
	histMan.H1D_BJetBinned(bothlepname)->Fill(mt, weight);
	if (charge > 0) {
		histMan.H1D_BJetBinned(plusname)->Fill(mt, weight);
		histMan.H1D_BJetBinned(pluslepname)->Fill(mt, weight);
	} else {
		histMan.H1D_BJetBinned(minusname)->Fill(mt, weight);
		histMan.H1D_BJetBinned(minuslepname)->Fill(mt, weight);
	}
} 

void Analysis::doTransverseMassPlots() {
	// Use -2 to get back to 4 jet step.
	bool elecEvt = tPlusXCandidates.passesSelectionStepUpTo(toplikeElectronSelSize - 2,
		toplikeElectronSelection);
	bool muEvt = tPlusXCandidates.passesSelectionStepUpTo(toplikeMuonSelSize - 2,
		toplikeMuonSelection);
	if (elecEvt || muEvt) {
		ElectronPointer electron;
		MuonPointer muon;
		if (elecEvt) {
			if(Event::usePFIsolation) {
					electron = tPlusXCandidates.GoodPFIsolatedElectrons().front();
			} else {
					electron = tPlusXCandidates.GoodIsolatedElectrons().front();
			}
		} else muon = tPlusXCandidates.GoodIsolatedMuons().front();
		for (unsigned int step = 0; step < 2; ++step) {
			if (step == 0 || tPlusXCandidates.hasAtLeastFiveGoodJets()) {
				TString leptSuffix;
				double mt4j, mt2j;
				float charge;
				if (elecEvt) {
					mt4j = tPlusXCandidates.transMass4j(electron);
					mt2j = tPlusXCandidates.transMass2j(electron);
					leptSuffix = "elec";
					charge = electron->charge();
				} else {
					mt4j = tPlusXCandidates.transMass4j(muon);
					mt2j = tPlusXCandidates.transMass2j(muon);
					leptSuffix = "mu";
					charge = muon->charge();
				}
				fillMTplots(mt4j, charge, "mt4j", leptSuffix, step);
				fillMTplots(mt2j, charge, "mt2j", leptSuffix, step);
			}
		}
	}
}


void Analysis::doTTBarAnalysis() {
    if (ttbarCandidate.passesNMinus1(TTbarEPlusJetsSelection::AtLeastFourGoodJets)) {
        histMan.H1D("numberOfJets")->Fill(ttbarCandidate.GoodJets().size());
    }

    if (ttbarCandidate.passesSelectionStepUpTo(TTbarEPlusJetsSelection::GoodPrimaryvertex)
            && ttbarCandidate.hasNoIsolatedMuon() && ttbarCandidate.isNotAZBosonEvent()
            && ttbarCandidate.hasAtLeastFourGoodJets()) {
        const ElectronCollection coll = ttbarCandidate.Electrons();
        for (unsigned int index = 0; index < coll.size(); ++index) {
            const ElectronPointer electron = coll.at(index);
            bool passesEt = electron->et() > 30;
            bool passesEta = fabs(electron->superClusterEta()) < 2.5
                    && electron->isInCrack() == false;
            bool passesID = electron->VBTF_W70_ElectronID();
            //            bool isNotEcalSpike = electron->isEcalSpike() == false;
            bool noConversion = electron->isFromConversion() == false;
            if (passesEt && passesEta && passesID && noConversion) {
                histMan.H1D("electronD0")->Fill(electron->d0(), weight);
            }
        }
    }
    /*
    if (tPlusXCandidates.passesFullEPlusJetSelection()) {
        if(Event::usePFIsolation) {
            // tPlusXCandidates.recoBestLeptoTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
            tPlusXCandidates.recoBestSingleTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
            twoTopsNonRes.recoNonResTops(ttbarCandidate.GoodPFIsolatedElectrons().front());
        } else {
            // tPlusXCandidates.recoBestLeptoTop(ttbarCandidate.GoodIsolatedElectrons().front());
            tPlusXCandidates.recoBestSingleTop(ttbarCandidate.GoodIsolatedElectrons().front());
            twoTopsNonRes.recoNonResTops(ttbarCandidate.GoodIsolatedElectrons().front());
        }
    }
    */
    
		ElectronCollection electrons = twoTopsNonRes.GoodElectrons();
		for (unsigned int index = 0; index < electrons.size(); ++index) {
			const ElectronPointer electron = electrons.at(index);
			histMan.H1D("eleciso")->Fill(electron->pfIsolation(), weight);
		}
		MuonCollection muons = twoTopsNonRes.GoodMuons();
		for (unsigned int index = 0; index < muons.size(); ++index) {
			const MuonPointer muon = muons.at(index);
			histMan.H1D("muiso")->Fill(muon->pfIsolation(), weight);
		}
		electrons = twoTopsNonRes.Electrons();
		for (unsigned int index = 0; index < electrons.size(); ++index) {
			const ElectronPointer electron = electrons.at(index);
			histMan.H1D("electron_pt")->Fill(electron->pt(), weight);
			if (electron->pt() > 45) {
				bool zdok = fabs(electron->ZDistanceToPrimaryVertex()) < 1;
				bool d0ok = fabs(electron->d0()) < 0.02;
				int val = 2;
				if (zdok)
					--val;
				if (d0ok)
					--val;
				histMan.H1D("elvertdist")->Fill(val, weight);
			}
		}
		muons = twoTopsNonRes.Muons();
		for (unsigned int index = 0; index < muons.size(); ++index) {
			const MuonPointer muon = muons.at(index);
			histMan.H1D("muon_pt")->Fill(muon->pt(), weight);
			if (muon->pt() > 35) {
				bool zdok = fabs(muon->ZDistanceToPrimaryVertex()) < 1.0;
				bool d0ok = fabs(muon->d0()) < 0.02;
				int val = 2;
				if (zdok)
					--val;
				if (d0ok)
					--val;
				histMan.H1D("muvertdist")->Fill(val, weight);
			}
		}
		// Step 6/4 is HT>700
		if (twoTopsNonRes.passesSelectionStepUpTo(6, toplikeElectronSelection) ||
				twoTopsNonRes.passesSelectionStepUpTo(4,
				toplikeMuonSelection))
		{
			int numjets = twoTopsNonRes.GoodJets().size();
			histMan.H1D("jetmultip")->Fill(numjets, weight);
			if (numjets > 0) {
				histMan.H1D_BJetBinned("leadingJetPt")->Fill(
					twoTopsNonRes.GoodJets().front()->pt(), weight);
				if (numjets > 1) {
					histMan.H1D_BJetBinned("nextLeadingJetPt")->Fill(
						twoTopsNonRes.GoodJets().at(1)->pt(), weight);
					for(unsigned int jetIndex = 2; jetIndex <
						twoTopsNonRes.GoodJets().size(); jetIndex++)
						histMan.H1D_BJetBinned("trailingJetPt")->Fill(
							twoTopsNonRes.GoodJets().at(jetIndex)->pt(), weight);
				}
			}
		}
		if (twoTopsNonRes.passesHTCut() && twoTopsNonRes.hasAtLeastOneGoodJet() &&
				twoTopsNonRes.hasAtLeastTwoGoodJets() &&
				twoTopsNonRes.hasAtLeastFiveGoodJets())
		{
			ElectronCollection electrons = twoTopsNonRes.LoosePFElectrons();
			// cout << " loose " << electrons.size();
			for (unsigned int index = 0; index < electrons.size(); ++index) {
				const ElectronPointer electron = electrons.at(index);
				histMan.H1D("looseeleciso")->Fill(electron->pfIsolation(), weight);
			}
			/*
			MuonCollection muons = twoTopsNonRes.GoodMuons();
			for (unsigned int index = 0; index < muons.size(); ++index) {
				const MuonPointer muon = muons.at(index);
				histMan.H1D("muiso")->Fill(muon->pfIsolation(), weight);
			}
				*/
		}
		bool elecEvt = twoTopsNonRes.passesFullEPlusJetSelection();
		bool muEvt =  twoTopsNonRes.passesFullMuPlusJetSelection();
    if (elecEvt || muEvt) {
				// twoTopsNonRes.getTrueWpChain();
				LeptonBin::value lepType = LeptonBin::Electron;
				if (muEvt)
					lepType = LeptonBin::Muon;
				histMan.setCurrentLepBin(lepType);
				LeptonPointer lepton;
				if (twoTopsNonRes.GoodPFIsolatedElectrons().size() >= 1) {
					lepton = twoTopsNonRes.GoodPFIsolatedElectrons().front();
					histMan.H1D("electron_et")->Fill(lepton->et(), weight);
					if (lepton->charge() > 0) {
						histMan.H1D("positron_pt")->Fill(lepton->pt(), weight);
					} else {
						histMan.H1D("electminus_pt")->Fill(lepton->pt(), weight);
					}
				}
				if (twoTopsNonRes.GoodIsolatedMuons().size() >= 1) {
					lepton = twoTopsNonRes.GoodIsolatedMuons().front();
					histMan.H1D("muon_et")->Fill(lepton->et(), weight);
				}
				goodTopBadTop.recoWpFrom1Top(lepton);
				twoTopsNonRes.recoNonResTops(lepton);
				tPlusXCandidates.recoBestSingleTop(lepton); // For microntuple
				histMan.H1D_LepBinned("transMasslm1b")->Fill(twoTopsNonRes.transMasslm1b(lepton)->mass(), weight);
				// histMan.H2D_BJetBinned("mass_vs_transMass")->Fill(twoTopsNonRes.massl1b(lepton)->mass(),
					// twoTopsNonRes.transMasslm1b(lepton)->mass(), weight);
				histMan.H1D_LepBinned("deltaPhiTops")->Fill(twoTopsNonRes.deltaPhiTops(lepton), weight);
				histMan.H1D_LepBinned("mLeptonicTop")->Fill(twoTopsNonRes.getLeptonicTop()->mass(), weight);
				histMan.H1D_LepBinned("mHadronicTop")->Fill(twoTopsNonRes.getHadronicTop()->mass(), weight);
				histMan.H1D_LepBinned("mLepTopLone")->Fill(goodTopBadTop.getLeptonicTop()->mass(), weight);
				histMan.H1D_LepBinned("mHadTopLone")->Fill(goodTopBadTop.getHadronicTop()->mass(), weight);
				histMan.H1D_LepBinned("numGenJets")->Fill(twoTopsNonRes.GenJets().size(),
					weight);
				histMan.H1D_LepBinned("missGenJets")->Fill(twoTopsNonRes.missGenJetsForRecoJets(),
					weight);
				histMan.H1D_LepBinned("missRecoJets")->Fill(twoTopsNonRes.missRecoJetsForGenJets(),
					weight);
				int numDiffJets = 
					twoTopsNonRes.GoodJets().size() - twoTopsNonRes.GenJets().size();
				histMan.H1D_LepBinned("numDiffJets")->Fill(numDiffJets, weight);
        histMan.H1D_LepBinned("MET_ltop")->Fill(twoTopsNonRes.MET()->et(),
        	weight);
        histMan.H1D_LepBinned("METpt_ltop")->Fill(twoTopsNonRes.MET()->pt(),
        	weight);
				if (lepton->charge() > 0) {
					histMan.H1D_LepBinned("mLepTopLonepos")->Fill(goodTopBadTop.getLeptonicTop()->mass(), weight);
					histMan.H1D_LepBinned("mHadTopLonepos")->Fill(goodTopBadTop.getHadronicTop()->mass(), weight);
					histMan.H1D_BJetBinned("wpmassplus")->Fill(twoTopsNonRes.getWprime()->mass(), weight);
					histMan.H1D_LepBinned("wpnegmasspos")->Fill(twoTopsNonRes.getNegWprime()->mass(), weight);
					histMan.H1D_LepBinned("wp1topmasspos")->Fill(goodTopBadTop.wpFrom1Top->mass(), weight);
					histMan.H1D_LepBinned("wp1topmassbadpos")->Fill(goodTopBadTop.wpBadFrom1Top->mass(), weight);
					histMan.H1D_LepBinned("deltaphilepjetpos")->Fill(twoTopsNonRes.deltaPhiLepJet(lepton), weight);
					histMan.H1D_LepBinned("deltaphilepdjetpos")->Fill(twoTopsNonRes.deltaPhiLepDJet(lepton), weight);
				} else {
					histMan.H1D_LepBinned("mLepTopLoneneg")->Fill(goodTopBadTop.getLeptonicTop()->mass(), weight);
					histMan.H1D_LepBinned("mHadTopLoneneg")->Fill(goodTopBadTop.getHadronicTop()->mass(), weight);
					histMan.H1D_BJetBinned("wpmassminus")->Fill(twoTopsNonRes.getWprime()->mass(), weight);
					histMan.H1D_LepBinned("wpnegmassneg")->Fill(twoTopsNonRes.getNegWprime()->mass(), weight);
					histMan.H1D_LepBinned("wp1topmassneg")->Fill(goodTopBadTop.wpFrom1Top->mass(), weight);
					histMan.H1D_LepBinned("wp1topmassbadneg")->Fill(goodTopBadTop.wpBadFrom1Top->mass(), weight);
					histMan.H1D_LepBinned("deltaphilepjetneg")->Fill(twoTopsNonRes.deltaPhiLepJet(lepton), weight);
					histMan.H1D_LepBinned("deltaphilepdjetneg")->Fill(twoTopsNonRes.deltaPhiLepDJet(lepton), weight);
				}

				/*
				histMan.H1D_LepBinned("wpTruPt")->Fill(twoTopsNonRes.getWpTru()->pt(), weight);
				histMan.H1D_LepBinned("wpTopTruPt")->Fill(twoTopsNonRes.getWpTopTru()->pt(), weight);
				int signVal = 2;
				if (twoTopsNonRes.getWpTopTru()->pdgId() < 0)
					signVal = 1;
				histMan.H1D_LepBinned("wpchgasymm")->Fill(signVal, weight);
				if (twoTopsNonRes.getSpecTop()->pdgId() > 0)
					signVal = 5;
				else signVal = 4;
				histMan.H1D_LepBinned("wpchgasymm")->Fill(signVal, weight);

				histMan.H1D_LepBinned("specTopTruPt")->Fill(twoTopsNonRes.getSpecTop()->pt(), weight);
				histMan.H1D_LepBinned("wpTopBTruPt")->Fill(twoTopsNonRes.getWpTopBTru()->pt(), weight);
				histMan.H1D_LepBinned("specTopBTruPt")->Fill(twoTopsNonRes.getSpecTopB()->pt(), weight);
				histMan.H1D_LepBinned("wpDDauTruPt")->Fill(twoTopsNonRes.getWpDDauTru()->pt(), weight);
				histMan.H1D_LepBinned("angleWpTopDTru")->Fill(twoTopsNonRes.getWpTopTru()->deltaPhi(twoTopsNonRes.getWpDDauTru()), weight);
				histMan.H1D_LepBinned("angleSpecTopDTru")->Fill(twoTopsNonRes.getSpecTop()->deltaPhi(twoTopsNonRes.getWpDDauTru()), weight);
        histMan.H1D_LepBinned("mass4jets")->Fill(tPlusXCandidates.mass4jets(), weight);
        histMan.H1D_LepBinned("HT_ltop")->Fill(tPlusXCandidates.fullHT(), weight);
        histMan.H1D_LepBinned("ST_ltop")->Fill(tPlusXCandidates.ST(), weight);
        // histMan.H1D_LepBinned("HT3jet_ltop")->Fill(tPlusXCandidates.jetsHTNoTop(), weight);
        // histMan.H1D_LepBinned("HT4jet_ltop")->Fill(tPlusXCandidates.jetsHTNoTop(4), weight);
        histMan.H1D_LepBinned("MET_ltop")->Fill(tPlusXCandidates.MET()->et(), weight);
        histMan.H1D_LepBinned("METpt_ltop")->Fill(tPlusXCandidates.MET()->pt(), weight);
        */
        // histMan.H1D("lone_electron_et")->Fill(tPlusXCandidates.getElectronFromWDecay()->et(), weight);
        // histMan.H1D_LepBinned("lone_neutrino_pz")->Fill(tPlusXCandidates.getNeutrinoFromWDecay()->pz(), weight);
        // histMan.H1D_LepBinned("pt_loneLepTop")->Fill(tPlusXCandidates.getLeptonicTop()->pt(), weight);

        // histMan.H1D_LepBinned("pt_lepTop")->Fill(twoTopsNonRes.getLeptonicTop()->pt(), weight);
        // histMan.H1D_LepBinned("pt_hadTop")->Fill(twoTopsNonRes.getHadronicTop()->pt(), weight);
        // histMan.H1D_LepBinned("angleHadTopD")->Fill(twoTopsNonRes.getHadTopDAngle(), weight);
        // histMan.H1D_LepBinned("mttbar")->Fill(twoTopsNonRes.mttbar(), weight);
				double goodWpMass = goodTopBadTop.wpFrom1Top->mass();
				double badWpMass = goodTopBadTop.wpBadFrom1Top->mass();
				if (goodWpMass < 0 || goodWpMass > 3000)
					goodWpMass = 0.0;
				if (badWpMass < 0 || badWpMass > 3000)
					badWpMass = 0.0;
				histMan.H1D_LepBinned("wp1topmass")->Fill(goodWpMass, weight);
				histMan.H1D_LepBinned("wp1topmassbad")->Fill(badWpMass, weight);
				histMan.H2D_BJetBinned("good_vs_badwp1topmass")->Fill(goodWpMass,
					badWpMass, weight);
        histMan.H1D_LepBinned("wpmass")->Fill(twoTopsNonRes.getWprime()->mass(), weight);
        histMan.H1D_LepBinned("wpnegmass")->Fill(twoTopsNonRes.getNegWprime()->mass(), weight);
        histMan.H1D_LepBinned("wpnegmassbad")->Fill(twoTopsNonRes.getNegWprimeBad()->mass(), weight);
        histMan.H1D_LepBinned("wpanglmass")->Fill(twoTopsNonRes.getBtBWprime()->mass(), weight);
        const ParticlePointer resonance = twoTopsNonRes.getResonance();
        histMan.H1D_BJetBinned("ttbar_pt")->Fill(resonance->pt(), weight);
				/*
				histMan.H2D_LepBinned("jet12vsjet13")->Fill(tPlusXCandidates.mass2jets(0, 1),
					tPlusXCandidates.mass2jets(0,2), weight);
				histMan.H2D_LepBinned("jet12vsjet23")->Fill(tPlusXCandidates.mass2jets(0, 1),
					tPlusXCandidates.mass2jets(1,2), weight);
				histMan.H2D_LepBinned("jet13vsjet23")->Fill(tPlusXCandidates.mass2jets(0, 2),
					tPlusXCandidates.mass2jets(1,2), weight);
				*/
				if (twoTopsNonRes.isRealData() == false) {
					histMan.H1D_LepBinned("tausInJets")->Fill(twoTopsNonRes.getTausInJets(),
						weight);
					histMan.H1D_LepBinned("numTaus")->Fill(twoTopsNonRes.numTaus(),
						weight);
					topTruth leptTop =
						twoTopsNonRes.getMCMatches(twoTopsNonRes.getLeptonicTop());
					topTruth hadTop =
						twoTopsNonRes.getMCMatches(twoTopsNonRes.getHadronicTop());
					histMan.H1D_BJetBinned("lepTopDeltaR")->Fill(leptTop.deltaR, weight);
					histMan.H1D_BJetBinned("hadTopDeltaR")->Fill(hadTop.deltaR, weight);
					histMan.H1D_BJetBinned("lepTopLep")->Fill(leptTop.leptonic, weight);
					histMan.H1D_BJetBinned("hadTopHad")->Fill(hadTop.hadronic, weight);
					histMan.H1D_BJetBinned("numHadTopMCMatches2top")->Fill(twoTopsNonRes.getNumMCMatchesHTop(), weight);
					histMan.H1D_BJetBinned("numHadTopCorrectID2top")->Fill(twoTopsNonRes.getNumCorrectIDHTop(), weight);
					histMan.H1D_LepBinned("numHadTopDauMCMatches2top")->Fill(twoTopsNonRes.goodHTopMatching(), weight);
					histMan.H1D_BJetBinned("numLepTopMCMatches2top")->Fill(twoTopsNonRes.getNumMCMatchesLTop(), weight);
					histMan.H1D_BJetBinned("numLepTopCorrectID2top")->Fill(twoTopsNonRes.getNumCorrectIDLTop(), weight);
					histMan.H1D_LepBinned("numLepTopDauMCMatches2top")->Fill(twoTopsNonRes.goodLTopMatching(), weight);
					/*
					
					int numStrict = goodTopBadTop.getNumMCMatchesHTop();
					int numLoose = goodTopBadTop.getNumCorrectIDHTop();
					histMan.H1D_BJetBinned("numHadTopMCMatches")->Fill(numStrict, weight);
					histMan.H1D_BJetBinned("numHadTopCorrectID")->Fill(numLoose, weight);
					histMan.H1D_LepBinned("numHadTopDauMCMatches")->Fill(goodTopBadTop.goodHTopMatching(), weight);
					numStrict = goodTopBadTop.getNumMCMatchesLTop();
					numLoose = goodTopBadTop.getNumCorrectIDLTop();
					histMan.H1D_BJetBinned("numLepTopMCMatches")->Fill(numStrict, weight);
					histMan.H1D_BJetBinned("numLepTopCorrectID")->Fill(numLoose, weight);		
					histMan.H1D_LepBinned("numLepTopDauMCMatches")->Fill(goodTopBadTop.goodLTopMatching(), weight);
					*/
				}
		}
    if (false && ttbarCandidate.passesFullTTbarEPlusJetSelection()) {
        histMan.H1D("numberOfBJets")->Fill(ttbarCandidate.GoodBJets().size(), weight);
        if(Event::usePFIsolation) {
            ttbarCandidate.reconstructTTbar(ttbarCandidate.GoodPFIsolatedElectrons().front());
            tPlusXCandidates.recoTprimeUsingChi2(ttbarCandidate.GoodPFIsolatedElectrons().front());
            tPlusXCandidates.recoBestSingleTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
            // loneTopsNoMassConstr.recoBestSingleTop(ttbarCandidate.GoodPFIsolatedElectrons().front());
        } else {
            ttbarCandidate.reconstructTTbar(ttbarCandidate.GoodIsolatedElectrons().front());
            tPlusXCandidates.recoTprimeUsingChi2(ttbarCandidate.GoodIsolatedElectrons().front());
            tPlusXCandidates.recoBestSingleTop(ttbarCandidate.GoodIsolatedElectrons().front());
            // loneTopsNoMassConstr.recoBestSingleTop(ttbarCandidate.GoodIsolatedElectrons().front());
			  }
        vector<TtbarHypothesisPointer> solutions = ttbarCandidate.Solutions();
        const ParticlePointer resonance = ttbarCandidate.getResonance();
        double mttbar = ttbarCandidate.mttbar();

        ParticlePointer leadingTop, nextToLeadingTop;

        if (ttbarCandidate.getHadronicTop()->pt() > ttbarCandidate.getLeptonicTop()->pt()) {
            leadingTop = ttbarCandidate.getHadronicTop();
            nextToLeadingTop = ttbarCandidate.getLeptonicTop();
        } else {
            leadingTop = ttbarCandidate.getLeptonicTop();
            nextToLeadingTop = ttbarCandidate.getHadronicTop();
        }
        double angleTops = leadingTop->angle(nextToLeadingTop);
        histMan.H1D_BJetBinned("angleTops")->Fill(angleTops, weight);
        histMan.H2D_BJetBinned("angleTops_vs_mttbar")->Fill(mttbar, angleTops, weight);
        histMan.H1D_BJetBinned("mLeptonicTop")->Fill(ttbarCandidate.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mHadronicTop")->Fill(ttbarCandidate.getHadronicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mAllTop")->Fill(ttbarCandidate.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mAllTop")->Fill(ttbarCandidate.getHadronicTop()->mass(), weight);
        histMan.H1D_BJetBinned("chiHadronicTop")->Fill(solutions.at(0)->hadronicChi2,
        	weight);
        histMan.H1D_BJetBinned("chiLeptonicTop")->Fill(solutions.at(0)->leptonicChi2,
        	weight);
        histMan.H1D_BJetBinned("chiGlobal")->Fill(solutions.at(0)->globalChi2,
        	weight);
        histMan.H1D_BJetBinned("chiTotal")->Fill(solutions.at(0)->totalChi2,
        	weight);
        histMan.H1D_BJetBinned("mLepTopLone")->Fill(tPlusXCandidates.getLeptonicTop()->mass(), weight);
        histMan.H1D_BJetBinned("mHadTopLone")->Fill(tPlusXCandidates.getHadronicTop()->mass(), weight);
        // histMan.H1D_BJetBinned("mLepTopLoneNM")->Fill(loneTopsNoMassConstr.getLeptonicTop()->mass(), weight);
        // histMan.H1D_BJetBinned("mHadTopLoneNM")->Fill(loneTopsNoMassConstr.getHadronicTop()->mass(), weight);

        histMan.H1D_BJetBinned("MET")->Fill(ttbarCandidate.MET()->et(), weight);
        histMan.H2D_BJetBinned("METvsMttbar")->Fill(mttbar, ttbarCandidate.MET()->et(), weight);
        histMan.H1D_BJetBinned("HT")->Fill(ttbarCandidate.fullHT(), weight);
        histMan.H2D_BJetBinned("HTvsMttbar")->Fill(mttbar, ttbarCandidate.fullHT(), weight);
				histMan.H1D_BJetBinned("leadingJetMass")->Fill(ttbarCandidate.GoodJets().front()->mass(), weight);
        if (ttbarCandidate.isRealData() == false) {
						histMan.H1D_BJetBinned("leadingGenJetMass")->Fill(ttbarCandidate.GenJets().front()->mass(), weight);
						topTruth leptTop = tPlusXCandidates.getMCMatches(tPlusXCandidates.getLeptonicTop());
						topTruth hadTop = tPlusXCandidates.getMCMatches(tPlusXCandidates.getHadronicTop());
						histMan.H1D_BJetBinned("lepTopDeltaR")->Fill(leptTop.deltaR, weight);
						histMan.H1D_BJetBinned("hadTopDeltaR")->Fill(hadTop.deltaR, weight);
						histMan.H1D_BJetBinned("lepTopLep")->Fill(leptTop.leptonic, weight);
						histMan.H1D_BJetBinned("hadTopHad")->Fill(hadTop.hadronic, weight);
						histMan.H1D_BJetBinned("numHadTopMCMatches")->Fill(tPlusXCandidates.getNumMCMatchesHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopMCMatches")->Fill(tPlusXCandidates.getNumMCMatchesLTop(), weight);
						histMan.H1D_BJetBinned("numHadTopCorrectID")->Fill(tPlusXCandidates.getNumCorrectIDHTop(), weight);
						histMan.H1D_BJetBinned("numLepTopCorrectID")->Fill(tPlusXCandidates.getNumCorrectIDLTop(), weight);
						// topTruth leptTopNM = loneTopsNoMassConstr.getMCMatches(loneTopsNoMassConstr.getLeptonicTop());
						// topTruth hadTopNM = loneTopsNoMassConstr.getMCMatches(loneTopsNoMassConstr.getHadronicTop());
						// histMan.H1D_BJetBinned("lepTopDeltaRNM")->Fill(leptTopNM.deltaR, weight);
						// histMan.H1D_BJetBinned("hadTopDeltaRNM")->Fill(hadTopNM.deltaR, weight);
						// histMan.H1D_BJetBinned("numHadTopMCMatchesNM")->Fill(loneTopsNoMassConstr.getNumMCMatchesHTop(), weight);
						// histMan.H1D_BJetBinned("numLepTopMCMatchesNM")->Fill(loneTopsNoMassConstr.getNumMCMatchesLTop(), weight);
						// histMan.H1D_BJetBinned("numHadTopCorrectIDNM")->Fill(loneTopsNoMassConstr.getNumCorrectIDHTop(), weight);
						// histMan.H1D_BJetBinned("numLepTopCorrectIDNM")->Fill(loneTopsNoMassConstr.getNumCorrectIDLTop(), weight);
				}

        if (Event::usePFIsolation)
            histMan.H1D_BJetBinned("mtW")->Fill(ttbarCandidate.transverseWmass(
                    ttbarCandidate.GoodPFIsolatedElectrons().front()), weight);
        else
            histMan.H1D_BJetBinned("mtW")->Fill(ttbarCandidate.transverseWmass(
                    ttbarCandidate.GoodIsolatedElectrons().front()), weight);
        histMan.H1D_BJetBinned("m3")->Fill(ttbarCandidate.M3(), weight);


//        histMan.H1D("mttbar")->Fill(mttbar, weight);
        histMan.H1D_BJetBinned("mttbar")->Fill(mttbar, weight);
				/*
        histMan.H1D_BJetBinned("tPrimeMass")->Fill(tPlusXCandidates.tpmass(), weight);
        histMan.H1D_BJetBinned("tPrimeHT")->Fill(tPlusXCandidates.TPrimeHTSystem(), weight);
        histMan.H1D_BJetBinned("tPrimepT")->Fill(tPlusXCandidates.PtTPrimeSystem(), weight);
        histMan.H1D_BJetBinned("tPrime_pt")->Fill((tPlusXCandidates.getResonance())->pt(), weight);
        histMan.H1D_BJetBinned("tPrime_px")->Fill((tPlusXCandidates.getResonance())->px(), weight);
        histMan.H1D_BJetBinned("tPrime_py")->Fill((tPlusXCandidates.getResonance())->py(), weight);
        histMan.H1D_BJetBinned("tPrime_pz")->Fill((tPlusXCandidates.getResonance())->pz(), weight);
        histMan.H2D_BJetBinned("tPrime_pt_vs_tPrimeM")->Fill(tPlusXCandidates.tpmass(), (tPlusXCandidates.getResonance())->pt(), weight);
				*/

        histMan.H1D_BJetBinned("ttbar_pt")->Fill(resonance->pt(), weight);

        histMan.H1D_BJetBinned("ttbar_px")->Fill(resonance->px(), weight);
        histMan.H1D_BJetBinned("ttbar_py")->Fill(resonance->py(), weight);
        histMan.H1D_BJetBinned("ttbar_pz")->Fill(resonance->pz(), weight);

        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar")->Fill(mttbar, resonance->pt(), weight);


        histMan.H1D("electron_et")->Fill(ttbarCandidate.getLeptonFromWDecay()->et(), weight);
        histMan.H1D_BJetBinned("neutrino_pz")->Fill(ttbarCandidate.getNeutrinoFromWDecay()->pz(), weight);

        histMan.H1D_BJetBinned("pt_leadingTop")->Fill(leadingTop->pt(), weight);
        histMan.H1D_BJetBinned("pt_NextToLeadingTop")->Fill(nextToLeadingTop->pt(), weight);
        histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar")->Fill(mttbar, leadingTop->pt(), weight);
        histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar")->Fill(mttbar, nextToLeadingTop->pt(), weight);

        if (ttbarCandidate.MET()->pt() > 20) {
            histMan.H1D_BJetBinned("angleTops_withMETCut")->Fill(angleTops, weight);
            histMan.H2D_BJetBinned("angleTops_vs_mttbar_withMETCut")->Fill(mttbar, angleTops, weight);
            histMan.H1D_BJetBinned("mttbar_withMETCut")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_withMETCut")->Fill(resonance->pt(), weight);
            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withMETCut")->Fill(mttbar, resonance->pt(), weight);

            histMan.H1D_BJetBinned("pt_leadingTop_withMETCut")->Fill(leadingTop->pt(), weight);
            histMan.H1D_BJetBinned("pt_NextToLeadingTop_withMETCut")->Fill(nextToLeadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETCut")->Fill(mttbar, leadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETCut")->Fill(mttbar, nextToLeadingTop->pt(), weight);

            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("angleTops_withMETAndAsymJets")->Fill(angleTops, weight);
                histMan.H2D_BJetBinned("angleTops_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, angleTops, weight);
                histMan.H1D_BJetBinned("mttbar_withMETAndAsymJets")->Fill(mttbar, weight);
                histMan.H1D_BJetBinned("ttbar_pt_withMETAndAsymJets")->Fill(resonance->pt(), weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, resonance->pt(), weight);

                histMan.H1D_BJetBinned("pt_leadingTop_withMETAndAsymJets")->Fill(leadingTop->pt(), weight);
                histMan.H1D_BJetBinned("pt_NextToLeadingTop_withMETAndAsymJets")->Fill(nextToLeadingTop->pt(), weight);
                histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, leadingTop->pt(), weight);
                histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets")->Fill(mttbar, nextToLeadingTop->pt(), weight);
            }
        }

        if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
            histMan.H1D_BJetBinned("angleTops_withAsymJetsCut")->Fill(angleTops, weight);
            histMan.H2D_BJetBinned("angleTops_vs_mttbar_withAsymJetsCut")->Fill(mttbar, angleTops, weight);
            histMan.H1D_BJetBinned("mttbar_withAsymJetsCut")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_withAsymJetsCut")->Fill(resonance->pt(), weight);

            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_withAsymJetsCut")->Fill(mttbar, resonance->pt(), weight);

            histMan.H1D_BJetBinned("pt_leadingTop_withAsymJetsCut")->Fill(leadingTop->pt(), weight);
            histMan.H1D_BJetBinned("pt_NextToLeadingTop_withAsymJetsCut")->Fill(nextToLeadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_leadingTop_vs_mttbar_withAsymJetsCut")->Fill(mttbar, leadingTop->pt(), weight);
            histMan.H2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut")->Fill(mttbar, nextToLeadingTop->pt(), weight);
        }

        for (unsigned int solutionIndex = 0; solutionIndex < solutions.size(); ++solutionIndex) {
            histMan.H1D_BJetBinned("mttbar_allSolutions")->Fill(solutions.at(solutionIndex)->resonance->mass(), weight);
            histMan.H1D_BJetBinned("ttbar_pt_allSolutions")->Fill(solutions.at(solutionIndex)->resonance->pt(), weight);
            histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions")->Fill(
                    solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                    weight);
            if (solutionIndex == 1) {
                histMan.H1D_BJetBinned("mttbar_2ndSolution")->Fill(solutions.at(solutionIndex)->resonance->mass(),
                        weight);

                histMan.H1D_BJetBinned("ttbar_pt_2ndSolution")->Fill(solutions.at(solutionIndex)->resonance->pt(),
                        weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);

            }
            if (solutionIndex == 2) {
                histMan.H1D_BJetBinned("mttbar_3rdSolution")->Fill(solutions.at(solutionIndex)->resonance->mass(),
                        weight);

                histMan.H1D_BJetBinned("ttbar_pt_3rdSolution")->Fill(solutions.at(solutionIndex)->resonance->pt(),
                        weight);

                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);
            }

            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), weight);
                histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->pt(), weight);
                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(),
                        weight);
                if (solutionIndex == 1) {
                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_2ndSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }
                if (solutionIndex == 2) {
                    histMan.H1D_BJetBinned("mttbar_3rdSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withMETCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }

                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_BJetBinned("mttbar_allSolutions_withMETAndAsymJets")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withMETAndAsymJets")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETAndAsymJets")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                }
            }

            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("mttbar_allSolutions_withAsymJetsCut")->Fill(
                        solutions.at(solutionIndex)->resonance->mass(), weight);
                histMan.H1D_BJetBinned("ttbar_pt_allSolutions_withAsymJetsCut")->Fill(
                        solutions.at(solutionIndex)->resonance->pt(), weight);
                histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withAsymJetsCut")->Fill(solutions.at(
                        solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);

                if (ttbarCandidate.MET()->et() > 20) {

                    if (solutionIndex == 1) {
                        histMan.H1D_BJetBinned("mttbar_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), weight);

                        histMan.H1D_BJetBinned("ttbar_pt_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->pt(), weight);
                        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    }

                    if (solutionIndex == 2) {
                        histMan.H1D_BJetBinned("mttbar_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), weight);

                        histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->pt(), weight);
                        histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETAndAsymJets")->Fill(solutions.at(
                                solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    }
                }

                if (solutionIndex == 1) {
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withAsymJetsCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);

                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                    histMan.H1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->mass(), weight);
                }

                if (solutionIndex == 2) {
                    histMan.H2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withAsymJetsCut")->Fill(solutions.at(
                            solutionIndex)->resonance->mass(), solutions.at(solutionIndex)->resonance->pt(), weight);
                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);

                    histMan.H1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut")->Fill(
                            solutions.at(solutionIndex)->resonance->pt(), weight);
                }

            }

            //            cout << "total Chi2 = " << solutions.at(solutionIndex)->totalChi2;
            //            cout << ", mass = " << solutions.at(solutionIndex)->resonance->mass() << endl;
        }

        if (ttbarCandidate.MET()->et() < 20) {
            histMan.H1D_BJetBinned("mttbar_QCDEnriched")->Fill(mttbar, weight);
            histMan.H1D_BJetBinned("ttbar_pt_QCDEnriched")->Fill(resonance->pt());
        }



//        cout << "Number of solutions: " << solutions.size() << " (" << ttbarCandidate.GoodJets().size() << " jets)"
//                << endl;
//        cout << "First solution, compare to chosen one" << endl;
//        cout << solutions.front()->totalChi2 << ", " << ttbarCandidate.getTotalChi2() << endl;

        if (mttbar != mttbar) {//isnan
            ttbarCandidate.inspectReconstructedEvent();
        }
        if (ttbarCandidate.isRealData()) {
            cout << "run " << ttbarCandidate.runnumber() << ", event " << ttbarCandidate.eventnumber() << ", lumi "
                    << ttbarCandidate.lumiblock();
            cout << ", top pair invariant mass = " << mttbar << " GeV" << endl;
            interestingEvents.push_back(InterestingEvent(ttbarCandidate, eventReader->getCurrentFile()));

            if (resonance->pt() > 100) {
                cout << "top pair pt = " << resonance->pt() << " GeV" << endl;
                ttbarCandidate.inspect();
                ttbarCandidate.inspectReconstructedEvent();
            }
        }

    }
}

void Analysis::doMicroNtuple() {
//        if( !goodTopBadTop.passesFullEPlusJetSelection() ) return;

        bool elecEvt = twoTopsNonRes.passesFullEPlusJetSelection();
        bool muEvt   = twoTopsNonRes.passesFullMuPlusJetSelection();

        if(!elecEvt && !muEvt) return;

        // make sure top reconstruction was ran 
        if( !goodTopBadTop.isReconstructed() ){
           cout<<"Reconstruction wasn't done: nothing to save into the micro ntuple"<<endl;
           return;
        }

        LeptonPointer lepton;
        if( twoTopsNonRes.GoodPFIsolatedElectrons().size() >= 1 )
            lepton = twoTopsNonRes.GoodPFIsolatedElectrons().front();
        if( twoTopsNonRes.GoodIsolatedMuons().size() >= 1 )
            lepton = twoTopsNonRes.GoodIsolatedMuons().front();

        // general information oabout this event:
        type          = currentEvent.getDataType();
        eventNumber   = currentEvent.eventnumber();
        runNumber     = currentEvent.runnumber();
        numberOfJets  = goodTopBadTop.GoodJets().size();
        numberOfBJets = goodTopBadTop.GoodBJets().size();
        leptoCharge   = (lepton->charge()<0?-1:1);
        ST            = goodTopBadTop.ST();
        HT            = goodTopBadTop.fullHT();
        MET           = goodTopBadTop.MET()->et();
        eSEL          = 0;
	nPileUpVtx    = currentEvent.numberOfGeneratedPileUpVertices();
        for(unsigned int cut=7; cut < toplikeElectronSelSize; ++cut) { // start with at least 1 b-tag
            if(goodTopBadTop.passesSelectionStepUpTo(cut, toplikeElectronSelection)) {
               eSEL = cut;
            } else break;
        }
        muSEL          = 0;
        for(unsigned int cut=7; cut < toplikeMuonSelSize; ++cut) { // start with at least 1 b-tag
            if(goodTopBadTop.passesSelectionStepUpTo(cut, toplikeMuonSelection)) {
               muSEL = cut;
            } else break;
        }

        // jets in the event:
        //  leading non-b-jet jet parameters and highest pT leftover jet (may not be the same jet):
        leadingJetPtRec = -1; leadingJetEtaRec = 0; leadingJetPhiRec = 0;
           freeJetPtRec = -1;    freeJetEtaRec = 0;    freeJetPhiRec = 0;
        int iDeanJet = -1;
        for( JetCollection::const_iterator jet  = goodTopBadTop.GoodJets().begin();
                                           jet != goodTopBadTop.GoodJets().end();
                                           jet ++ ){
           if( *jet != goodTopBadTop.getLeptonicBJet() ){
              if( leadingJetPtRec  < (*jet)->pt() ){
                  leadingJetPtRec  = (*jet)->pt();
                  leadingJetEtaRec = (*jet)->eta();
                  leadingJetPhiRec = (*jet)->phi();
              }
           }
           ++iDeanJet;
           JetPx[iDeanJet] = (*jet)->px();
           JetPy[iDeanJet] = (*jet)->py();
           JetPz[iDeanJet] = (*jet)->pz();
           JetE[iDeanJet] = (*jet)->energy();
           JetBTag[iDeanJet] = (*jet == goodTopBadTop.getLeptonicBJet() ? 1 : 0);
	   std::cout<<(*jet)->pt()<<endl;

        }
	std::cout<<currentEvent.GoodBJets().size()<<"  "<<goodTopBadTop.getLeptonicBJet()->pt()<<std::endl;
	for (unsigned int k=0; k < currentEvent.GoodBJets().size(); k++){
	  std::cout<<currentEvent.GoodBJets().at(k)->pt()<<std::endl;
	  BJetPx[k] = currentEvent.GoodBJets().at(k)->px();
	  BJetPy[k] = currentEvent.GoodBJets().at(k)->py();
	  BJetPz[k] = currentEvent.GoodBJets().at(k)->pz();
	  BJetE[k] = currentEvent.GoodBJets().at(k)->energy();

	}

	// std::cout << goodTopBadTop.fullHT() << std::endl;
	std::cout<<"----------------------------"<<endl;
        freeJetPtRec  = goodTopBadTop.dJetFromWp->pt();
        freeJetEtaRec = goodTopBadTop.dJetFromWp->eta();
        freeJetPhiRec = goodTopBadTop.dJetFromWp->phi();

        leptonPtRec  = lepton->pt();
        leptonEtaRec = lepton->eta();
        leptonPhiRec = lepton->phi();

        bJetTlPtRec  = goodTopBadTop.getLeptonicBJet()->pt();
        bJetTlEtaRec = goodTopBadTop.getLeptonicBJet()->eta();
        bJetTlPhiRec = goodTopBadTop.getLeptonicBJet()->phi();

        bJetThPtRec  = goodTopBadTop.getHadronicBJet()->pt();
        bJetThEtaRec = goodTopBadTop.getHadronicBJet()->eta();
        bJetThPhiRec = goodTopBadTop.getHadronicBJet()->phi();

        jet1WhPtRec  = goodTopBadTop.getJet1FromHadronicW()->pt();
        jet1WhEtaRec = goodTopBadTop.getJet1FromHadronicW()->eta();
        jet1WhPhiRec = goodTopBadTop.getJet1FromHadronicW()->phi();
        jet2WhPtRec  = goodTopBadTop.getJet2FromHadronicW()->pt();
        jet2WhEtaRec = goodTopBadTop.getJet2FromHadronicW()->eta();
        jet2WhPhiRec = goodTopBadTop.getJet1FromHadronicW()->phi();

        metPtRec     = goodTopBadTop.MET()->pt();
        metPhiRec    = goodTopBadTop.MET()->phi();
        metEtaRec    = goodTopBadTop.getNeutrinoFromWDecay()->eta();

        whMrec   = goodTopBadTop.getHadronicW()->mass();
        whPtRec  = goodTopBadTop.getHadronicW()->pt();
        whEtaRec = goodTopBadTop.getHadronicW()->eta();
        whPhiRec = goodTopBadTop.getHadronicW()->phi();

        thMrec   = goodTopBadTop.getHadronicTop()->mass();
        thPtRec  = goodTopBadTop.getHadronicTop()->pt();
        thEtaRec = goodTopBadTop.getHadronicTop()->eta();
        thPhiRec = goodTopBadTop.getHadronicTop()->phi();

        wlMrec   = goodTopBadTop.getLeptonicW()->mass();
        wlPtRec  = goodTopBadTop.getLeptonicW()->pt();
        wlEtaRec = goodTopBadTop.getLeptonicW()->eta();
        wlPhiRec = goodTopBadTop.getLeptonicW()->phi();

        tlMrec   = goodTopBadTop.getLeptonicTop()->mass();
        tlPtRec  = goodTopBadTop.getLeptonicTop()->pt();
        tlEtaRec = goodTopBadTop.getLeptonicTop()->eta();
        tlPhiRec = goodTopBadTop.getLeptonicTop()->phi();

        tPosJetMassRec = goodTopBadTop.wpBadFrom1Top->mass();
        tPosJetPtRec   = goodTopBadTop.wpBadFrom1Top->pt();
        tPosJetEtaRec  = goodTopBadTop.wpBadFrom1Top->eta();
        tPosJetPhiRec  = goodTopBadTop.wpBadFrom1Top->phi();
        
        tNegJetMassRec = goodTopBadTop.wpFrom1Top->mass();
        tNegJetPtRec   = goodTopBadTop.wpFrom1Top->pt();
        tNegJetEtaRec  = goodTopBadTop.wpFrom1Top->eta();
        tNegJetPhiRec  = goodTopBadTop.wpFrom1Top->phi();

        microTuple->Fill();
}


void Analysis::doNotePlots() {
    if (ttbarCandidate.GoodElectrons().size() >= 1 && ttbarCandidate.Jets().size() >= 2) {
        const ElectronCollection electrons = ttbarCandidate.GoodElectrons();
        ElectronCollection nonConversionElectrons;
        for (unsigned int index = 0; index < electrons.size(); ++index) {
            const ElectronPointer electron = electrons.at(index);
            if (electron->isFromConversion() == false && electron->isTaggedAsConversion(0.02,0.02) == false) {
//                ConversionTagger tagger = ConversionTagger();
//                tagger.calculateConversionVariables(electron, ttbarCandidate.Tracks(), 3.8, 0.45);
//                if (tagger.isFromConversion(0.02, 0.02) == false)
                    nonConversionElectrons.push_back(electron);
            }
        }
        if (nonConversionElectrons.size() == 1) {
            const ElectronPointer electron = nonConversionElectrons.front();
            JetCollection goodjets;
            for (unsigned index = 0; index < ttbarCandidate.Jets().size(); ++index) {
                if (ttbarCandidate.Jets().at(index)->isGood())
                    goodjets.push_back(ttbarCandidate.Jets().at(index));
            }
            if (goodjets.size() >= 2) {
                unsigned int closestID = electron->getClosestJetIndex(goodjets);
                float minDR = electron->deltaR(goodjets.at(closestID));
                float ptRel = electron->relativePtTo(goodjets.at(closestID));
                histMan.H2D("ptRel_vs_DRmin")->Fill(minDR, ptRel, weight);
                if (ttbarCandidate.MET()->et() < 20 && ttbarCandidate.transverseWmass(electron) < 35) {
                    histMan.H1D("DRmin_QCDenriched")->Fill(minDR, weight);
                    histMan.H1D("ptRel_QCDenriched")->Fill(ptRel, weight);
                } else if (ttbarCandidate.MET()->et() > 30 && ttbarCandidate.transverseWmass(electron) > 50) {
                    histMan.H1D("DRmin_WZenriched")->Fill(minDR, weight);
                    histMan.H1D("ptRel_WZenriched")->Fill(ptRel, weight);
                }
            }
        }

    }
}

void Analysis::doQCDStudy() {
    if (ttbarCandidate.passesRelIsoSelection()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_CombRelIso")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_CombRelIso_1btag")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_CombRelIso_2btag")->Fill(electron->relativeIsolation(), weight);
    }

    if (ttbarCandidate.passesPFIsoSelection() && NTupleEventReader::electronAlgorithm
            == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_PFIsolation")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_PFIsolation_1btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_PFIsolation_2btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.MET()->pt() > 20) {
            histMan.H1D_JetBinned("QCDest_PFIsolation_WithMETCut")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodJets().size() >= 2) {
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_JetBinned("QCDest_PFIsolation_WithMETCutAndAsymJetCuts")->Fill(electron->pfIsolation(),
                            weight);
                }
            }
        }

        if (ttbarCandidate.GoodJets().size() >= 2) {
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_JetBinned("QCDest_PFIsolation_WithAsymJetCuts")->Fill(electron->pfIsolation(), weight);

            }
        }
    }


    if (ttbarCandidate.passesRelIsoControlSelection()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion_1btag")->Fill(electron->relativeIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_CombRelIso_controlRegion_2btag")->Fill(electron->relativeIsolation(), weight);

        if (NTupleEventReader::electronAlgorithm == ElectronAlgorithm::ParticleFlow) {
            const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodBJets().size() >= 1)
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_1btag")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.GoodBJets().size() >= 2)
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_2btag")->Fill(electron->pfIsolation(), weight);

            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCut")->Fill(electron->pfIsolation(), weight);

                if (ttbarCandidate.GoodJets().size() >= 2) {
                    if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                        histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCutAndAsymJetCuts")->Fill(
                                electron->pfIsolation(), weight);
                    }
                }
            }

            if (ttbarCandidate.GoodJets().size() >= 2) {
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithAsymJetCuts")->Fill(electron->pfIsolation(), weight);

                }
            }
        }

    }

    if (ttbarCandidate.passesPFIsoControlSelection() && NTupleEventReader::electronAlgorithm
            == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostPFIsolatedElectron();
        histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 1)
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion_1btag")->Fill(electron->pfIsolation(), weight);

        if (ttbarCandidate.GoodBJets().size() >= 2)
            histMan.H1D_JetBinned("QCDest_PFIsolation_controlRegion_2btag")->Fill(electron->pfIsolation(), weight);
    }

    if (ttbarCandidate.passesRelIsoSelection() && ttbarCandidate.hasAtLeastFourGoodJets()) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron(false);
        if (electron->isIsolated() == false && !isnan(electron->relativeIsolation()) && !isinf(
                electron->relativeIsolation())) {
            try {
                ttbarCandidate.reconstructTTbar(electron);
                const ParticlePointer resonance = ttbarCandidate.getResonance();
                histMan.H1D_BJetBinned("mttbar_controlRegion")->Fill(resonance->mass(), weight);
                if (ttbarCandidate.MET()->pt() > 20) {
                    histMan.H1D_BJetBinned("mttbar_controlRegion_withMETCut")->Fill(resonance->mass(), weight);
                    if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                        histMan.H1D_BJetBinned("mttbar_controlRegion_withMETAndAsymJets")->Fill(resonance->mass(),
                                weight);
                }

                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                    histMan.H1D_BJetBinned("mttbar_controlRegion_withAsymJetsCut")->Fill(resonance->mass(), weight);

            } catch (ReconstructionException& e) {
                cout << "Could not reconstruct event: " << e.what() << endl;
            }
        }
    }

    if (ttbarCandidate.passesConversionSelection()) {
        ElectronPointer electron;
        if (Event::usePFIsolation)
            electron = currentEvent.GoodPFIsolatedElectrons().front();
        else
            electron = currentEvent.GoodIsolatedElectrons().front();
        try {
            ttbarCandidate.reconstructTTbar(electron);
            const ParticlePointer resonance = ttbarCandidate.getResonance();
            histMan.H1D_BJetBinned("mttbar_conversions")->Fill(resonance->mass(), weight);
            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_conversions_withMETCut")->Fill(resonance->mass(), weight);
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                    histMan.H1D_BJetBinned("mttbar_conversions_withMETAndAsymJets")->Fill(resonance->mass(), weight);
            }
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50)
                histMan.H1D_BJetBinned("mttbar_conversions_withAsymJetsCut")->Fill(resonance->mass(), weight);
        } catch (ReconstructionException& e) {
            cout << "Could not reconstruct event: " << e.what() << endl;
        }

    }

    if (ttbarCandidate.Electrons().size() > 0 && ttbarCandidate.GoodPFIsolatedElectrons().size() < 2
            && NTupleEventReader::electronAlgorithm == ElectronAlgorithm::ParticleFlow) {
        const ElectronPointer electron = ttbarCandidate.MostIsolatedElectron(true);
        histMan.H1D_JetBinned("MostPFIsolatedElectron_dPhiIn")->Fill(electron->dPhiIn(), weight);
        histMan.H1D_JetBinned("MostPFIsolatedElectron_dEtaIn")->Fill(electron->dEtaIn(), weight);
    }

    if (ttbarCandidate.passesAntiIsolationSelection()) {
        ElectronPointer electron = ttbarCandidate.GoodElectrons().front();
        try {
            ttbarCandidate.reconstructTTbar(electron);
            float mttbar = ttbarCandidate.mttbar();
            histMan.H1D_BJetBinned("mttbar_antiIsolated")->Fill(mttbar, weight);
            if (ttbarCandidate.MET()->pt() > 20) {
                histMan.H1D_BJetBinned("mttbar_antiIsolated_withMETCut")->Fill(mttbar, weight);
                if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                    histMan.H1D_BJetBinned("mttbar_antiIsolated_withMETAndAsymJets")->Fill(mttbar, weight);
                }
            }
            if (ttbarCandidate.GoodJets().front()->pt() > 70 && ttbarCandidate.GoodJets().at(1)->pt() > 50) {
                histMan.H1D_BJetBinned("mttbar_antiIsolated_withAsymJetsCut")->Fill(mttbar, weight);
            }

        } catch (ReconstructionException& e) {
            cout << "Could not reconstruct event: " << e.what() << endl;
        }
    }

}

void Analysis::doMCttbarReconstruction() {
	MCParticlePointer top, antitop, b_from_top, b_from_antitop, W_plus, W_minus, electron, neutrino, quark_from_W, antiquark_from_W;
	JetCollection genJets = currentEvent.GenJets();
	JetPointer topBjet, antitopBjet, jet1fromW, jet2fromW;
	TtbarHypothesis MCttbarEvent;
	bool ejets_event = false;
	bool leptonic_Wplus_found = false, leptonic_Wminus_found = false;
	bool hadronic_Wplus_found = false, hadronic_Wminus_found = false;
	bool fully_hadronic_event = false, fully_leptonic_event = false;
	bool non_electron_leptonic_channel = false;
	int index = 0;
	int top_index = -100, antitop_index = -100, W_plus_index = -100, W_minus_index = -100, electron_index = -100, neutrino_index = -100,
			b_from_top_index = -100, b_from_antitop_index = -100, quark_from_W_index = -100, antiquark_from_W_index = -100;

	// MC ttbar reconstruction
	for (MCParticleCollection::const_iterator mc_particle = currentEvent.GenParticles().begin(); mc_particle != currentEvent.GenParticles().end(); ++mc_particle, ++index) {

		if ((*mc_particle)->status() != 3) continue;
		//top quark
		if ((*mc_particle)->pdgId() == 6) {
			top = *mc_particle;
			top_index = index;
			continue;
		}

		//anti-top quark
		if ((*mc_particle)->pdgId() == -6) {
			antitop = *mc_particle;
			antitop_index = index;
			continue;
		}

		//W� bosons
		if (((*mc_particle)->pdgId() == 24) && ((*mc_particle)->motherIndex() == top_index)) {
			W_plus = *mc_particle;
			W_plus_index = index;
			continue;
		}

		if (((*mc_particle)->pdgId() == -24) && ((*mc_particle)->motherIndex() == antitop_index)) {
			W_minus = *mc_particle;
			W_minus_index = index;
			continue;
		}

		//b-quarks
		if (((*mc_particle)->pdgId() == 5) && ((*mc_particle)->motherIndex() == top_index)) {
			b_from_top = *mc_particle;
			b_from_top_index = index;
			continue;
		}
		if (((*mc_particle)->pdgId() == -5) && ((*mc_particle)->motherIndex() == antitop_index)) {
			b_from_antitop = *mc_particle;
			b_from_antitop_index = index;
			continue;
		}

		//W+ decay products
		if ((*mc_particle)->motherIndex()==W_plus_index) {
			if ((*mc_particle)->pdgId() == -11) {
				electron = *mc_particle;
				electron_index = index;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->pdgId() == 12) {
				neutrino = *mc_particle;
				neutrino_index = index;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->isLepton()) {
				non_electron_leptonic_channel = true;
				leptonic_Wplus_found = true;
			}

			else if ((*mc_particle)->isQuark()  && ((*mc_particle)->pdgId()>0)) {
				quark_from_W = *mc_particle;
				quark_from_W_index = index;
				hadronic_Wplus_found = true;
			}

			else if ((*mc_particle)->isQuark() && ((*mc_particle)->pdgId()<0)) {
				antiquark_from_W = *mc_particle;
				antiquark_from_W_index = index;
				hadronic_Wplus_found = true;
			}

			else {
				cout << "Something went wrong: W+ has unusual decay products." << endl;
			}

		}

		//W- decay products
		if ((*mc_particle)->motherIndex()==W_minus_index) {
			if ((*mc_particle)->pdgId() == 11) {
				electron = *mc_particle;
				electron_index = index;
				leptonic_Wminus_found = true;
			}

			else if ((*mc_particle)->pdgId() == -12) {
				neutrino = *mc_particle;
				neutrino_index = index;
				leptonic_Wminus_found = true;
			}

			else if ((*mc_particle)->isLepton()) {
				leptonic_Wminus_found = true;
				non_electron_leptonic_channel = true;
			}

			else if ((*mc_particle)->isQuark()  && ((*mc_particle)->pdgId()>0)) {
				quark_from_W = *mc_particle;
				quark_from_W_index = index;
				hadronic_Wminus_found = true;
			}

			else if ((*mc_particle)->isQuark() && ((*mc_particle)->pdgId()<0)) {
				antiquark_from_W = *mc_particle;
				antiquark_from_W_index = index;
				hadronic_Wminus_found = true;
			}

			else {
				cout << "Something went wrong: W- has unusual decay products." << endl;
			}
		}
	}

	//classify the event
	if (((leptonic_Wplus_found) || (leptonic_Wminus_found)) && ((hadronic_Wplus_found) || (hadronic_Wminus_found))
			&& (!non_electron_leptonic_channel)) { ejets_event = true; }
	if (((leptonic_Wplus_found) || (leptonic_Wminus_found)) && (!hadronic_Wplus_found)
			&& (!hadronic_Wminus_found)) { fully_leptonic_event = true; }
	if (((hadronic_Wplus_found) || (hadronic_Wminus_found)) && (!leptonic_Wplus_found)
			&& (!leptonic_Wminus_found)) { fully_hadronic_event = true; }

	if (ejets_event) {
		//matching genJets and partons

		if (genJets.size()>0) {
			int closestJetQuarkFromWIndex = quark_from_W->getClosestJetIndex(genJets);
			float minDR_quarkW = quark_from_W->deltaR(genJets.at(closestJetQuarkFromWIndex));
			jet1fromW = genJets.at(closestJetQuarkFromWIndex);

			int closestJetAntiQuarkFromWIndex = antiquark_from_W->getClosestJetIndex(genJets);
			float minDR_antiquarkW = antiquark_from_W->deltaR(genJets.at(closestJetAntiQuarkFromWIndex));
			jet2fromW = genJets.at(closestJetAntiQuarkFromWIndex);

			int closestJetBfromTopIndex = b_from_top->getClosestJetIndex(genJets);
			float minDR_BfromTop = b_from_top->deltaR(genJets.at(closestJetBfromTopIndex));
			topBjet = genJets.at(closestJetBfromTopIndex);

			int closestJetBfromAntiTopIndex = b_from_antitop->getClosestJetIndex(genJets);
			float minDR_BfromAntiTop = b_from_antitop->deltaR(genJets.at(closestJetBfromAntiTopIndex));
			antitopBjet = genJets.at(closestJetBfromAntiTopIndex);

			//delta R between genJets and partons histograms
			histMan.H1D("deltaRjet1")->Fill(minDR_quarkW);
			histMan.H1D("deltaRjet2")->Fill(minDR_antiquarkW);
			histMan.H1D("deltaRjet3")->Fill(minDR_BfromTop);
			histMan.H1D("deltaRjet4")->Fill(minDR_BfromAntiTop);

			histMan.H1D("deltaRjet_sum")->Fill(minDR_quarkW);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_antiquarkW);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_BfromTop);
			histMan.H1D("deltaRjet_sum")->Fill(minDR_BfromAntiTop);
		}

		if (leptonic_Wplus_found) {
			MCttbarEvent.leptonicTop = (ParticlePointer) top;
			MCttbarEvent.hadronicTop = (ParticlePointer) antitop;
			MCttbarEvent.leptonicW = (ParticlePointer) W_plus;
			MCttbarEvent.hadronicW = (ParticlePointer) W_minus;
			MCttbarEvent.leptonicBjet = topBjet;
			MCttbarEvent.hadronicBJet = antitopBjet;
			MCttbarEvent.jet1FromW = jet1fromW;
			MCttbarEvent.jet2FromW = jet2fromW;
			MCttbarEvent.neutrinoFromW = (ParticlePointer) neutrino;
			ElectronPointer e(new Electron(electron->energy(), electron->px(), electron->py(), electron->pz()));
			MCttbarEvent.leptonFromW = e;
		}
		else if (hadronic_Wplus_found) {
			MCttbarEvent.leptonicTop = (ParticlePointer) antitop;
			MCttbarEvent.hadronicTop = (ParticlePointer) top;
			MCttbarEvent.leptonicW = (ParticlePointer) W_minus;
			MCttbarEvent.hadronicW = (ParticlePointer) W_plus;
			MCttbarEvent.leptonicBjet = antitopBjet;
			MCttbarEvent.hadronicBJet = topBjet;
			MCttbarEvent.jet1FromW = jet1fromW;
			MCttbarEvent.jet2FromW = jet2fromW;
			MCttbarEvent.neutrinoFromW = (ParticlePointer) neutrino;
			ElectronPointer e(new Electron(electron->energy(), electron->px(), electron->py(), electron->pz()));
			MCttbarEvent.leptonFromW = e;
		}
		else cout << "ERROR: no hadronic or leptonic W's in semileptonic event (nonsense).\n";

		//comparing deltaR between genJets from W and closest partons
		histMan.H2D("deltaR_genJets_partons")->Fill(MCttbarEvent.jet1FromW->deltaR(MCttbarEvent.jet2FromW),quark_from_W->deltaR(antiquark_from_W));

		//invariant mass histograms
		histMan.H1D("W_inv_mass_from_truth_partons")->Fill(quark_from_W->invariantMass(antiquark_from_W));
		histMan.H1D("W_inv_mass_from_genJets")->Fill(MCttbarEvent.jet1FromW->invariantMass(MCttbarEvent.jet2FromW));
		histMan.H1D("top_leptonic_inv_mass_from_truth")->Fill(MCttbarEvent.leptonicW->invariantMass(MCttbarEvent.leptonicBjet));
		histMan.H1D("top_hadronic_inv_mass_from_truth")->Fill(MCttbarEvent.hadronicW->invariantMass(MCttbarEvent.hadronicBJet));

		histMan.H1D("m3_mc")->Fill(MCttbarEvent.M3());

		// comparing truth and reco objects
		if (ttbarCandidate.passesFullTTbarEPlusJetSelection()) {
			histMan.H1D("m3_diff")->Fill(fabs(MCttbarEvent.M3()-ttbarCandidate.M3()));

			histMan.H1D("deltaRElectron")->Fill(MCttbarEvent.leptonFromW->deltaR(ttbarCandidate.getLeptonFromWDecay()));
			histMan.H1D("deltaRLeptonicBjet")->Fill(MCttbarEvent.leptonicBjet->deltaR(ttbarCandidate.getLeptonicBJet()));
			histMan.H1D("deltaRHadronicBjet")->Fill(MCttbarEvent.hadronicBJet->deltaR(ttbarCandidate.getHadronicBJet()));
			histMan.H1D("deltaRjet1fromW")->Fill(MCttbarEvent.jet1FromW->deltaR(ttbarCandidate.getJet1FromHadronicW()));
			histMan.H1D("deltaRjet2fromW")->Fill(MCttbarEvent.jet2FromW->deltaR(ttbarCandidate.getJet2FromHadronicW()));
		}
	}

}

void Analysis::printInterestingEvents() {
    cout << "Interesting events:" << endl;
    for (unsigned int index = 0; index < interestingEvents.size(); ++index) {
        interestingEvents.at(index).print();
    }
}

void Analysis::printSummary() {
    EventTablePrinter::printWprimeCutFlow(cutflowPerSample, weights,
			toplikeElectronSelection, toplikeElectronSelSize, "Electron", luminosity);
    EventTablePrinter::printWprimeCutFlowUnwt(cutflowPerSample, weights,
			toplikeElectronSelection, toplikeElectronSelSize, "Electron");
    EventTablePrinter::printWprimeCutFlow(mucutflowPerSample, weights,
			toplikeMuonSelection, toplikeMuonSelSize, "Muon", luminosity);
    EventTablePrinter::printWprimeCutFlowUnwt(mucutflowPerSample, weights,
			toplikeMuonSelection, toplikeMuonSelSize, "Muon");
    // EventTablePrinter::printCutFlowLatexTable(cutflowPerSample);
    // EventTablePrinter::printUnweightedCutFlowLatexTable(cutflowPerSample);

    cout << "total number of processed events: " << eventReader->getNumberOfProccessedEvents() << endl;
    cout << endl;
    for (unsigned int cut = 0; cut < toplikeElectronSelSize; ++cut) {
        cout << "Selection step '" << ToplikeSelectionSteps::StringSteps[cut] << "'" << endl;
        cout << "passed events (single cut): " << singleCuts.at(cut) << endl;
        if (cut < toplikeElectronSelSize - 2)
            cout << "passed events (up to this cut):" << cutflow.at(cut) << endl;
        else
            cout << "passed events (full selection):" << cutflow.at(cut) << endl;
        cout << endl;
    }

    cout << "number of events without electrons: " << brokenEvents.size() << endl;
    cout << "\n\n Number of positive, negative W' tops: " << wpTopPos << ", " << wpTopNeg << endl;
    cout << "\n\n Number of positive, negative spectator tops: " << specTopPos << ", " << specTopNeg << endl;
}

void Analysis::createHistograms() {
    histMan.setCurrentLumi(Analysis::luminosity);
    histMan.prepareForSeenDataTypes(eventReader->getSeenDatatypes());
/*
    //histograms for Jet study
    histMan.addH1D_BJetBinned("AllJetMass", "AllJetMass", 500, 0, 500);
    histMan.addH1D_BJetBinned("AllGoodJetMass", "AllGoodJetMass", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastOneJets", "GoodJetMass_atLeastOneJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastTwoJets", "GoodJetMass_atLeastTwoJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastThreeJets", "GoodJetMass_atLeastThreeJets", 500, 0, 500);
    histMan.addH1D_BJetBinned("GoodJetMass_atLeastFourJets", "GoodJetMass_atLeastFourJets", 500, 0, 500);
    //MC histograms
    histMan.addH1D("deltaRElectron", "delta R between truth and reco electron", 100, 0, 0.2);
    histMan.addH1D("deltaRLeptonicBjet", "delta R between truth and reco b-jet on leptonic side", 100, 0, 0.5);
    histMan.addH1D("deltaRHadronicBjet", "delta R between truth and reco b-jet on hadronic side", 100, 0, 0.5);
    histMan.addH1D("deltaRjet1fromW", "delta R between truth and reco jet1 from W decay", 100, 0, 0.5);
    histMan.addH1D("deltaRjet2fromW", "delta R between truth and reco jet2 from W decay", 100, 0, 0.5);

    histMan.addH1D("deltaRjet1", "delta R between quark from W and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet2", "delta R between antiquark from W and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet3", "delta R between b quark from top and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet4", "delta R between b quark from antitop and closest genJet", 100, 0, 0.5);
    histMan.addH1D("deltaRjet_sum", "summarized delta R between partons and genJets", 100, 0, 0.5);

    histMan.addH2D("deltaR_genJets_partons", "delta R between genJets from W as opposed to partons", 100, 0, 5, 100, 0, 5);

    histMan.addH1D("W_inv_mass_from_truth_partons", "W inv. mass from truth partons", 100, 0, 120);
    histMan.addH1D("W_inv_mass_from_genJets", "W inv. mass from genJets", 100, 0, 120);
    histMan.addH1D("top_leptonic_inv_mass_from_truth", "Leptonic top inv. mass from truth partons", 100, 100, 220);
    histMan.addH1D("top_hadronic_inv_mass_from_truth", "Haronic top inv. mass from truth partons", 100, 100, 220);
    histMan.addH1D("m3_mc", "M3 for truth event", 500, 0, 500);
    histMan.addH1D("m3_diff", "M3 difference between truth and reco", 500, 0, 500);

    //end of MC histograms

		*/
    histMan.addH1D("electron_et", "electron_et", 600, 0, 600);
    histMan.addH1D("electron_pt", "electron_pt", 600, 0, 600);
    histMan.addH1D("electminus_pt", "electminus_pt", 600, 0, 600);
    histMan.addH1D("positron_pt", "positron_pt", 600, 0, 600);
    histMan.addH1D("muon_et", "muon_et", 600, 0, 600);
    histMan.addH1D("muon_pt", "muon_pt", 600, 0, 600);
    histMan.addH1D("lone_electron_et", "lone_electron_et", 600, 0, 600);
/*    histMan.addH1D_JetBinned("MostPFIsolatedElectron_dPhiIn", "MostPFIsolatedElectron_dPhiIn", 50, 0, 0.1);
    histMan.addH1D_JetBinned("MostPFIsolatedElectron_dEtaIn", "MostPFIsolatedElectron_dEtaIn", 50, 0, 0.02 );
    histMan.addH1D_JetBinned("diElectronMass", "diElectronMass", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("tPrimeMass", "tPrimeMass", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrimepT", "tPrimepT", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrimeHT", "tPrimeHT", 600, 0, 2.0);

    histMan.addH1D_BJetBinned("mttbar_conversions", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_QCDEnriched", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withMETCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withMETCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withMETCut", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withMETAndAsymJets", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withMETAndAsymJets", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withMETAndAsymJets", "mttbar", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_conversions_withAsymJetsCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_antiIsolated_withAsymJetsCut", "mttbar", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_controlRegion_withAsymJetsCut", "mttbar", 5000, 0, 5000);
*/
    // histMan.addH1D_BJetBinned("mttbar", "mttbar", 6000, 0, 6000);
    histMan.addH1D("eleciso", "eleciso", 600, 0, 3);
    histMan.addH1D("looseeleciso", "looseeleciso", 600, 0, 3);
    histMan.addH1D("muiso", "muiso", 600, 0, 3);
    histMan.addH1D("elvertdist", "elvertdist", 3, 0, 3);
    histMan.addH1D("muvertdist", "muvertdist", 3, 0, 3);
    histMan.addH1D_LepBinned("wpmass", "wpmass", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpnegmass", "wpnegmass", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpnegmassbad", "wpnegmassbad", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpnegmassneg", "wpnegmassneg", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpnegmasspos", "wpnegmasspos", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmass", "wp1topmass", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmassbad", "wp1topmassbad", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmasspos", "wp1topmasspos", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmassneg", "wp1topmassneg", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmassbadpos", "wp1topmassbadpos", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wp1topmassbadneg", "wp1topmassbadneg", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpanglmass", "wpanglmass", 6000, 0, 6000);
    histMan.addH1D_LepBinned("transMasslm1b", "transMasslm1b", 6000, 0, 6000);
    histMan.addH1D("jetmultip", "jetmultip", 10, 0, 10);
    histMan.addH1D_LepBinned("tausInJets", "number of taus in jets", 10, 0, 10);
    histMan.addH1D_LepBinned("numTaus", "number of taus in event", 20, 0, 20);
    histMan.addH1D_LepBinned("deltaPhiTops", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_LepBinned("deltaphilepjetpos", "angle between lepton & jet", 400, 0, 4);
    histMan.addH1D_LepBinned("deltaphilepdjetpos", "angle between lepton & jet", 400, 0, 4);
    histMan.addH1D_LepBinned("deltaphilepjetneg", "angle between lepton & jet", 400, 0, 4);
    histMan.addH1D_LepBinned("deltaphilepdjetneg", "angle between lepton & jet", 400, 0, 4);
    histMan.addH1D_LepBinned("wpchgasymm", "wpchgasymm", 7, 0, 6);
    histMan.addH1D_LepBinned("angleWpTopDTru", "angle between top and d", 400, 0, 4);
    histMan.addH1D_LepBinned("angleSpecTopDTru", "angle between top and d", 400, 0, 4);
    histMan.addH1D_LepBinned("wpTruPt", "wpTruPt", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpTopTruPt", "wpTopTruPt", 6000, 0, 6000);
    histMan.addH1D_LepBinned("specTopTruPt", "specTopTruPt", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpTopBTruPt", "wpTopTruPt", 6000, 0, 6000);
    histMan.addH1D_LepBinned("specTopBTruPt", "specTopTruPt", 6000, 0, 6000);
    histMan.addH1D_LepBinned("wpDDauTruPt", "wpDDauTruPt", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("wpmassplus", "wpmassplus", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("wpmassminus", "wpmassminus", 6000, 0, 6000);
/*    histMan.addH1D_BJetBinned("mttbar_withMETCut", "mttbar_withMETCut", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_withMETAndAsymJets", "mttbar_withMETAndAsymJets", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_withAsymJetsCut", "mttbar_withAsymJetsCut", 5000, 0, 5000);
//
    histMan.addH1D_BJetBinned("mttbar_2ndSolution", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions", "mttbar_allSolutions", 5000, 0, 5000);
//
    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withMETCut", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withMETCut", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withMETCut", "mttbar_allSolutions", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withMETAndAsymJets", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withMETAndAsymJets", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withMETAndAsymJets", "mttbar_allSolutions", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mttbar_2ndSolution_withAsymJetsCut", "mttbar_2ndSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_3rdSolution_withAsymJetsCut", "mttbar_3rdSolution", 5000, 0, 5000);
    histMan.addH1D_BJetBinned("mttbar_allSolutions_withAsymJetsCut", "mttbar_allSolutions", 5000, 0, 5000);
*/
    histMan.addH1D_LepBinned("mLeptonicTop", "mLeptonicTop", 1200, 0, 1200);
    histMan.addH1D_LepBinned("mHadronicTop", "mHadronicTop", 1200, 0, 1200);
    histMan.addH1D_BJetBinned("pt_hadTop", "pt_hadTop", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("pt_lepTop", "pt_lepTop", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("leadingJetPt", "leadingJetPt", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("nextLeadingJetPt", "nextLeadingJetPt", 6000, 0,
    	6000);
    histMan.addH1D_BJetBinned("trailingJetPt", "trailingJetPt", 6000, 0, 6000);

/*    histMan.addH1D_BJetBinned("mAllTop", "mAllTop", 500, 0, 500);
    histMan.addH1D_BJetBinned("m3", "m3", 5000, 0, 5000);

    histMan.addH1D_BJetBinned("mLepTopLoneNM", "mLepTopLoneNM", 600, 0, 600);
*/
    histMan.addH1D_LepBinned("mLepTopLone", "mLepTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mLepTopLoneneg", "mLepTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mLepTopLonepos", "mLepTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mHadTopLone", "mHadTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mHadTopLoneneg", "mLepTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mHadTopLonepos", "mHadTopLone", 600, 0, 600);
    histMan.addH1D_LepBinned("mass4jets", "mass4jets", 6000, 0, 6000);
    histMan.addH1D_LepBinned("numGenJets", "numGenJets", 20, 0, 20);
    histMan.addH1D_LepBinned("missGenJets", "missing GenJets for Reco", 20, 0,
    	20);
    histMan.addH1D_LepBinned("missRecoJets", "missing reco jets for GenJets",
    	20, 0, 20);
    histMan.addH1D_LepBinned("numDiffJets", "numDiffJets", 20, -10, 10);
    histMan.addH1D_BJetBinned("ttbar_pt", "ttbar_pt", 1200, 0, 1200);
/*    histMan.addH1D_BJetBinned("mHadTopLoneNM", "mHadTopLoneNM", 600, 0, 600);
    histMan.addH1D_BJetBinned("mAllTop", "mAllTop", 600, 0, 600);
    histMan.addH1D_BJetBinned("chiHadronicTop", "chiHadronicTop", 400, 0, 200);
    histMan.addH1D_BJetBinned("chiLeptonicTop", "chiLeptonicTop", 400, 0, 200);
    histMan.addH1D_BJetBinned("chiGlobal", "chiGlobal", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("chiTotal", "chiTotal", 6000, 0, 6000);
    */

    histMan.addH1D_BJetBinned("lepTopDeltaR", "lepTopDeltaR", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("hadTopDeltaR", "hadTopDeltaR", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("lepTopLep", "lepTopLep", 3, 0, 3);
    histMan.addH1D_BJetBinned("hadTopHad", "hadTopHad", 3, 0, 3);
    histMan.addH1D_BJetBinned("numHadTopMCMatches", "numHadTopMCMatches", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopMCMatches", "numLepTopMCMatches", 6, 0, 6);
    
    histMan.addH1D_LepBinned("numLepTopDauMCMatches", "numLepTopDauMCMatches", 6, 0, 6);
    histMan.addH1D_LepBinned("numHadTopDauMCMatches", "numHadTopDauMCMatches", 6, 0, 6);
    
    histMan.addH1D_BJetBinned("numHadTopCorrectID", "numHadTopCorrectID", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopCorrectID", "numLepTopCorrectID", 6, 0, 6);
		
    histMan.addH1D_BJetBinned("numHadTopMCMatches2top", "numHadTopMCMatches2top", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopMCMatches2top", "numLepTopMCMatches2top", 6, 0, 6);
    
    histMan.addH1D_LepBinned("numLepTopDauMCMatches2top", "numLepTopDauMCMatches2top", 6, 0, 6);
    histMan.addH1D_LepBinned("numHadTopDauMCMatches2top", "numHadTopDauMCMatches2top", 6, 0, 6);
    
    histMan.addH1D_BJetBinned("numHadTopCorrectID2top", "numHadTopCorrectID2top", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopCorrectID2top", "numLepTopCorrectID2top", 6, 0, 6);
    /*
    histMan.addH1D_BJetBinned("lepTopDeltaRNM", "lepTopDeltaRNM", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("hadTopDeltaRNM", "hadTopDeltaRNM", 900, 0, 3.1);
    histMan.addH1D_BJetBinned("numHadTopMCMatchesNM", "numHadTopMCMatchesNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopMCMatchesNM", "numLepTopMCMatchesNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numHadTopCorrectIDNM", "numHadTopCorrectIDNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("numLepTopCorrectIDNM", "numLepTopCorrectIDNM", 6, 0, 6);
    histMan.addH1D_BJetBinned("tPrime_pt", "tPrime_pt", 6000, 0, 6000);
    histMan.addH2D_BJetBinned("tPrime_pt_vs_tPrimeM", "tPrime_pt_vs_tPrimeM", 600, 0, 6000, 600, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_px", "tPrime_px", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_py", "tPrime_py", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("tPrime_pz", "tPrime_pz", 6000, 0, 6000);

    histMan.addH1D_BJetBinned("ttbar_pt_withMETCut", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_withMETAndAsymJets", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_withAsymJetsCut", "ttbar_pt", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withMETCut", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withMETCut", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withMETCut", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withMETAndAsymJets", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withMETAndAsymJets", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withMETAndAsymJets", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH1D_BJetBinned("ttbar_pt_2ndSolution_withAsymJetsCut", "ttbar_pt_2ndSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_3rdSolution_withAsymJetsCut", "ttbar_pt_3rdSolution", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_allSolutions_withAsymJetsCut", "ttbar_pt_allSolutions", 1000, 0, 1000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution", "ttbar_pt_vs_mttbar_2ndSolution", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution", "ttbar_pt_vs_mttbar_3rdSolution", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions", "ttbar_pt_vs_mttbar_allSolutions", 500, 0, 5000, 500,
            0, 5000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500,
            0, 5000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000, 500, 0,
            5000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withMETAndAsymJets", "ttbar_pt_vs_mttbar", 500, 0,
            5000, 100, 0, 1000);

    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_2ndSolution_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_3rdSolution_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
    histMan.addH2D_BJetBinned("ttbar_pt_vs_mttbar_allSolutions_withAsymJetsCut", "ttbar_pt_vs_mttbar", 500, 0, 5000,
            100, 0, 1000);
//
    histMan.addH1D_BJetBinned("ttbar_px", "ttbar_px", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_py", "ttbar_py", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pz", "ttbar_pz", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("ttbar_pt_QCDEnriched", "ttbar_pt", 1000, 0, 1000);
*/
    histMan.addH1D_LepBinned("HT_ltop", "HT_ltop", 6000, 0, 6000);
    histMan.addH1D_LepBinned("MET_ltop", "MET_ltop", 6000, 0, 6000);
    histMan.addH1D_LepBinned("METpt_ltop", "METpt_ltop", 6000, 0, 6000);
    histMan.addH1D_LepBinned("ST_ltop", "ST_ltop", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("HT3jet_ltop", "HT3jet_ltop", 6000, 0, 6000);
    histMan.addH1D_BJetBinned("HT4jet_ltop", "HT4jet_ltop", 6000, 0, 6000);
/*
    histMan.addH1D_BJetBinned("HT", "HT", 5000, 0, 5000);
    histMan.addH2D_BJetBinned("HTvsMttbar", "HT vs mttbar", 500, 0, 5000, 500, 0, 5000);
*/    histMan.addH1D("numberOfJets", "numberOfJets", 10, 0, 10);
    histMan.addH1D("numberOfBJets", "numberOfBJets", 10, 0, 10);
/*    histMan.addH1D_BJetBinned("MET", "MET", 200, 0, 1000);
    histMan.addH2D_BJetBinned("METvsMttbar", "MET vs mttbar", 500, 0, 5000, 200, 0, 1000);
    histMan.addH1D_BJetBinned("leadingJetMass", "leadingJetMass", 200, 0, 200);
    histMan.addH1D_BJetBinned("leadingGenJetMass", "leadingGenJetMass", 200, 0, 200);
    histMan.addH1D_BJetBinned("mtW", "mtW", 600, 0, 600);
*/
		const std::string mtprefixes[] = {"mt4j", "mt2j"};
		const std::string mtmidnames[] = {"bothchg", "plus", "minus"};
		const std::string mtlepsuffixes[] = {"lep", "elec", "mu"};
		const std::string mtjetsuffixes[] = {"4j", "5j"};
		TString mthistname;
		for (unsigned int pcnt = 0; pcnt < 2; ++pcnt) {
			mthistname = mtprefixes[pcnt];
			for (unsigned int mcnt = 0; mcnt < 3; ++mcnt) {
				TString mthistname1 = mthistname;
				mthistname1 += mtmidnames[mcnt];
				for (unsigned int lcnt = 0; lcnt < 3; ++lcnt) {
					TString mthistname2 = mthistname1;
					mthistname2 += mtlepsuffixes[lcnt];
					for (unsigned int jcnt = 0; jcnt < 2; ++jcnt) {
						TString mthistname3 = mthistname2;
						mthistname3 += mtjetsuffixes[jcnt];
						histMan.addH1D_BJetBinned((const char *) mthistname3, (const char *) mthistname3, 6000.0, 0.0, 6000.0);
					}
				}
			}
		}
    histMan.addH1D("electronD0", "electronD0", 1000, 0, 0.2);
    histMan.addH1D_BJetBinned("neutrino_pz", "neutrino_pz", 1000, -500, 500);
    histMan.addH1D_BJetBinned("lone_neutrino_pz", "lone_neutrino_pz", 1000, -600, 600);
/*    histMan.addH2D("ptRel_vs_DRmin", "ptRel_vs_DRmin", 100, 0, 1, 300, 0, 300);
    histMan.addH1D("ptRel_QCDenriched", "ptRel_QCDenriched", 300, 0, 300);
    histMan.addH1D("DRmin_QCDenriched", "DRmin_QCDenriched", 100, 0, 1);
    histMan.addH1D("ptRel_WZenriched", "ptRel_WZenriched", 300, 0, 300);
    histMan.addH1D("DRmin_WZenriched", "DRmin_WZenriched", 100, 0, 1);
    histMan.addH1D_JetBinned("QCDest_CombRelIso", "RelIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion", "RelIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_1btag", "RelIso (>=1 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion_1btag", "RelIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_2btag", "RelIso (>=2 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_CombRelIso_controlRegion_2btag", "RelIso control region", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithMETCut", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithMETCutAndAsymJetCuts", "PFIso", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_WithAsymJetCuts", "PFIso", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_1btag", "PFIso (>=1 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_1btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_2btag", "PFIso (>=2 btag)", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_2btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithMETCut", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithMETCutAndAsymJetCuts", "PFIso control region",
            1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion_WithAsymJetCuts", "PFIso control region", 1000, 0, 10);

    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCut", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithMETCutAndAsymJetCuts", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_WithAsymJetCuts", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_1btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_JetBinned("QCDest_PFIsolation_controlRegion2_2btag", "PFIso control region", 1000, 0, 10);
    histMan.addH1D_BJetBinned("pt_leadingTop", "pt_leadingTop", 1000, 0, 1000);
*/    histMan.addH1D_BJetBinned("pt_loneLepTop", "pt_loneLepTop", 1000, 0, 1000);
/*    histMan.addH1D_BJetBinned("pt_loneHadTop", "pt_loneHadTop", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop", "pt_NextToLeadingTop", 1000, 0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar", "pt_leadingTop_vs_mttbar", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar", "pt_NextToLeadingTop_vs_mttbar", 500, 0, 5000, 100, 0,
            1000);
    histMan.addH1D_BJetBinned("pt_leadingTop_withMETCut", "pt_leadingTop_withMETCut", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withMETCut", "pt_NextToLeadingTop_withMETCut", 1000, 0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETCut", "pt_leadingTop_vs_mttbar_withMETCut", 500, 0,
            5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETCut", "pt_NextToLeadingTop_vs_mttbar_withMETCut",
            500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("pt_leadingTop_withMETAndAsymJets", "pt_leadingTop_withMETAndAsymJets", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withMETAndAsymJets", "pt_NextToLeadingTop_withMETAndAsymJets", 1000,
            0, 1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withMETAndAsymJets",
            "pt_leadingTop_vs_mttbar_withMETAndAsymJets", 500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets",
            "pt_NextToLeadingTop_vs_mttbar_withMETAndAsymJets", 500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("pt_leadingTop_withAsymJetsCut", "pt_leadingTop_withMETAndAsymJets", 1000, 0, 1000);
    histMan.addH1D_BJetBinned("pt_NextToLeadingTop_withAsymJetsCut", "pt_NextToLeadingTop_withMETAndAsymJets", 1000, 0,
            1000);
    histMan.addH2D_BJetBinned("pt_leadingTop_vs_mttbar_withAsymJetsCut", "pt_leadingTop_vs_mttbar_withAsymJetsCut",
            500, 0, 5000, 100, 0, 1000);
    histMan.addH2D_BJetBinned("pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut",
            "pt_NextToLeadingTop_vs_mttbar_withAsymJetsCut", 500, 0, 5000, 100, 0, 1000);

    histMan.addH1D_BJetBinned("angleTops", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withMETCut", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withMETAndAsymJets", "angle between top quarks", 400, 0, 4);
    histMan.addH1D_BJetBinned("angleTops_withAsymJetsCut", "angle between top quarks", 400, 0, 4);

    histMan.addH2D_BJetBinned("angleTops_vs_mttbar", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withMETCut", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withMETAndAsymJets", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
    histMan.addH2D_BJetBinned("angleTops_vs_mttbar_withAsymJetsCut", "angleTops_vs_mttbar", 500, 0, 5000, 400, 0, 4);
*/
    histMan.addH1D_BJetBinned("angleHadTopD", "angle between hadronic top and light jet", 400, 0, 4);
    histMan.addH2D_BJetBinned("mass_vs_transMass", "mass_vs_transMass", 500, 0, 1000, 500, 0, 1000);
    histMan.addH2D_BJetBinned("good_vs_badwp1topmass", "good_vs_badwp1topmass", 500, 0, 1300, 500, 0, 1300);
}

void Analysis::initMicroNtuple() {
    microTuple->Branch("type",   &type,   "type/I");
    microTuple->Branch("weight", &weight, "weight/D");
    microTuple->Branch("PileUpWeight", &PileUpWeight, "PileUpWeight/D");
    microTuple->Branch("eventNumber",   &eventNumber,   "eventNumber/I");
    microTuple->Branch("runNumber",     &runNumber,     "runNumber/I");
    microTuple->Branch("numberOfJets",  &numberOfJets,  "numberOfJets/I");
    microTuple->Branch("numberOfBJets", &numberOfBJets, "numberOfBJets/I");
    microTuple->Branch("leptoCharge",   &leptoCharge,   "leptoCharge/I");
    microTuple->Branch("eSEL",          &eSEL,          "eSEL/I");
    microTuple->Branch("muSEL",         &muSEL,         "muSEL/I");
    microTuple->Branch("ST",            &ST,            "ST/D");
    microTuple->Branch("HT",            &HT,            "HT/D");
    microTuple->Branch("MET",           &MET,           "MET/D");
    microTuple->Branch("nPileUpVtx",           &nPileUpVtx,           "nPileUpVtx/I");
    microTuple->Branch("leadingJetPdgId", &leadingJetPdgId, "leadingJetPdgId/I");
    microTuple->Branch("leadingJetIndGen",&leadingJetIndGen,"leadingJetIndGen/I");
    microTuple->Branch("leadingJetPtGen", &leadingJetPtGen, "leadingJetPtGen/D");
    microTuple->Branch("leadingJetEtaGen",&leadingJetEtaGen,"leadingJetEtaGen/D");
    microTuple->Branch("leadingJetPhiGen",&leadingJetPhiGen,"leadingJetPhiGen/D");
    microTuple->Branch("leadingJetPtRec", &leadingJetPtRec, "leadingJetPtRec/D");
    microTuple->Branch("leadingJetEtaRec",&leadingJetEtaRec,"leadingJetEtaRec/D");
    microTuple->Branch("leadingJetPhiRec",&leadingJetPhiRec,"leadingJetPhiRec/D");

    microTuple->Branch("JetPx[numberOfJets]",JetPx);
    microTuple->Branch("JetPy[numberOfJets]",JetPy);
    microTuple->Branch("JetPz[numberOfJets]",JetPz);
    microTuple->Branch("JetE[numberOfJets]",JetE);
    microTuple->Branch("JetBTag[numberOfJets]",JetBTag);
   
    microTuple->Branch("BJetPx[numberOfBJets]",BJetPx);
    microTuple->Branch("BJetPy[numberOfBJets]",BJetPy);
    microTuple->Branch("BJetPz[numberOfBJets]",BJetPz);
    microTuple->Branch("BJetE[numberOfBJets]",BJetE);


    microTuple->Branch("freeJetPdgId",    &freeJetPdgId,   "freeJetPdgId/I");
    microTuple->Branch("freeJetIndGen",   &freeJetIndGen,  "freeJetIndGen/I");
    microTuple->Branch("freeJetPtGen",    &freeJetPtGen,   "freeJetPtGen/D");
    microTuple->Branch("freeJetEtaGen",   &freeJetEtaGen,  "freeJetEtaGen/D");
    microTuple->Branch("freeJetPhiGen",   &freeJetPhiGen,  "freeJetPhiGen/D");
    microTuple->Branch("freeJetPtRec",    &freeJetPtRec,   "freeJetPtRec/D");
    microTuple->Branch("freeJetEtaRec",   &freeJetEtaRec,  "freeJetEtaRec/D");
    microTuple->Branch("freeJetPhiRec",   &freeJetPhiRec,  "freeJetPhiRec/D");

    microTuple->Branch("bJetTlPtGen",    &bJetTlPtGen,  "bJetTlPtGen/D");
    microTuple->Branch("bJetTlEtaGen",   &bJetTlEtaGen, "bJetTlEtaGen/D");
    microTuple->Branch("bJetTlPhiGen",   &bJetTlPhiGen, "bJetTlPhiGen/D");
    microTuple->Branch("bJetTlPtRec",    &bJetTlPtRec,  "bJetTlPtRec/D");
    microTuple->Branch("bJetTlEtaRec",   &bJetTlEtaRec, "bJetTlEtaRec/D");
    microTuple->Branch("bJetTlPhiRec",   &bJetTlPhiRec, "bJetTlPhiRec/D");

    microTuple->Branch("bJetThPtGen",    &bJetThPtGen,  "bJetThPtGen/D");
    microTuple->Branch("bJetThEtaGen",   &bJetThEtaGen, "bJetThEtaGen/D");
    microTuple->Branch("bJetThPhiGen",   &bJetThPhiGen, "bJetThPhiGen/D");
    microTuple->Branch("bJetThPtRec",    &bJetThPtRec,  "bJetThPtRec/D");
    microTuple->Branch("bJetThEtaRec",   &bJetThEtaRec, "bJetThEtaRec/D");
    microTuple->Branch("bJetThPhiRec",   &bJetThPhiRec, "bJetThPhiRec/D");

    microTuple->Branch("jet1WhPdgId",    &jet1WhPdgId,  "jet1WhPdgId/I");
    microTuple->Branch("jet2WhPdgId",    &jet2WhPdgId,  "jet2WhPdgId/I");
    microTuple->Branch("jet1WhIndGen",   &jet1WhIndGen, "jet1WhIndGen/I");
    microTuple->Branch("jet2WhIndGen",   &jet2WhIndGen, "jet2WhIndGen/I");
    microTuple->Branch("jet1WhPtGen",    &jet1WhPtGen,  "jet1WhPtGen/D");
    microTuple->Branch("jet1WhEtaGen",   &jet1WhEtaGen, "jet1WhEtaGen/D");
    microTuple->Branch("jet1WhPhiGen",   &jet1WhPhiGen, "jet1WhPhiGen/D");
    microTuple->Branch("jet1WhPtRec",    &jet1WhPtRec,  "jet1WhPtRec/D");
    microTuple->Branch("jet1WhEtaRec",   &jet1WhEtaRec, "jet1WhEtaRec/D");
    microTuple->Branch("jet1WhPhiRec",   &jet1WhPhiRec, "jet1WhPhiRec/D");
    microTuple->Branch("jet2WhPtGen",    &jet2WhPtGen,  "jet2WhPtGen/D");
    microTuple->Branch("jet2WhEtaGen",   &jet2WhEtaGen, "jet2WhEtaGen/D");
    microTuple->Branch("jet2WhPhiGen",   &jet2WhPhiGen, "jet2WhPhiGen/D");
    microTuple->Branch("jet2WhPtRec",    &jet2WhPtRec,  "jet2WhPtRec/D");
    microTuple->Branch("jet2WhEtaRec",   &jet2WhEtaRec, "jet2WhEtaRec/D");
    microTuple->Branch("jet2WhPhiRec",   &jet2WhPhiRec, "jet2WhPhiRec/D");

    microTuple->Branch("metPtGen",       &metPtGen,     "metPtGen/D");
    microTuple->Branch("metPhiGen",      &metPhiGen,    "metPhiGen/D");
    microTuple->Branch("metEtaGen",      &metEtaGen,    "metEtaGen/D");
    microTuple->Branch("metPtRec",       &metPtRec,     "metPtRec/D");
    microTuple->Branch("metPhiRec",      &metPhiRec,    "metPhiRec/D");
    microTuple->Branch("metEtaRec",      &metEtaRec,    "metEtaRec/D");

    microTuple->Branch("leptonPdgId",    &leptonPdgId,  "leptonPdgId/I");
    microTuple->Branch("leptonPtGen",    &leptonPtGen,  "leptonPtGen/D");
    microTuple->Branch("leptonEtaGen",   &leptonEtaGen, "leptonEtaGen/D");
    microTuple->Branch("leptonPhiGen",   &leptonPhiGen, "leptonPhiGen/D");
    microTuple->Branch("leptonPtRec",    &leptonPtRec,  "leptonPtRec/D");
    microTuple->Branch("leptonEtaRec",   &leptonEtaRec, "leptonEtaRec/D");
    microTuple->Branch("leptonPhiRec",   &leptonPhiRec, "leptonPhiRec/D");

    microTuple->Branch("wlPdgId",        &wlPdgId,      "wlPdgId/I");
    microTuple->Branch("wlPtGen",        &wlPtGen,      "wlPtGen/D");
    microTuple->Branch("wlEtaGen",       &wlEtaGen,     "wlEtaGen/D");
    microTuple->Branch("wlPhiGen",       &wlPhiGen,     "wlPhiGen/D");
    microTuple->Branch("wlMgen",         &wlMgen,       "wlMgen/D");
    microTuple->Branch("wlPtRec",        &wlPtRec,      "wlPtRec/D");
    microTuple->Branch("wlEtaRec",       &wlEtaRec,     "wlEtaRec/D");
    microTuple->Branch("wlPhiRec",       &wlPhiRec,     "wlPhiRec/D");
    microTuple->Branch("wlMrec",         &wlMrec,       "wlMrec/D");

    microTuple->Branch("tlPdgId",        &tlPdgId,      "tlPdgId/I");
    microTuple->Branch("tlPtGen",        &tlPtGen,      "tlPtGen/D");
    microTuple->Branch("tlEtaGen",       &tlEtaGen,     "tlEtaGen/D");
    microTuple->Branch("tlPhiGen",       &tlPhiGen,     "tlPhiGen/D");
    microTuple->Branch("tlMgen",         &tlMgen,       "tlMgen/D");
    microTuple->Branch("tlPtRec",        &tlPtRec,      "tlPtRec/D");
    microTuple->Branch("tlEtaRec",       &tlEtaRec,     "tlEtaRec/D");
    microTuple->Branch("tlPhiRec",       &tlPhiRec,     "tlPhiRec/D");
    microTuple->Branch("tlMrec",         &tlMrec,       "tlMrec/D");

    microTuple->Branch("whPdgId",        &whPdgId,      "whPdgId/I");
    microTuple->Branch("whPtGen",        &whPtGen,      "whPtGen/D");
    microTuple->Branch("whEtaGen",       &whEtaGen,     "whEtaGen/D");
    microTuple->Branch("whPhiGen",       &whPhiGen,     "whPhiGen/D");
    microTuple->Branch("whMgen",         &whMgen,       "whMgen/D");
    microTuple->Branch("whPtRec",        &whPtRec,      "whPtRec/D");
    microTuple->Branch("whEtaRec",       &whEtaRec,     "whEtaRec/D");
    microTuple->Branch("whPhiRec",       &whPhiRec,     "whPhiRec/D");
    microTuple->Branch("whMrec",         &whMrec,       "whMrec/D");

    microTuple->Branch("thPdgId",        &thPdgId,      "thPdgId/I");
    microTuple->Branch("thPtGen",        &thPtGen,      "thPtGen/D");
    microTuple->Branch("thEtaGen",       &thEtaGen,     "thEtaGen/D");
    microTuple->Branch("thPhiGen",       &thPhiGen,     "thPhiGen/D");
    microTuple->Branch("thMgen",         &thMgen,       "thMgen/D");
    microTuple->Branch("thPtRec",        &thPtRec,      "thPtRec/D");
    microTuple->Branch("thEtaRec",       &thEtaRec,     "thEtaRec/D");
    microTuple->Branch("thPhiRec",       &thPhiRec,     "thPhiRec/D");
    microTuple->Branch("thMrec",         &thMrec,       "thMrec/D");

    microTuple->Branch("wpPdgId",        &wpPdgId,      "wpPdgId/I");
    microTuple->Branch("wpPtGen",        &wpPtGen,      "wpPtGen/D");
    microTuple->Branch("wpEtaGen",       &wpEtaGen,     "wpEtaGen/D");
    microTuple->Branch("wpPhiGen",       &wpPhiGen,     "wpPhiGen/D");

    microTuple->Branch("dRtlBjet",&dRtlBjet,"dRtlBjet/D");
    microTuple->Branch("dRthBjet",&dRthBjet,"dRthBjet/D");
    microTuple->Branch("dRwhJet1",&dRwhJet1,"dRwhJet1/D");
    microTuple->Branch("dRwhJet2",&dRwhJet2,"dRwhJet2/D");

    microTuple->Branch("tPosJetMassGen", &tPosJetMassGen, "tPosJetMassGen/D");
    microTuple->Branch("tPosJetPtGen",   &tPosJetPtGen,   "tPosJetPtGen/D");
    microTuple->Branch("tPosJetEtaGen",  &tPosJetEtaGen,  "tPosJetEtaGen/D");
    microTuple->Branch("tPosJetPhiGen",  &tPosJetPhiGen,  "tPosJetPhiGen/D");

    microTuple->Branch("tNegJetMassGen", &tNegJetMassGen, "tNegJetMassGen/D");
    microTuple->Branch("tNegJetPtGen",   &tNegJetPtGen,   "tNegJetPtGen/D");
    microTuple->Branch("tNegJetEtaGen",  &tNegJetEtaGen,  "tNegJetEtaGen/D");
    microTuple->Branch("tNegJetPhiGen",  &tNegJetPhiGen,  "tNegJetPhiGen/D");

    microTuple->Branch("tPosJetMassRec", &tPosJetMassRec, "tPosJetMassRec/D");
    microTuple->Branch("tPosJetPtRec",   &tPosJetPtRec,   "tPosJetPtRec/D");
    microTuple->Branch("tPosJetEtaRec",  &tPosJetEtaRec,  "tPosJetEtaRec/D");
    microTuple->Branch("tPosJetPhiRec",  &tPosJetPhiRec,  "tPosJetPhiRec/D");

    microTuple->Branch("tNegJetMassRec", &tNegJetMassRec, "tNegJetMassRec/D");
    microTuple->Branch("tNegJetPtRec",   &tNegJetPtRec,   "tNegJetPtRec/D");
    microTuple->Branch("tNegJetEtaRec",  &tNegJetEtaRec,  "tNegJetEtaRec/D");
    microTuple->Branch("tNegJetPhiRec",  &tNegJetPhiRec,  "tNegJetPhiRec/D");

//    microTuple->Branch("matchRtlJet",  &matchRtlJet,  "matchRtlJet/D");
//    microTuple->Branch("matchRthJet",  &matchRthJet,  "matchRthJet/D");
//    microTuple->Branch("matchRwhJet1", &matchRwhJet1, "matchRwhJet1/D");
//    microTuple->Branch("matchRwhJet2", &matchRwhJet2, "matchRwhJet2/D");
}

Analysis::Analysis() :
    eventReader(new NTupleEventReader()),
//    eventFilter(Filter::makeTopPairEPlusJetsFilter()),
    currentEvent(),
    ttbarCandidate(),
    histMan(),
    cutflow(),
    singleCuts(),
    cutflowPerFile(),
    singleCutsPerFile(),
    mucutflow(),
    musingleCuts(),
    mucutflowPerFile(),
    musingleCutsPerFile(),
    interestingEvents(),
    brokenEvents(),
    eventCheck(),
    weights(Analysis::luminosity/*current lumi*/, 7, "Cert_160404_180252_7TeV_Collisions11_JSON.pileup_v2.root"),
    weight(0),
    cutflowPerSample(DataType::NUMBER_OF_DATA_TYPES, toplikeElectronSelSize,
                    JetBin::NUMBER_OF_JET_BINS),
    mucutflowPerSample(DataType::NUMBER_OF_DATA_TYPES, toplikeMuonSelSize,
                    JetBin::NUMBER_OF_JET_BINS),
		wpTopPos(0),
		wpTopNeg(0),
		specTopPos(0),
		specTopNeg(0)
{
    //    outputfile->SetCompressionLevel(7);
    for (unsigned int cut = 0; cut < toplikeElectronSelSize; ++cut) {
        cutflow[cut] = 0;
        singleCuts[cut] = 0;
    }
    for (unsigned int cut = 0; cut < toplikeMuonSelSize; ++cut) {
        mucutflow[cut] = 0;
        musingleCuts[cut] = 0;
    }
}

Analysis::~Analysis() {
    histMan.writeToDisk();
    outputFile->cd();
    microTuple->Write();
    outputFile->Close();
}

void Analysis::addInputFile(const char* fileName) {
    eventReader->addInputFile(fileName);
}

void Analysis::setMaximalNumberOfEvents(long maxEvents) {
    if (maxEvents > 0) {
        eventReader->setMaximumNumberOfEvents(maxEvents);
    }
}

void Analysis::setUsedNeutrinoSelectionForTopPairReconstruction(NeutrinoSelectionCriterion::value selection) {
    TopPairEventCandidate::usedNeutrinoSelection = selection;
}

void Analysis::setUsedTTbarReconstructionCriterion(TTbarReconstructionCriterion::value selection) {
	TopPairEventCandidate::usedTTbarReconstruction = selection;
}

void Analysis::checkForDuplicatedEvents(){
    map<unsigned long, std::vector<unsigned long> >::const_iterator iter;
    std::vector<pair<unsigned long, unsigned long> > duplicateEvents;

    for(iter = eventCheck.begin(); iter != eventCheck.end(); ++iter){
        std::vector<unsigned long> events = (*iter).second;
        std::sort(events.begin(), events.end());
        for(unsigned long ev = 0; ev < events.size() -1; ++ev){
            if(events.at(ev) == events.at(ev +1)){
                duplicateEvents.push_back(make_pair((*iter).first, events.at(ev)));
            }
        }
    }

    if (duplicateEvents.size() > 0){
        cout << "found duplicate events" << endl;
        for(unsigned long ev = 0; ev < duplicateEvents.size() -1; ++ev){
            cout << "run: " << duplicateEvents.at(ev).first << " event: " << duplicateEvents.at(ev).second << endl;
        }
    }
}

void Analysis::checkForBrokenEvents(){
    if(ttbarCandidate.Electrons().size() == 0){
        brokenEvents.push_back(InterestingEvent(ttbarCandidate, eventReader->getCurrentFile()));
    }

    if(ttbarCandidate.eventnumber() == 1019245){
        cout << "broken event" << endl;
        ttbarCandidate.inspect();
    }
}
