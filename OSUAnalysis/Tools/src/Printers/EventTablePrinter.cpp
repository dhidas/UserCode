#include "../../interface/Printers/EventTablePrinter.h"
#include <iostream>
#include <string>
#include <iomanip>
#include "../../interface/TopPairEventCandidate.h"
#include "../../interface/ToplikeCandidate.h"
#include "../../interface/EventWeightProvider.h"

namespace BAT {
namespace EventTablePrinter {

using namespace std;

static const DataType::value sampleList[] = {DataType::DATA, DataType::WprimeTToTTD_M400,
			DataType::WprimeTToTTD_M600,
			DataType::WprimeTToTTD_M800, DataType::WprimeTToTTD_M1000,
			DataType::ttjets, // DataType::Wjets, // DataType::Zjets,
			DataType::QCD_BCtoE_Pt30to80,
			DataType::QCD_BCtoE_Pt80to170, DataType::QCD_EMEnriched_Pt80to170,
			DataType::QCD_MuEnriched_Pt20Pt15
			// DataType::QCD_Flat_15to3000
};	// Set numCols to give only 1 QCD column.

static const unsigned int listSiz = sizeof(sampleList) / sizeof(DataType::value);
static const unsigned int numCols = listSiz - 3; 	// 4 QCD samples


static string cnv_(string inp) {
	size_t loc_ = string::npos;
	while ((loc_ = inp.find('_')) != string::npos) {
		inp[loc_] = '-';
	}
	return (inp);
}


/*
static void stobByCut(const unsigned long nTotMCbyCut[], const Counter& cutFlow) {
    cout << "\\begin{sidewaystable} \n" << endl;
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "%       " << "S/B table (weighted)";
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "\\begin{tabular}{|l|rrr|}" << endl;
    cout << "\\hline" << endl;
    cout << "\\multicolumn{4}{|l|}";
    cout << "{S/B by Jet Cut}";
    cout << "\\\\\\hline" << endl;
    cout << " Signal mode " << setw(40) << left << right;
    //    cout << "Cut" << setw(40);
    
    for (unsigned int cut = ToplikeEPlusJetsSelection::AtLeastThreeGoodJets; cut < cutFlow.getSizeOfSecondDimension() &&
					cut < ToplikeEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
				if (cut == ToplikeEPlusJetsSelection::AtLeastOneBtag ||
						cut == ToplikeEPlusJetsSelection::AtLeastTwoBtags)
					continue;
        cout << " & " << setw(13) << ToplikeEPlusJetsSelection::StringSteps[cut];
		}
    cout << " \\\\\n\\hline" << endl;
		for (unsigned int type = DataType::WprimeTToTTD_M600; type <= DataType::WprimeToTBbar_M1000; ++type) {
				cout << cnv_(DataType::names.at(type));
				for (unsigned int cut = ToplikeEPlusJetsSelection::AtLeastThreeGoodJets; cut < cutFlow.getSizeOfSecondDimension() &&
						cut < ToplikeEPlusJetsSelection::NUMBER_OF_SELECTION_STEPS; ++cut) {
						if (cut == ToplikeEPlusJetsSelection::AtLeastOneBtag ||
								cut == ToplikeEPlusJetsSelection::AtLeastTwoBtags)
							continue;
						double totMC = nTotMCbyCut[cut];
						if (totMC <= 0)
							totMC = 1.0;
						cout << " & " << setw(12) << (double) cutFlow.sumThirdDimension(type, cut)/totMC;
				}
				cout << "  \\\\\n\\hline" << endl;
		}
    cout << "\\end{tabular}\\\\[5mm] \n" << endl;
    cout << "\\end{sidewaystable} \n" << endl;
    cout << endl;
}
*/


void printWprimeCutFlow(const Counter& cutFlow, EventWeightProvider &weights,
		const ToplikeSelectionSteps::Step cutArray[], unsigned int cutArraySiz,
		const char *const selType, const float luminosity)
{
    using namespace std;

    //print header
		static bool firstTime = true;
		if (firstTime) {
			cout << "Printing cut flow table\n" << endl;
			cout << "\\documentclass[legalpaper]{article} \n" << endl;
			cout << "\\usepackage{rotating} \n" << endl;
			cout << "\\begin{document} \n" << endl;
			firstTime = false;
		}
    cout << "\\begin{sidewaystable} \n" << endl;
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "%       " << selType << " selection cut flow table (weighted)";
    cout << "\n%------------------------------------------------------------------------\n";

    TString beginStmt =  "\\begin{tabular}{|l|r";
    for (unsigned int cnt = 1; cnt < numCols; ++cnt) {
			beginStmt += "r";
		}
		beginStmt += "|}";
    cout << beginStmt << endl;
    // cout << "\\begin{tabular}{|l|rrrrrrr|r|}" << endl;
    cout << "\\hline" << endl;
    TString multicol = "\\multicolumn{";
    multicol += (numCols + 1);
    multicol += "}{|l|}";
    cout << multicol;
    // cout << "\\multicolumn{9}{|l|}";
		float lumi = luminosity / 1000.0;
    cout << "{" << selType << " selection -- Expected number of events for ";
		cout << setw(2) << lumi << "/fb luminosity}";
    cout << "\\\\\\hline" << endl;
    cout << " Cut " << setw(40) << left << right;
    //    cout << "Cut" << setw(40);
    
    for (unsigned int type = 0; type < listSiz; ++type) {
        if (type == DataType::QCD_BCtoE_Pt30to80)  //first QCD data type
            cout << " &" << setw(13) << "QCD";
				else if (type < DataType::QCD_BCtoE_Pt20to30 || type > DataType::PhotonJets_Pt200toInf)
					cout << " & " << setw(13) << cnv_(DataType::names.at(sampleList[type]));
    }
    // cout << " & " << setw(25) << "Total Background MC \\\\\n\\hline" << endl;
    cout << "\\\\\n\\hline" << endl;
		unsigned long nTotMCbyCut[50];
    for (int cut = -5; cut < (int) cutFlow.getSizeOfSecondDimension() &&
				cut < (int) cutArraySiz; ++cut) {
				unsigned long nQCD = 0;
				unsigned long nTotalMC = 0;
				if (cut == -5)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Number generated" << right;
				else if (cut == -4)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "First skim step (PATtuple or generator filter)" << right;
				else if (cut == -3)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler number of electrons skimmed" << right;
				else if (cut == -2)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler number of muons skimmed" << right;
				else if (cut == -1)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler total number of events skimmed" << right;
				else cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< ToplikeSelectionSteps::StringSteps[cutArray[cut]] << right;
				for (unsigned int typeInd = 0; typeInd < listSiz; ++typeInd) {
						DataType::value type = sampleList[typeInd];
						double rawNum = -1.0;
						unsigned long numEvts = 0;
						switch (cut) {
							case -5:
									rawNum = weights.numberOfProcessedEvents[type];
									break;
							case -4:
									rawNum = weights.numberOfPattplSkimEvents[type];
									break;
							case -3:
									rawNum = weights.numberOfElectronSkimEvents[type];
									break;
							case -2:
									rawNum = weights.numberOfMuonSkimEvents[type];
									break;
							case -1:
									rawNum = weights.numberOfNtplSkimEvents[type];
									break;
								default:
									numEvts = cutFlow.sumThirdDimension(type, cut);
						}
						if (rawNum >= 0.0)
							numEvts = (rawNum * weights.getWeight(type)) + 0.5;
						if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf) {
									nQCD += numEvts;
									if (typeInd == (listSiz - 1) || type == DataType::PhotonJets_Pt200toInf)
										cout << " & " << setw(12) << fixed << nQCD;
						}
						if (type > DataType::DATA && type < DataType::WprimeTToTTD_M600)
								nTotalMC += numEvts;
							
						if (type >= DataType::QCD_BCtoE_Pt20to30 &&
									type <= DataType::PhotonJets_Pt200toInf)
								continue;

						cout << " & " << setw(12) << fixed << numEvts;
			}
			// cout << " & " << setw(12) << fixed << nTotalMC;
			if (cut >= 0)
				nTotMCbyCut[cut] = nTotalMC;
			cout << "  \\\\\n\\hline" << endl;
			/*
			
			cout << " Stage " << setw(2) << cut << " " << setw(40) << left
							<< ToplikeEPlusJetsSelection::StringSteps[cut] << right;
			for (unsigned int type = DataType::WprimeTToTTD_M600; type <= DataType::WprimeToTBbar_M1000; ++type) {
					double totMC = nTotalMC;
					if (totMC <= 0)
						totMC = 1.0;
					cout << " & " << setw(12) << (double) cutFlow.sumThirdDimension(type, cut)/totMC;
			}
			cout << " & & & & ";
			cout << "  \\\\\n\\hline" << endl;
			*/
    }
    cout << "\\end{tabular}\\\\[5mm] \n" << endl;
    cout << "\\end{sidewaystable} \n" << endl;
		// stobByCut(nTotMCbyCut, cutFlow);
    // cout << "\\end{document} \n" << endl;
    cout << endl;
}


void printWprimeCutFlowUnwt(const Counter& cutFlow, EventWeightProvider &weights,
		const ToplikeSelectionSteps::Step cutArray[], unsigned int cutArraySiz,
		const char *const selType)
{
    using namespace std;


    //print header
    cout << "\\begin{sidewaystable} \n" << endl;
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "%       " << selType << " selection cut flow table (unweighted)";
    cout << "\n%------------------------------------------------------------------------\n";

    TString beginStmt =  "\\begin{tabular}{|l|r";
    for (unsigned int cnt = 1; cnt < numCols; ++cnt) {
			beginStmt += "r";
		}
		beginStmt += "|}";
    cout << beginStmt << endl;
    // cout << "\\begin{tabular}{|l|rrrrrrr|r|}" << endl;
    cout << "\\hline" << endl;
    TString multicol = "\\multicolumn{";
    multicol += (numCols + 1);
    multicol += "}{|l|}";
    cout << multicol;
    // cout << "\\multicolumn{9}{|l|}";
    cout << "{" << selType << " selection -- Unscaled number of events}";
    cout << "\\\\\\hline" << endl;
    cout << " Cut " << setw(40) << left << right;
    //    cout << "Cut" << setw(40);
    
    for (unsigned int type = 0; type < listSiz; ++type) {
        if (type == DataType::QCD_BCtoE_Pt30to80)  //first QCD data type
            cout << " &" << setw(13) << "QCD";
				else if (type < DataType::QCD_BCtoE_Pt20to30 || type > DataType::PhotonJets_Pt200toInf)
					cout << " & " << setw(13) << cnv_(DataType::names.at(sampleList[type]));
    }
    // cout << " & " << setw(25) << "Total Background MC \\\\\n\\hline" << endl;
    cout << "\\\\\n\\hline" << endl;
    for (int cut = -5; cut < (int) cutFlow.getSizeOfSecondDimension() &&
				cut < (int) cutArraySiz; ++cut) {
				// if (cut == ToplikeEPlusJetsSelection::AtLeastOneBtag)
					// continue;
				unsigned long nQCD = 0;
				// unsigned long nSingleTop = 0;
				unsigned long nTotalMC = 0;
				if (cut == -5)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Number generated" << right;
				else if (cut == -4)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "First skim step (PATtuple or generator filter)" << right;
				else if (cut == -3)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler number of electrons skimmed" << right;
				else if (cut == -2)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler number of muons skimmed" << right;
				else if (cut == -1)
					cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< "Ntupler total number of events skimmed" << right;
				else cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< ToplikeSelectionSteps::StringSteps[cutArray[cut]] << right;
				for (unsigned int typeInd = 0; typeInd < listSiz; ++typeInd) {
						DataType::value type = sampleList[typeInd];
						unsigned long numEvts = 0;
						switch (cut) {
							case -5:
									numEvts = weights.numberOfProcessedEvents[type];
									break;
							case -4:
									numEvts = weights.numberOfPattplSkimEvents[type];
									break;
							case -3:
									numEvts = weights.numberOfElectronSkimEvents[type];
									break;
							case -2:
									numEvts = weights.numberOfMuonSkimEvents[type];
									break;
							case -1:
									numEvts = weights.numberOfNtplSkimEvents[type];
									break;
								default:
									numEvts = cutFlow.sumThirdDimensionUnweighted(type, cut);
						}
						if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf) {
									nQCD += numEvts;
									if (typeInd == (listSiz - 1) || type == DataType::PhotonJets_Pt200toInf)
										cout << " & " << setw(12) << fixed << nQCD;
						}
						if (type > DataType::DATA && type < DataType::WprimeTToTTD_M600)
								nTotalMC += numEvts;
							
						if (type >= DataType::QCD_BCtoE_Pt20to30 &&
									type <= DataType::PhotonJets_Pt200toInf)
								continue;

						cout << " & " << setw(12) << fixed << numEvts;
				}
				// cout << " & " << setw(12) << fixed << nTotalMC;
				cout << "  \\\\\n\\hline" << endl;
				
				/*
				cout << " Stage " << setw(2) << cut << " " << setw(40) << left
								<< ToplikeEPlusJetsSelection::StringSteps[cut] << right;
				for (unsigned int type = DataType::WprimeTToTTD_M600; type <= DataType::WprimeToTBbar_M1000; ++type) {
						double totMC = nTotalMC;
						if (totMC <= 0)
							totMC = 1.0;
						cout << " & " << setw(12) << (double) cutFlow.sumThirdDimensionUnweighted(type, cut)/totMC;
				}
				cout << " & & & & ";
				cout << "  \\\\\n\\hline" << endl;
				*/
    }

    cout << "\\end{tabular}\\\\[5mm] \n" << endl;
    cout << "\\end{sidewaystable} \n" << endl;
		static bool firstTime = true;
		if (firstTime)
			firstTime = false;
		else cout << "\\end{document} \n" << endl;
    cout << endl;
}


void printCutFlowLatexTable(const Counter& cutFlow) {
    using namespace std;
    cout << "Printing cut flow table" << endl;

    //print header
    cout << "\\newpage\n" << endl;
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "%       " << "cut flow table (weighted)";
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "\\begin{tabular}{|l|rrrrrr|r|}" << endl;
    cout << "\\hline" << endl;
    cout << "\\multicolumn{8}{|l|}";
    cout << "{Expected number of events for " << "36.145" << "/pb}";
    cout << "\\\\\\hline" << endl;
    cout << " Cut " << setw(40) << left << right;
    //    cout << "Cut" << setw(40);

    for (unsigned int type = 0; type < cutFlow.getSizeOfFirstDimension(); ++type) {
        if (type == DataType::QCD_BCtoE_Pt20to30) {//first QCD data type, this needs to be changed
            cout << " &" << setw(13) << "QCD";
        }
        if (type == DataType::singleTop_And_W)
            cout << " &" << setw(13) << "singeTop";

        if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
            continue;

        if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
            continue;
        if (type == DataType::VQQ || type >= DataType::Zprime_M500GeV_W5GeV)//don't use them now
            continue;
        cout << " & " << setw(13) << cnv_(DataType::names.at(type));
    }
    cout << " &" << setw(25) << "Total MC \\\\\n\\hline" << endl;
    for (unsigned int cut = 0; cut < cutFlow.getSizeOfSecondDimension(); ++cut) {
        unsigned long nQCD = 0;
        unsigned long nSingleTop = 0;
        unsigned long nTotalMC = 0;
        for (unsigned int type = 0; type < cutFlow.getSizeOfFirstDimension(); ++type) {
            if (type == 0) {
                cout << " Stage " << setw(2) << cut << " " << setw(40) << left
                        << TTbarEPlusJetsSelection::StringSteps[cut] << right;
            }
            if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
                nQCD += cutFlow.sumThirdDimension(type, cut);

            if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
                nSingleTop += cutFlow.sumThirdDimension(type, cut);

            if (type == DataType::PhotonJets_Pt200toInf)
                cout << " & " << setw(12) << fixed << nQCD;

            if (type == DataType::singleTopSChannel)
                cout << " & " << setw(12) << fixed << nSingleTop;

            if (type > DataType::DATA)
                nTotalMC += cutFlow.sumThirdDimension(type, cut);

            if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
                continue;
            if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
                continue;
            if (type == DataType::VQQ || type >= DataType::Zprime_M500GeV_W5GeV)//don't use them now
                continue;

            cout << " & " << setw(12) << fixed << cutFlow.sumThirdDimension(type, cut);
        }
        cout << " & " << setw(12) << fixed << nTotalMC;
        cout << "  \\\\\n\\hline" << endl;
    }

    cout << "\\end{tabular}\\\\[5mm]" << endl;
    cout << endl;
}


void printUnweightedCutFlowLatexTable(const Counter& cutFlow) {
    using namespace std;
    cout << "Printing cut flow table" << endl;

    //print header
    cout << "\\newpage\n" << endl;
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "%       " << "cut flow table (unweighted)";
    cout << "\n%------------------------------------------------------------------------\n";
    cout << "\\begin{tabular}{|l|rrrrrrr|r|}" << endl;
    cout << "\\hline" << endl;
    cout << "\\multicolumn{8}{|l|}";
    cout << "{Expected number of events for " << "36.145" << "/pb}";
    cout << "\\\\\\hline" << endl;
    cout << " Cut " << setw(40) << left << right;
    //    cout << "Cut" << setw(40);

    for (unsigned int type = 0; type < cutFlow.getSizeOfFirstDimension(); ++type) {
        if (type == DataType::QCD_BCtoE_Pt20to30) {//first QCD data type, this needs to be changed
            cout << " &" << setw(13) << "QCD";
        }
        if (type == DataType::singleTop_And_W)
            cout << " &" << setw(13) << "singeTop";

        if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
            continue;

        if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
            continue;
        if (type == DataType::VQQ || type >= DataType::Zprime_M500GeV_W5GeV)//don't use them now
            continue;
        cout << " & " << setw(13) << cnv_(DataType::names.at(type));
    }
    cout << " &" << setw(25) << "Total MC \\\\\n\\hline" << endl;
    for (unsigned int cut = 0; cut < cutFlow.getSizeOfSecondDimension(); ++cut) {
        unsigned long nQCD = 0;
        unsigned long nSingleTop = 0;
        unsigned long nTotalMC = 0;
        for (unsigned int type = 0; type < cutFlow.getSizeOfFirstDimension(); ++type) {
            if (type == 0) {
                cout << " Stage " << setw(2) << cut << " " << setw(40) << left
                        << TTbarEPlusJetsSelection::StringSteps[cut] << right;
            }
            if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
                nQCD += cutFlow.sumThirdDimensionUnweighted(type, cut);

            if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
                nSingleTop += cutFlow.sumThirdDimensionUnweighted(type, cut);

            if (type == DataType::PhotonJets_Pt200toInf)
                cout << " & " << setw(12) << fixed << nQCD;

            if (type == DataType::singleTopSChannel)
                cout << " & " << setw(12) << fixed << nSingleTop;

            if (type > DataType::DATA)
                nTotalMC += cutFlow.sumThirdDimensionUnweighted(type, cut);

            if (type >= DataType::singleTop_And_W && type <= DataType::singleTopSChannel)
                continue;
            if (type >= DataType::QCD_BCtoE_Pt20to30 && type <= DataType::PhotonJets_Pt200toInf)
                continue;
            if (type == DataType::VQQ || type >= DataType::Zprime_M500GeV_W5GeV)//don't use them now
                continue;

            cout << " & " << setw(12) << fixed << cutFlow.sumThirdDimensionUnweighted(type, cut);
        }
        cout << " & " << setw(12) << fixed << nTotalMC;
        cout << "  \\\\\n\\hline" << endl;
    }

    cout << "\\end{tabular}\\\\[5mm]" << endl;
    cout << endl;
}

}
}
