/*
 * EventPrinter.h
 *
 *  Created on: 22 Aug 2010
 *      Author: kreczko
 */

#ifndef EVENTTABLEPRINTER_H_
#define EVENTTABLEPRINTER_H_
#include "../EventCounter.h"
#include "../EventWeightProvider.h"
#include "../ToplikeCandidate.h"

namespace BAT {

namespace EventTablePrinter {
void printCutFlowCSVTable(const Counter& cutFlow);
void printCutFlowLatexTable(const Counter& cutFlow);
void printUnweightedCutFlowLatexTable(const Counter& cutFlow);
void printWprimeCutFlowUnwt(const Counter& cutFlow, EventWeightProvider &weights,
	const ToplikeSelectionSteps::Step cutArray[], unsigned int cutArraySiz,
	const char *const selType);
void printWprimeCutFlow(const Counter& cutFlow, EventWeightProvider &weights,
	const ToplikeSelectionSteps::Step cutArray[], unsigned int cutArraySiz,
	const char *const selType, const float luminosity = 5.0);
}
}
#endif /* EVENTPRINTER_H_ */
