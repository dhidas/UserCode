/*
 * RecoCandidateBase.h
 *
 *  Base class for reconstructed candidates
 *
 *	To create a new reconstructed candidate, called "mycandidate", for example,
 *	do the following in a new header file:

#include "RecoCandidateBase.h"


namespace OAT {

class mycandidate: public RecoCandidateBase {

public:
    mycandidate()
		{}

    mycandidate(const Event& event) :
			RecoCandidateBase(event)
		{}

    virtual bool passesFullSelection() const;
		virtual void reconstructCandidate(const ElectronPointer &electron);
		
		// The following is an example member function that would return the mass
		// of the candidate.  You could create other functions instead.
		virtual double getMass() const;
		
		
		// Add other functions here that return the values you want to plot.

protected:

//  Your data members,if any, should be declared here.

};
}

Then you have to code the last three functions in a new file you create in the
"src" directory.  Make sure to add the following lines to the top of your
source code file:

#include "../interface/RecoCandidateBase.h"

using namespace OAT;

After completing the code for "mycandidate", you can put calls to your class
in IntroAnalysis.h and IntroAnalysis.cpp.  Add the declaration of your
candidate to the "private" section of the IntroAnalysis class where other
candidates are declared.  Then, in IntroAnalysis.cpp, instantiate your
candidate in the initiateEvent member function where the other candidates
are instantiated.

Then you should write your own analysis function that will
check if the current event should be selected by calling passesFullSelection,
and, if so, reconstruct it with reconstructCandidate and then add its mass to
a histogram with getMass.  See doTTBarAnalysis in IntroAnalysis.cpp for an
example of an analysis function.  You need to add a call to your analysis
function in the main event loop of the "analyze" function.

 */
 

#ifndef RECOCANDIDATEBASE_H_
#define RECOCANDIDATEBASE_H_

#include <boost/shared_ptr.hpp>

#include "Event.h"
#include "RecoObjects/Particle.h"
#include "RecoObjects/Electron.h"
#include "RecoObjects/Jet.h"


using namespace BAT;

namespace OAT {

class RecoCandidateBase: public Event {

public:

    RecoCandidateBase()
		{}

    RecoCandidateBase(const Event& event) :
			Event(event)
		{}

    virtual ~RecoCandidateBase()
		{}

    virtual bool passesFullSelection() const = 0;

		/*
		Sample reconstruction function:

		virtual void reconstructCandidate(const ElectronPointer &electron);

		Note your reconstruction function may not take a parameter or may
		take different parameters from this example.
		*/

};

}

#endif /* RECOCANDIDATEBASE_H_ */
