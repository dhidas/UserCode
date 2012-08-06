#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>



#include "pulseHeightSpectrum.h"


using namespace std;
//============================================================================
//
// Efficeny


int main(int argc, char **argv){


    //char filename[50];
    int runNr;

    if(argc < 2) {

        // in case non suffient commandline input was given
        cout<<"Usage: ./pulseHeightMain runNr \n";
        return 0;
    }

    runNr = atoi(argv[1]);

    // constructor defines measurement number
    pulseHeightSpectrum meas(runNr);
    // loads files (returns 0 in case of error)
    if (meas.loadFile() == 0) return 0;
    
    // 1st event loop: create the pixel mask with different error codes
    meas.createMask();

    // 2nd Eventloop fill create and fill clusters
    meas.clusterLoop();
    
    meas.fitLandau();

    //sprintf(filename,"eff_cce_hist_%d.root",runNr);
    meas.writeHistograms();

    return 1;
}

