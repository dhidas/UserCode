#ifndef GUARD_PLTEvent_h
#define GUARD_PLTEvent_h

#include "PLTBinaryFileReader.h"
#include "PLTTelescope.h"
#include "PLTGainCal.h"


#include <map>


class PLTEvent
{
  public:
    PLTEvent ();
    PLTEvent (std::string const);
    PLTEvent (std::string const, std::string const);
    ~PLTEvent ();


    std::vector<PLTPlane*> fPlanes;
    std::vector<PLTTelescope*> fTelescopes;
    std::vector<PLTHit> fHits;


    size_t NPlanes ();
    PLTPlane* Plane(size_t);

    size_t NTelescopes ();
    PLTTelescope* Telescope (size_t);

    void Clear ();
    void AddHit (PLTHit);
    void MakeEvent ();
    unsigned long EventNumber ()
    { 
      return fEvent;
    }
    size_t NHits ()
    {
      return fHits.size();
    }

    int GetNextEvent ();

  private:
    unsigned long fRun;
    unsigned long fRunSection;
    unsigned long fEvent;

    PLTGainCal fGainCal;
    PLTBinaryFileReader fBinFile;

    std::map<int, PLTTelescope> fTelescopeMap;
    std::map<std::pair<int, int>, PLTPlane> fPlaneMap;

};



#endif
