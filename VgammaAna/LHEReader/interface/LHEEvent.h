#ifndef GUARD_LHEEvent_h
#define GUARD_LHEEvent_h

#include <vector>
#include <fstream>

#include "TString.h"
#include "TLorentzVector.h"

class LHEParticle : public TLorentzVector
{
  public:
    LHEParticle () {};
    ~LHEParticle () {};

    int Id;
    int Mo1;
    int Mo2;
};


class LHEEvent
{
  public:
    LHEEvent (TString const);
    ~LHEEvent ();

    float Weight;
    std::vector<LHEParticle> Particles;


    int NextEvent ();

  private:
    TString LHEFileName;
    std::ifstream* LHEFile;

};



























#endif

