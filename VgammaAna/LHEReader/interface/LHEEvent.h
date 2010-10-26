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
    LHEEvent (std::vector<TString> const&);
    ~LHEEvent ();

    float Weight;
    std::vector<LHEParticle> Particles;


    int NextEvent ();

  private:
    bool OpenNextFile ();
    TString LHEFileName;
    std::ifstream* LHEFile;
    std::vector<TString> LHEFileNames;
    size_t ifile;

};



























#endif

