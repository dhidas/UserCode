#include "DHidasLJAna/LeptonPlusJets/interface/DHRunTracker.h"

DHRunTracker::DHRunTracker ()
{
}


DHRunTracker::~DHRunTracker ()
{
}


bool DHRunTracker::IsDuplicate(unsigned int run, unsigned int event)
{
  static std::pair<unsigned int, unsigned int> pair;
  pair = std::make_pair<unsigned int, unsigned int>(run, event);
  if (fSet.find(pair) == fSet.end()) {
    fSet.insert(pair);
    return false;
  }

  return true;
}




