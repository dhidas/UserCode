#include "DHidasLJAna/LeptonPlusJets/interface/DHRunTracker.h"

DHRunTracker::DHRunTracker ()
{
}


DHRunTracker::~DHRunTracker ()
{
}


bool DHRunTracker::IsDuplicate(unsigned int run, unsigned long event)
{
  static std::pair<unsigned int, unsigned long> pair;
  pair = std::make_pair<unsigned int, unsigned long>(run, event);
  if (fSet.count(pair)) {
    return true;
  }

  fSet.insert(pair);
  return false;
}




