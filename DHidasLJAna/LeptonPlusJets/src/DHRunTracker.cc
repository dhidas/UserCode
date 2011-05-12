#include "DHidasLJAna/LeptonPlusJets/interface/DHRunTracker.h"

DHRunTracker::DHRunTracker ()
{
}


DHRunTracker::~DHRunTracker ()
{
}


bool DHRunTracker::IsDuplicate(unsigned int run)
{
  if (std::binary_search(fVector.begin(), fVector.end(), run)) {
    return true;
  }

  fVector.push_back(run);
  std::sort(fVector.begin(), fVector.end());
  return false;
}




