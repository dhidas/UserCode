// Author: Dean Andrew Hidas <http://www-cdf.fnal.gov/~dhidas/>

////////////////////////////////////////////////////////////////
//
// TDUtility
//
// The TDUtility Class!
//
// This is a class of globally useful functions that are
// essentially analysis independant.
//
////////////////////////////////////////////////////////////////


#include "DHidasAna/Dtuple/interface/TDUtility.h"

// Root Likes ClassDef and ClassImp.
// Comment them out if you don't need them.
// There should NOT be a ; since this is a root macro and not a function
ClassImp(TDUtility)



TDUtility::TDUtility ()
{
  // Default constructor
}



TDUtility::~TDUtility ()
{
  // Destructor
}



bool TDUtility::IsDuplicateEvent (int const run, int const event)
{
  // Returns true if this is a duplicate event, false if the run-event 
  // pair is not already in the list.

  // Static set of run-event pairs
  static std::set< std::pair<int, int> > RunEventSet;

  // Make a pair of the input
  std::pair<int, int> ThisEvent(run, event);

  // If we find this pair in the set return true(duplicate found)
  if (RunEventSet.find(ThisEvent) != RunEventSet.end()) {
    return true;
  }

  // It's not a duplicate so save it to the set and return false
  RunEventSet.insert(ThisEvent);
  return false;
}




TString TDUtility::GetLeptonFlavorsString (std::vector<TLepton>& Leptons)
{
  TString Electrons, Muons, Taus;

  for (std::vector<TLepton>::iterator Lep = Leptons.begin(); Lep != Leptons.end(); ++Lep) {
    if (Lep->IsFlavor(TLepton::kLeptonFlavor_Electron)) {
      Electrons += Lep->GetFlavorString();
    } else if (Lep->IsFlavor(TLepton::kLeptonFlavor_Muon)) {
      Muons += Lep->GetFlavorString();
    } else if (Lep->IsFlavor(TLepton::kLeptonFlavor_Tau)) {
      Taus += Lep->GetFlavorString();
    }
  }

  return Electrons+Muons+Taus;
}
