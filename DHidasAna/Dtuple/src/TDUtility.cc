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




float TDUtility::GetConversionR (TLorentzVector& trk1_p4, 
							      int trk1_q, float trk1_d0, 
                    TLorentzVector& trk2_p4, 
							      int trk2_q, float trk2_d0, 
							      float bFieldAtOrigin) {


  std::pair<float, float> xy = GetConversionXY(trk1_p4, trk1_q, trk1_d0, trk2_p4, trk2_q, trk2_d0, bFieldAtOrigin);
  return sqrt(xy.first*xy.first + xy.second*xy.second);
}


std::pair<float, float> TDUtility::GetConversionXY(
                    TLorentzVector& trk1_p4, 
							      int trk1_q, float trk1_d0, 
                    TLorentzVector& trk2_p4, 
							      int trk2_q, float trk2_d0, 
							      float bFieldAtOrigin) {
  
  
  double tk1Curvature = -0.3*bFieldAtOrigin*(trk1_q/trk1_p4.Perp())/100.;
  double rTk1 = TMath::Abs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_p4.Phi());
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_p4.Phi());

  double tk2Curvature = -0.3*bFieldAtOrigin*(trk2_q/trk2_p4.Perp())/100.;
  double rTk2 = TMath::Abs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_p4.Phi());
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_p4.Phi());

  double dist = sqrt(pow(xTk2-xTk1, 2) + pow(yTk2-yTk1 , 2));
  dist = dist - (rTk1 + rTk2);


  double CosTheta = (xTk2-xTk1)/sqrt( pow(xTk2-xTk1, 2) + pow(yTk2-yTk1, 2) );
  double SinTheta = (yTk2-yTk1)/sqrt( pow(xTk2-xTk1, 2) + pow(yTk2-yTk1, 2) );
  double MyX = xTk1 + (rTk1 + 0.5*dist)*CosTheta;
  double MyY = yTk1 + (rTk1 + 0.5*dist)*SinTheta;



  //printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", xTk1, yTk1, xTk2, yTk2, CosTheta, SinTheta);
  return std::make_pair(MyX, MyY);
  
}


void TDUtility::PrintMapIntInt (std::map<int, int>& MyMap, TString const Name)
{
  for (std::map<int, int>::iterator It = MyMap.begin(); It != MyMap.end(); ++It) {
    printf("%s  %15i %15i\n", Name.Data(), It->first, It->second);
  }
  return;
}
