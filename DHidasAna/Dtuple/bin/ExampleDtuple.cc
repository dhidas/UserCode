#include <iostream>
#include <vector>

#include "DHidasAna/Dtuple/interface/TDtuple.h"

#include "TRandom.h"


//
// This program will make a small sample dtuple filled with random numbers
//

int main (int argc, char* argv[])
{
  // Get a Dtuple object
  TDtuple dtuple("ExampleDtuple");

  // If you want to see multiple output files
  // uncomment the line below
  //dtuple.SetMaxTreeSize(1000000);

  // make vectors for the TLeptons and TJets
  std::vector<TLepton> Leptons;
  std::vector<TJet> Jets;
  std::vector<TPhoton> Photons;

  
  // Loop over many events
  for (int i=0; i != 5000; ++i) {

    // Just informational output
    if (i % 1000 == 0) {
      std::cout << "Events Processed: " << i << std::endl;
    }


    // Clear the TLeptons, TZVertex, and TFakeCEM vectors.
    // This resets their size to zero too.
    Leptons.clear();
    Jets.clear();
    Photons.clear();


    // Loop to fill a few TLepton objects
    for (int j=0; j != 3; ++j) {

      // Create a TLepton object
      TLepton lep;

      // Fill all variables in the new TLepton object
      lep.SetCharge(gRandom->Gaus(0,1) > 0 ? 1 : -1);
      lep.SetTrkIso(gRandom->Rndm(2));
      lep.SetPx(gRandom->Gaus(0,50));
      lep.SetPy(gRandom->Gaus(0,50));
      lep.SetPz(gRandom->Gaus(0,50));
      lep.SetE(lep.P()+0.1);
      lep.SetTrkPt(gRandom->Gaus(0,50));
      lep.SetEmE((lep.Perp()+0.0001)*2/3);
      lep.SetHadE((lep.Perp()+0.0001)/3);
      lep.SetDetEta(gRandom->Gaus(0,2));
      lep.SetZ0(gRandom->Gaus(0,20));

      // Add the new TLepton object to our vector of TLepton objects
      Leptons.push_back(lep);

    }


    // Loop to fill a few TJet objects
    for (int j=0; j != 5; ++j) {

      // Create a TJet object
      TJet jet;

      // Fill all variables in the new TJet object
      jet.SetPx(gRandom->Gaus(0,50));
      jet.SetPy(gRandom->Gaus(0,50));
      jet.SetPz(gRandom->Gaus(0,50));
      jet.SetE(jet.P() + 5);
      jet.SetDetEta(gRandom->Gaus(0,2));
      jet.SetEScaleFactor(TMath::Abs(gRandom->Gaus(0,0.1)));
      jet.SetEmF(TMath::Abs((gRandom->Gaus(0,1))));
      jet.SetHadF(TMath::Abs((gRandom->Gaus(0,1))));

      // Add the new TJet object to our vector of TJet objects
      Jets.push_back(jet);
    }

    for (int j=0; j != 2; ++j) {

      // Create a TJet object
      TPhoton photon;

      // Add the new TJet object to our vector of TJet objects
      Photons.push_back(photon);
    }






    // Before we add anything to the Dtuple object we must clear it.
    // This resets counters, clears objects, and sets some dummy values.
    // This must be done for each event(once only per event)
    dtuple.Clear();

    // Add the TLeptons, TJets, TZVertices, and TFakeCEMs to the Dtuple object
    dtuple.AddLeptons(Leptons);
    dtuple.AddJets(Jets);
    dtuple.AddPhotons(Photons);

    // Fill the event variables in the Dtuple object
    dtuple.SetMetX(100.0 * gRandom->Rndm(3));
    dtuple.SetMetY(100.0 * gRandom->Rndm(4));
    dtuple.SetRawMetX(100.0 * gRandom->Rndm(5));
    dtuple.SetRawMetX(100.0 * gRandom->Rndm(5));
    dtuple.SetRun(i);
    dtuple.SetEvent((int) (100000.0 * gRandom->Rndm(6)));
    dtuple.SetSumEt(100.0 * gRandom->Rndm(7));
    dtuple.SetLum(gRandom->Rndm(8));

    // Add this event to the Dtuple object tree.
    // If you don't call Fill this event won't be saved
    dtuple.Fill();

  }

  // You MUST call Write() if you want to save the information!!
  dtuple.Write();

  return 0;
}


