////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Mar  1 15:32:23 CET 2012
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>

#include "PLTEvent.h"
#include "PLTU.h"

#include "TH1F.h"
#include "TH2F.h"



int FindNoisyPixels (std::string const DataFileName, std::string const AlignmentFileName)
{
  // Set some basic style
  PLTU::SetStyle();

  std::cout << "DataFileName:    " << DataFileName << std::endl;
  std::cout << "AlignmentFileName: " << AlignmentFileName << std::endl;

  // Grab the plt event reader
  PLTEvent Event(DataFileName, "", AlignmentFileName);

  PLTPlane::FiducialRegion MyFiducialRegion = PLTPlane::kFiducialRegion_Diamond;
  Event.SetPlaneClustering(PLTPlane::kClustering_AllTouching, MyFiducialRegion);
  Event.SetPlaneFiducialRegion(MyFiducialRegion);


  std::map<int, TH2F*> Map;

  // Loop over all events in file
  for (int ientry = 0; Event.GetNextEvent() >= 0; ++ientry) {
    if (ientry % 10000 == 0) {
      std::cout << "Processing event: " << ientry << std::endl;
    }

    for (size_t itele = 0; itele != Event.NTelescopes(); ++itele) {
      PLTTelescope* Telescope = Event.Telescope(itele);
      //std::cout << Telescope->NTracks() << std::endl;
      //Telescope->DrawTracksAndHits("hello.gif");
      //exit(0);
    }

    // Loop over all planes with hits in event
    for (size_t ip = 0; ip != Event.NPlanes(); ++ip) {

      // THIS plane is
      PLTPlane* Plane = Event.Plane(ip);

      int const ID = 10 * Plane->Channel() + Plane->ROC();

      if (!Map.count(ID)) {
        TString const Name = TString::Format("Occ_%i", ID);
        Map[ID] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW,PLTU::FIRSTROW, PLTU::LASTROW+1);
      }

      for (size_t ih = 0; ih != Plane->NHits(); ++ih) {
        PLTHit* Hit = Plane->Hit(ih);

        Map[ID]->Fill(Hit->Column(), Hit->Row());
      }
    }
  }

  for (std::map<int, TH2F*>::iterator It = Map.begin(); It != Map.end(); ++It) {
    int const Channel = It->first / 10;
    int const ROC     = It->first % 10;


    // Grab a 1-D hist from the 2D occupancy
    TH1F* Hist = PLTU::HistFrom2D(It->second, "", 200, true);
    

    Double_t Median[1];
    Double_t Prob[1] = { 0.5 };
    Hist->GetQuantiles(1, Median, Prob);
    float const Threshold = Median[0] + 4.0 * Hist->GetRMS();
    delete Hist;

    TH2F* Hist2D = It->second;

    for (int i = 1; i <= Hist2D->GetNbinsX(); ++i) {
      for (int j = 1; j <= Hist2D->GetNbinsY(); ++j) {
        if (Hist2D->GetBinContent(i, j) > Threshold) {
          printf("Channel: %2i ROC: %1i  Hot Pixel: %2i %2i  Hit: %15i    Median = %15.3E  Threshold = %15.3E\n",
              Channel, ROC, PLTU::FIRSTCOL + i, PLTU::FIRSTROW + j, (int) Hist2D->GetBinContent(i, j), Median[0], Threshold);
        }
      }
    }



    delete It->second;

  }


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [DataFileName] [AlignmentFileName]" << std::endl;
    return 1;
  }

  std::string const DataFileName = argv[1];
  std::string const AlignmentFileName = argv[2];
  FindNoisyPixels(DataFileName, AlignmentFileName);

  return 0;
}
