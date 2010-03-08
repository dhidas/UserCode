#include "DHidasAna/Dtuple/interface/TDtupleReader.h"


TDtupleReader::TDtupleReader (TChain* Chain) : TDtuple(Chain)
{
}


TDtupleReader::~TDtupleReader ()
{
}


void TDtupleReader::Loop (long unsigned int Max)
{

  if (Max == 0) {
    for (long unsigned int ientry = 0; GetEntry(ientry); ++ientry) {
      if (ientry % 1000 == 0) {
        std::cout << "Events Processed: " << ientry << std::endl;
      }
      Analyze(ientry);
    }
  } else {
    for (long unsigned int ientry = 0; GetEntry(ientry) && (ientry < Max); ++ientry) {
      if (ientry % 1000 == 0) {
        std::cout << "Events Processed: " << ientry << std::endl;
      }
      Analyze(ientry);
    }
  }

  return;
}



