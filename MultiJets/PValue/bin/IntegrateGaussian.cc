////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Mar 31 19:28:20 CDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TMath.h"
#include "TF1.h"


int IntegrateGaussian (float const PVal)
{
  // Define Min and max
  float const Min = -40;
  float const Max =  40;

  // This function is normalized
  TF1 f("gaus", "TMath::Gaus(x, 0, 1, 1)", Min, Max);


  // Calculate the sigma from -inf (-min)
  float Sigma = 0;
  for ( ;  f.Integral(Min,  Sigma) < 1.0 - PVal; Sigma += 0.01) {
  }
  float const FromInfinity = Sigma;

  // Calculate the sigma from the erf
  float const FromCenter   = TMath::Sqrt(2)*TMath::ErfcInverse(PVal);

  // Print
  printf("PValues=%8E  -inf=%5.2f   mid=%5.3f\n", PVal, FromInfinity, FromCenter);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [p-value]" << std::endl;
    return 1;
  }

  float const PVal = atof(argv[1]);

  IntegrateGaussian(PVal);

  return 0;
}
