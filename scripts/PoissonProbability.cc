////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@fnal.gov>
//
// Created on: Fri Jun  5 08:29:37 CEST 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TMath.h"


int PoissonProbability (float const Expected, int const Observed)
{

  long double Sum = 0.0;
  for (int i = Observed; i < Observed + 200000; ++i) {
    Sum += TMath::Poisson(i, Expected);
  }
  printf("Probability of observing >= %i where %9.3f is expected is %3.20f\n", Observed, Expected, (double) Sum);

  Sum = 0.0;
  for (int i = Observed; i > Observed - 200000 && i >= 0; --i) {
    Sum += TMath::Poisson(i, Expected);
  }
  printf("Probability of observing <= %i where %9.3f is expected is %3.20f\n", Observed, Expected, (double) Sum);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [Expected] [Observed]" << std::endl;
    return 1;
  }

  PoissonProbability(atof(argv[1]), atoi(argv[2]));

  return 0;
}
