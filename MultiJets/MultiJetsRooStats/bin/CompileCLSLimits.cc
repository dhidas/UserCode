////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Sun Sep 18 07:57:40 EDT 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

class LimitObj
{
  public:
    LimitObj () {};
    ~LimitObj () {};

    float Mass;
    float Observed;
    float Expected[5];

    int operator< (const LimitObj& rhs) const
    {
      return this->Mass < rhs.Mass;
    }
};

int CompileCLSLimits (std::vector<std::string>& InFileNames)
{
  // Vector for holding all results
  std::vector<LimitObj> Limits;

  // Open each file
  for (size_t ifile = 0; ifile != InFileNames.size(); ++ifile) {

    LimitObj L;

    // Open file
    std::ifstream f(InFileNames[ifile].c_str());
    f >> L.Mass;
    if (f.eof()) {
      continue;
    }

    f >> L.Expected[0]
      >> L.Expected[1]
      >> L.Expected[2]
      >> L.Expected[3]
      >> L.Expected[4]
      >> L.Observed;

    Limits.push_back(L);
    if (f.eof()) {
      std::cerr << "ERROR: file is at end: " << InFileNames[ifile] << std::endl;
      throw;
    }
  }

  // Sort them
  std::sort(Limits.begin(), Limits.end());

  // Open file for output
  FILE* f = fopen("Limits.dat", "w");
  fprintf(f, "Test\n");

  for (std::vector<LimitObj>::iterator it = Limits.begin(); it != Limits.end(); ++it) {
    fprintf(f, " %12.3E", it->Mass);
  }
  fprintf(f, "\n");

  for (int ie = 0; ie != 5; ++ie) {
    for (std::vector<LimitObj>::iterator it = Limits.begin(); it != Limits.end(); ++it) {
      fprintf(f, " %12.3E", it->Expected[ie]);
    }
    fprintf(f, "\n");
  }

  for (std::vector<LimitObj>::iterator it = Limits.begin(); it != Limits.end(); ++it) {
    fprintf(f, " %12.3E", it->Observed);
  }
  fprintf(f, "\n");
  fprintf(f, "1");

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]s" << std::endl;
    return 1;
  }

  std::vector<std::string> InFileNames;
  for (int i = 1; i < argc; ++i) {
    InFileNames.push_back(argv[i]);
  }

  CompileCLSLimits(InFileNames);

  return 0;
}
