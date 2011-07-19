////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Mar  7 16:07:36 EST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>


bool cmp (std::pair<float, float> A, std::pair<float, float> B)
{
  return (A.first < B.first);
}


int DataMe (std::vector<std::string> const& Files)
{
  float A,B;
  std::vector<std::pair<float, float> > C;
  for (size_t i = 0; i != Files.size(); ++i) {
    std::ifstream In(Files[i].c_str());
    if (!In) {
      std::cerr << "ERROR: cannot open file" << std::endl;
      throw;
    }

    In >> A >> B;
    C.push_back( std::make_pair<float, float>(A,B) );
  }

  std::sort(C.begin(), C.end(), cmp);

  for (size_t i = 0; i != C.size(); ++i) {
    printf(" %10.3f", C[i].first);
  }
  printf("\n");
  for (size_t i = 0; i != C.size(); ++i) {
    printf(" %10E", C[i].second);
  }
  printf("\n");


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [Files]" << std::endl;
    return 1;
  }

  std::vector<std::string> Files;
  for (int i = 1; i < argc; ++i) {
    Files.push_back(argv[i]);
  }

  DataMe(Files);

  return 0;
}
