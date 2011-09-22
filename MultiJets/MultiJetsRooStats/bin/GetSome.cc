////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Thu Sep 22 11:37:34 CEST 2011
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


int GetSome (std::string const InFileName)
{
  std::ifstream f(InFileName.c_str());

  FILE* of = fopen("SomeData.dat", "w");

  for (std::string l; std::getline(f, l); ) {
    float v;
    std::istringstream ls;
    ls.str(l);
    for (int i = 0; ls >> v; ++i) {
      if (i % 2 == 0) {
        continue;
      }
      fprintf(of, " %E", v);
    }
    fprintf(of, "\n");



  }
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [InFileName]" << std::endl;
    return 1;
  }

  GetSome(argv[1]);

  return 0;
}
