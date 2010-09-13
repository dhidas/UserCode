////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Sep 13 14:25:22 CDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <vector>
#include <fstream>


int CombineLHEFiles (std::string OutFileName, std::vector<std::string> InFileNames)
{
  std::ofstream OutFile(OutFileName.c_str());
  if (!OutFile.is_open()) {
    std::cerr << "ERROR: cannot open output file " << OutFileName << std::endl;
    return 1;
  }

  for (size_t ifile = 0; ifile != InFileNames.size(); ++ifile) {
    std::cout << "File: " << ifile << "  " << InFileNames[ifile] << std::endl;
    std::ifstream InFile(InFileNames[ifile].c_str());
    if (!InFile.is_open()) {
      std::cerr << "ERROR: cannot open input file " << InFileNames[ifile]<< std::endl;
      return 1;
    }

    // If you haven't gotten the header, get it now or skip it
    if (ifile == 0) {
      for (std::string oneline; std::getline(InFile, oneline); ) {
        if (oneline.find("</LesHouchesEvents>", 0) == std::string::npos) {
          OutFile << oneline << "\n";
        }
      }
    } else {
      for (std::string oneline; std::getline(InFile, oneline); ) {
        if (oneline.find("<event>", 0) != std::string::npos) {
          OutFile << oneline << "\n";
          for ( ; std::getline(InFile, oneline); ) {
            if (oneline.find("</event>") != std::string::npos) {
              OutFile << oneline << "\n";
              break;
            } else {
              OutFile << oneline << "\n";
            }
          }
        }
      }
    }

    InFile.close();
  }

  OutFile << "</LesHouchesEvents>";
  OutFile.close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName] [InFileName]s" << std::endl;
    return 1;
  }

  std::string const OutFileName = argv[1];

  std::vector<std::string> InFileNames;
  for (int i = 2; i < argc; ++i) {
    InFileNames.push_back(argv[i]);
  }

  return CombineLHEFiles(OutFileName, InFileNames);

  return 0;
}
