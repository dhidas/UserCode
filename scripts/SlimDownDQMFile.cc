////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hida@cern.ch>
//
// Created on: Thu Nov  5 14:49:53 CET 2009
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TObjString.h"

int CopyObject(TFile& InFile, TFile& OutFile, TString const ObjName)
{
  OutFile.cd();
  TObjArray* ObjArray = ObjName.Tokenize("/");
  for (int i = 0; i != ObjArray->GetEntries(); ++i) {
    TString const NewDirName = ( (TObjString*) ObjArray->At(i) )->String();
    // check this dir first...TODO
    if (!OutFile.mkdir(NewDirName)) {
      std::cerr << "ERROR: cannot make dir in output file" << std::endl;
      return 1;
    }
    OutFile.cd(NewDirName);
  }
  //InFile.Get(ObjName)->Clone()->Write();
  OutFile.cd();
  return 0;
}


int CopyDirectory(TFile& InFile, TFile& OutFile, TString const DirName)
{
  std::cout << "Looking at Directory: " << DirName << std::endl;
  // Get contents and loop over them
  TDirectory* Dir = InFile.GetDirectory(DirName);
  if (Dir == 0x0) {
    std::cerr << "ERROR: Directory does not exist: " << DirName << std::endl;
    return 1;
  }

  TIter MyIter(Dir->GetListOfKeys());
  TKey* Key;
  while ( (Key = (TKey*) MyIter()) ) {
    TObject* Obj = Key->ReadObj();
    TString const Name = Obj->GetName();
    if (TString(Obj->ClassName()).BeginsWith("TH")) {
      std::cout << "Trying to Copy: " << Name << std::endl;
      CopyObject(InFile, OutFile, Name);
    } else if (TString(Obj->ClassName()).Contains("Directory")) {
      CopyDirectory(InFile, OutFile, DirName+Name+"/");
    }
  }
  return 0;
}




int SlimDownDQMFile (TString const InFileName, TString const OutFileName)
{
  // Copy only the selected histograms from InFile to OutFile.  OutFile is
  // assumed to be a new file which does not already exist.  If the output
  // file already exists then I will not overwrite it.  This is just to be safe

  // Open the input file
  TFile InFile(InFileName, "read");
  if (!InFile.IsOpen()) {
    std::cerr << "ERROR: cannot open input file " << InFileName << std::endl;
    return 1;
  }

  // Open the output file
  TFile OutFile(OutFileName, "create");
  if (!OutFile.IsOpen()) {
    std::cerr << "ERROR: cannot open output file " << OutFileName << std::endl;
    return 1;
  }

  std::vector<TString> Names;
  Names.push_back("DQMData/");

  // Loop over all of the hist names
  for (std::vector<TString>::iterator NameIt = Names.begin(); NameIt != Names.end(); ++NameIt) {
    TString const ThisName = *NameIt;

    if (ThisName.EndsWith("/")) {
      CopyDirectory(InFile, OutFile, ThisName);
    } else {
      CopyObject(InFile, OutFile, ThisName);
    }
  }


  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [InFile] [OutFile]" << std::endl;
    return 1;
  }

  SlimDownDQMFile(argv[1], argv[2]);

  return 0;
}
