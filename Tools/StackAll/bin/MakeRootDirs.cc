////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Sat Sep 25 14:25:48 PDT 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"


TDirectoryFile* GetOrMakeDir (TString const& Path, TFile* InFile)
{
  TString const MyPath = Path(Path.First(":") + 1, Path.Length() - Path.First(":"));
  std::cout << "MyPath: " << MyPath << std::endl;

  TDirectoryFile* ThisDir = (TDirectoryFile*) InFile->GetDirectory(MyPath);
  if (ThisDir) {
    return (TDirectoryFile*) ThisDir;
  }

  TDirectoryFile* BackDir = GetOrMakeDir( MyPath(0, MyPath.Last('/')), InFile );
  if (BackDir) {
    return (TDirectoryFile*) BackDir->mkdir( MyPath(MyPath.Last('/')+1, MyPath.Length()-MyPath.Last('/')-1).Data() );
  }

  return GetOrMakeDir( Path(0, Path.Last('/')), InFile );
}





int MakeRootDirs ()
{
  TFile MyFile("TestRootFile.root", "recreate");

  TH1D* Hist = new TH1D("hist0", "hist0", 10, 0, 10);
  Hist->SetDirectory(&MyFile);
  Hist->FillRandom("gaus", 20);

  TDirectory* Dir = MyFile.mkdir("Dir1");
  Hist = new TH1D("hist1", "hist1", 10, 0, 10);
  Hist->SetDirectory(Dir);
  Hist->FillRandom("gaus", 20);

  Dir = Dir->mkdir("Dir1B");
  Hist = new TH1D("hist2", "hist2", 10, 0, 10);
  Hist->SetDirectory(Dir);
  Hist->FillRandom("gaus", 20);

  GetOrMakeDir("blah:/one/two/three/four", &MyFile);
  GetOrMakeDir("blah:/one/two/three/four/FIVE", &MyFile);

  Hist = new TH1D("hist3", "hist3", 10, 0, 10);
  Hist->SetDirectory( GetOrMakeDir("blah:/one/two/three/four/FIVE/Six", &MyFile) );
  Hist->FillRandom("gaus", 20);

  MyFile.Write();
  MyFile.Close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  MakeRootDirs();

  return 0;
}
