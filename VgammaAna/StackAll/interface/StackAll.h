#ifndef GUARD_StackAll_h
#define GUARD_StackAll_h

#include <vector>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "THStack.h"


class FileProp
{
  public:
    FileProp () {};
    ~FileProp () {};

    TString Name;
    int Color;
    int IsData;
    float Scale;
    TString FileName;

    TFile* File;
};


class StackAll
{
  public:
    StackAll (TString const, TString const);
    ~StackAll ();


    void ReadInputFile (TString const);
    void OpenOutFile (TString const);
    void CloseOutFile ();
    void Run ();
    void MakeAllStacks ();

  private:
    TFile* fOutFile;
    std::vector<FileProp*> fFiles;

    void ReadDirectory (TDirectoryFile*, FileProp const*);

    template <typename T> void AddToHistMap (T*, TString const&, FileProp const*);
    TDirectoryFile* GetOrMakeDir (TString const&);

    std::multimap<TString, FileProp*> fPropMap;
    std::map<TString, std::map<TString, std::map<TString, std::map<TString, void *> > > > fHistMap;
    std::map<TString, std::map<TString, std::map<TString, THStack*> > > fStackMap;
    //fHistMap[Dir][Type][Name][Proc];
    //fStackMap[Dir][Type][Name];

};





#endif
