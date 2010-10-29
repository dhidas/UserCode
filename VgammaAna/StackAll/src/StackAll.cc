#include "StackAll.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "TKey.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFrame.h"

StackAll::StackAll (TString const InputFileName, TString const OutFileName)
{
  ReadInputFile(InputFileName);
  OpenOutFile(OutFileName);
}




StackAll::~StackAll ()
{
  CloseOutFile();
}




void StackAll::ReadInputFile (TString const InFileName)
{
  // Open file
  std::ifstream InFile(InFileName.Data());
  if (!InFile.is_open()) {
    std::cerr << "ERROR: cannot open input file; " << InFileName << std::endl;
    exit (1);
  }


  // Loop over lines in the input file.  If the line begins with a % skip it
  // or if it is blank skip it.  Don't start with a space or tab
  for ( std::string OneLine; std::getline(InFile, OneLine); ) {
    if (OneLine == "" || OneLine.at(0) == '%' || OneLine.at(0) == ' ' || OneLine.at(0) == '\t') {
      continue;
    }

    std::cout << "ReadInputFile: " << OneLine << std::endl;
    std::istringstream LineStream;
    LineStream.str(OneLine);

    FileProp* ThisFile = new FileProp;
    LineStream >> ThisFile->Name
               >> ThisFile->Color
               >> ThisFile->IsData
               >> ThisFile->Scale
               >> ThisFile->FileName;

    ThisFile->File       = new TFile(ThisFile->FileName, "read");


    // Check file is open
    if (!ThisFile->File->IsOpen()) {
      std::cerr << "ERROR: cannot open file: " << ThisFile->FileName << std::endl;
      exit(1);
    }

    fFiles.push_back(ThisFile);
    fPropMap.insert( std::pair<TString, FileProp*>(ThisFile->Name, ThisFile) );
  }

  return;
}




void StackAll::OpenOutFile (TString const Name)
{
  fOutFile = new TFile(Name, "recreate");
  if (!fOutFile->IsOpen()) {
    std::cerr << "ERROR: cannot open output file: " << Name << std::endl;
    exit(1);
  }

  return;
}




void StackAll::CloseOutFile ()
{
  fOutFile->Write();
  fOutFile->Close();
  return;
}




void StackAll::Run ()
{
  for (size_t i = 0; i != fFiles.size(); ++i) {
    this->ReadDirectory( (TDirectoryFile*) fFiles[i]->File, fFiles[i]);
  }

  this->MakeAllStacks();

  return;
}




void StackAll::MakeAllStacks ()
{
  for (std::map< TString, std::map<TString, std::map<TString, std::map<TString, void *> > > >::iterator iDir = fHistMap.begin(); iDir != fHistMap.end(); ++iDir) {
    TString const Dir = iDir->first;
    for (std::map< TString, std::map<TString, std::map<TString, void *> > >::iterator iType = iDir->second.begin(); iType != iDir->second.end(); ++iType) {
      TString const Type = iType->first;
      for (std::map< TString, std::map<TString, void *> >::iterator iName = iType->second.begin(); iName != iType->second.end(); ++iName) {
        TString const Name = iName->first;

        // Make the Legend
        TLegend Legend(0.6, 0.7, 0.9, 0.9, "");
        Legend.SetNColumns(2);
        Legend.SetFillColor(0);

        // For summing the stack
        double StackSum = 0.0;

        for (std::map< TString, void *>::iterator iProc = iName->second.begin(); iProc != iName->second.end(); ++iProc) {
          TString const Proc = iProc->first;

          printf("Stacking: %s  %s  %s  %s\n", Dir.Data(), Type.Data(), Name.Data(), Proc.Data());

          // Make a new Stack
          if (iProc == iName->second.begin()) {
            fStackMap[Dir][Type][Name] = new THStack(Name, Name);
          }

          // Fill the stack!

          FileProp* FP = fPropMap.find(Proc)->second;
          std::cout << "IsData: " << FP->IsData << "  Proc: " << Proc << std::endl;

          // This is because it works...
          ((TH1*)   iProc->second)->SetLineColor( FP->Color );
          ((TH1*)   iProc->second)->SetFillColor( FP->Color );
          std::cout << "Entries: " << ((TH1*)   iProc->second)->GetEntries() << std::endl;

          // Let's do some scaling if need be
          if (FP->Scale > 0) {
            ((TH1*)   iProc->second)->Scale(FP->Scale);
          } else if (FP->Scale < 0) {
            ((TH1*)   iProc->second)->Scale(-FP->Scale/ ((TH1*)   iProc->second)->Integral() );
          } else {
            // Do nothing!  ie zero means no scaling
          }
          std::cout << "Entries: " << ((TH1*)   iProc->second)->GetEntries() << std::endl;

          double ThisHistIntegral = 0.0;
          if (FP->IsData == 0) {
            Legend.AddEntry((TH1*) iProc->second, Proc, "f");
            if (Type.BeginsWith("TH1F")) {
              fStackMap[Dir][Type][Name]->Add( (TH1F*) iProc->second, "hist");
              ThisHistIntegral = ((TH1F*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH1D")) {
              fStackMap[Dir][Type][Name]->Add( (TH1D*) iProc->second, "hist");
              ThisHistIntegral = ((TH1D*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH2F")) {
              fStackMap[Dir][Type][Name]->Add( (TH2F*) iProc->second, "hist");
              ThisHistIntegral = ((TH2F*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH2D")) {
              fStackMap[Dir][Type][Name]->Add( (TH2D*) iProc->second, "hist");
              ThisHistIntegral = ((TH2D*) iProc->second)->Integral();
            }
            printf("Integral: %15s %12.4f   %s\n", Proc.Data(), ThisHistIntegral, (Dir+"/"+Name).Data());
            StackSum += ThisHistIntegral;
          }

        }
        printf("Integral: %15s %12.4f   %s\n", "StackSum", StackSum, (Dir+"/"+Name).Data());

        GetOrMakeDir(Dir)->cd();
        TCanvas Canvas(Name, Name);
        Canvas.cd();
        fStackMap[Dir][Type][Name]->Draw();

        for (std::map< TString, void *>::iterator iProc = iName->second.begin(); iProc != iName->second.end(); ++iProc) {
          TString const Proc = iProc->first;
          FileProp* FP = fPropMap.find(Proc)->second;

          // Set the marker size
          ( (TH1*) iProc->second)->SetMarkerStyle(8);
          ( (TH1*) iProc->second)->SetMarkerSize(0.7);
          ( (TH1*) iProc->second)->SetLineWidth(2.0);

          double ThisHistIntegral = 0.0;
          if (FP->IsData == 1) {
            Legend.AddEntry((TH1*) iProc->second, Proc, "lp");
            if (Type.BeginsWith("TH1F")) {
              ( (TH1F*) iProc->second)->Draw("epsame");
              ThisHistIntegral = ((TH1F*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH1D")) {
              ( (TH2D*) iProc->second)->Draw("epsame");
              ThisHistIntegral = ((TH1D*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH2F")) {
              ( (TH1F*) iProc->second)->Draw("epsame");
              ThisHistIntegral = ((TH2F*) iProc->second)->Integral();
            } else if (Type.BeginsWith("TH2D")) {
              ( (TH2D*) iProc->second)->Draw("epsame");
              ThisHistIntegral = ((TH2D*) iProc->second)->Integral();
            }
            printf("Integral: %15s %12.4f   %s\n", Proc.Data(), ThisHistIntegral, (Dir+"/"+Name).Data());
          }

        }

        Legend.Draw("same");
        Canvas.SetBorderMode(0);
        Canvas.SetHighLightColor(0);
        Canvas.SetFillColor(0);
        //Canvas.GetPad(0)->GetFrame()->SetBorderMode(0);
        Canvas.Write();
        Canvas.SaveAs(Name+".eps");
        //fStackMap[Dir][Type][Name]->Write();
      }
    }
  }


  return;
}




void StackAll::ReadDirectory (TDirectoryFile* Dir, FileProp const* Prop)
{
  TIter MyIter(Dir->GetListOfKeys());
  TKey* Key;
  while ( (Key = (TKey*) MyIter()) ) {
    TObject *Object = (TObject*) Key->ReadObj();
    std::cout << "ClassName: " << Object->ClassName()
              << "     Name: " << Object->GetName() << std::endl;

    if (Object == 0x0) {
      continue;
    }

    TString const MyClassName = Object->ClassName();
    TString const PathName = Dir->GetPath();
    TString const DirName = PathName(PathName.First(":") + 1, PathName.Length() - PathName.First(":") - 1);

    if (MyClassName.BeginsWith("TH1F")) {
      AddToHistMap( (TH1F*) Object, DirName, Prop);
    } else if (MyClassName.BeginsWith("TH1D")) {
      AddToHistMap( (TH1D*) Object, DirName, Prop);
    } else if (MyClassName.BeginsWith("TH2F")) {
      AddToHistMap( (TH2F*) Object, DirName, Prop);
    } else if (MyClassName.BeginsWith("TH2D")) {
      AddToHistMap( (TH2D*) Object, DirName, Prop);
    }

    if (TString(Object->ClassName()).BeginsWith("TDirectory")) {
      this->ReadDirectory( (TDirectoryFile*) Object, Prop);
    }
  }

  return;
}




template <typename T>
void StackAll::AddToHistMap (T* Hist, TString const& Dir, FileProp const* Prop)
{
  TString const Proc = Prop->Name;
  std::cout << "Proc: " << Proc << std::endl;
  if (fHistMap[Dir][Hist->ClassName()][Hist->GetName()].find(Proc) == fHistMap[Dir][Hist->ClassName()][Hist->GetName()].end()) {
    fHistMap[Dir][Hist->ClassName()][Hist->GetName()][Proc] = (T*) Hist->Clone();
    ( (T*) fHistMap[Dir][Hist->ClassName()][Hist->GetName()][Proc])->SetDirectory(0x0);
  } else {
    ( (T*) fHistMap[Dir][Hist->ClassName()][Hist->GetName()][Proc])->Add(Hist);
  }

  return;
}




TDirectoryFile* StackAll::GetOrMakeDir (TString const& DirName)
{
  std::cout << "DirName: " << DirName << std::endl;

  TDirectoryFile* ThisDir = (TDirectoryFile*) fOutFile->GetDirectory(DirName);
  if (ThisDir) {
    return ThisDir;
  }

  TDirectoryFile* BackDir = GetOrMakeDir( DirName(0, DirName.Last('/')));
  if (BackDir) {
    return (TDirectoryFile*) BackDir->mkdir( DirName(DirName.Last('/')+1, DirName.Length()-DirName.Last('/')-1).Data() );
  }

  return GetOrMakeDir( DirName(0, DirName.Last('/')) );
}

