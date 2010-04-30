// Author: Dean Andrew Hidas <http://www-cdf.fnal.gov/~dhidas/>

////////////////////////////////////////////////////////////////
//
// TDtuple
//
// The TDtuple class has been designed as a common analysis small
// ntuple (Dtuple!).  The TDtuple consists of TClonesArray(s) of
// TLepton objects, TJet objects, counters for each of those objects
// and some event variables which belong to the TDtuple class.
//
////////////////////////////////////////////////////////////////

#include "DHidasAna/Dtuple/interface/TDtuple.h"

#include "TSystem.h"





// Root Likes ClassDef and ClassImp.
// Comment them out if you don't need them.
// There should NOT be a ; since this is a root macro and not a function
ClassImp(TDtuple)



//
// Default constructor
//
TDtuple::TDtuple ()
{
  // Default constructor.
  // the default output file is set to dtuple.root

  SetRootFile("dtuple.root");
}




//
// Constructor: takes the name of the desired output root file
//
TDtuple::TDtuple (std::string const InputFileName)
{
  // Constructor: takes the name of the desired output root file
  // If the name does not end in .root it will be appended for you

  // Ntuple output name
  std::string NtupleName;

  // Does this name end in .root already or not.  if not let's make it.
  if (InputFileName.size() >= 5 && std::string(InputFileName.end()-5, InputFileName.end()) == ".root") {
    NtupleName = InputFileName;
  } else {
    NtupleName = InputFileName + ".root";
  }

  // Setup the output root file
  SetRootFile(NtupleName);
}




//
// Constructor: takes two strings, the dataset name and a string which you wish
//              to append to it.  usually a number: ie ttop0z_24.root where
//              ttop0z is the data set name and 24 is the number yuo want to append
//
TDtuple::TDtuple (std::string const InputDataSetName, std::string const InputFileNumber)
{
  // Constructor: takes to strings as inputs.  This is usually the dataset name
  // and a number in string format) you wish to append to it:
  // ttop0z_24.root where ttop0z is the dataset name and 24 is the
  // number you may want to append.  This is useful for submitting caf jobs

  // The output rootfile name
  std::string const NtupleName = InputDataSetName + "_" + InputFileNumber + ".root";

  // Setup the rootfile
  SetRootFile(NtupleName);

}


//
//// Constructor: Takes a TFile* which would be the output file assuming you have already
//                defined it.
TDtuple::TDtuple (TFile *fi)
{
  // Will take as input a pointer to a TFile.  This will setup the branches
  // in a rootfile that you have already defined.  This option is provided
  // because you may have already defined an output TFile.  This is a very good
  // constructor to use.

  // Setup the root output file with the given TFile*
  SetRootFile(fi);

}




//
// Constructor: Takes a TChain*.  This should be used if you want to read an existing Dtuple tree
// 
TDtuple::TDtuple (TChain* chain)
{
  // Constructor!  Takes a TChain* that you wish to read.  This is the constructor that you should
  //               use if you are planning on reading an existing ntuple which contains the dtuple
  //               tree.  This constructor will call SetBranchAddresses().

  // Set the TTree fDtupleTree using the input chain
  fDtupleTree = chain;

  // Set the branch addresses
  SetBranchAddresses();
}




//
// Destructon!
//
TDtuple::~TDtuple ()
{
  // Destructor: If you have not called Write() then you may have
  // NOT written the last file.  You MUST call Write() yourself.
  // This will delete the TTree* fDtupleTree

  // delete the TDtuple TTree
  //delete fDtupleTree;
}




TTree* TDtuple::GetTree ()
{
  // Get a pointer to the dtuple tree object

  return fDtupleTree;
}






void TDtuple::SetTree (TTree* inTree)
{
  // Get a pointer to the dtuple tree object

  fDtupleTree = inTree;
  return;
}






//
// Function: SetRootFile
//
// Purpose; To setup the output rootfile
//
// Arguments: string - the name of the output rootfile
//
// Return: void
// 
void TDtuple::SetRootFile (std::string const NtupleName)
{
  // Will setup the output root file given a string for the name of the file

  // Get a TFile object pointer with the given name
  fRootFile = new TFile(NtupleName.c_str(), "recreate");
  fRootFile->cd();

  // Setup branches and branch addresses
  SetBranches ();

  return;
}



//
// Function: SetRootFile
//
// Purpose; To setup the output rootfile
//
// Arguments: TFile* - pointer to the root output file
//
// Return: void
// 
void TDtuple::SetRootFile (TFile* fi)
{
  // Will setup the output root file given a TFile*.  This binds the Dtuple
  // TTree to the TFile

  // Set the root file pointer.
  fRootFile = fi;
  fRootFile->cd();

  // Setup branches and branch addresses
  SetBranches ();

  return;
}



//
// Function: SetMaxTreeSize
//
// Purpose: To set the max file size for our putput root file
//
// Arguments: int - the max size
//
// Return: void
//
void TDtuple::SetMaxTreeSize (int const size)
{
  // Set the max output file size for the output rootfile

  // Set the max treesize
  TTree::SetMaxTreeSize(size);

  return;
}


//
// Function: SetBranches
//
// Purpose: to set the branch addresses in the tree
//
// Arguments: none
//
// Return: void
//
void TDtuple::SetBranches ()
{
  // Creates TClonesArrays for our objects as well as sets the branches
  // for all of the event variables we wish to store to the dtuple.  This
  // binds the TDtuple tree "dtuple" to the output TFile

  // Make our dtuple tree
  fDtupleTree = new TTree("dtuple", "Common Analysis TDtuple");

  // Set the location of the dtuple tree
  fDtupleTree->SetDirectory(fRootFile);

  std::cout << "TDtuple: Setting TTree and Branches" << std::endl;
  std::cout << "TDtuple: fRootFile " << fRootFile->GetName() << std::endl;

  // Get Tclones array for each composite object
  fLepton = new TClonesArray("TLepton");
  fJet = new TClonesArray("TJet");
  fPhoton = new TClonesArray("TPhoton");

  // Set branches for event variables
  fDtupleTree->Branch("Run", &fRun, "Run/I");
  fDtupleTree->Branch("Event", &fEvent, "Event/I");
  fDtupleTree->Branch("RunSection", &fRunSection, "RunSection/I");
  fDtupleTree->Branch("EventFlags", &fEventFlags, "EventFlags/I");
  fDtupleTree->Branch("TriggerBits", &fTriggerBits, "TriggerBits/I");
  fDtupleTree->Branch("MetX", &fMetX, "MetX/F");
  fDtupleTree->Branch("MetY", &fMetY, "MetY/F");
  fDtupleTree->Branch("RawMetX", &fRawMetX, "RawMetX/F");
  fDtupleTree->Branch("RawMetY", &fRawMetY, "RawMetY/F");
  fDtupleTree->Branch("RawSumEt", &fRawSumEt, "RawSumEt/F");
  fDtupleTree->Branch("SumEt", &fSumEt, "SumEt/F");
  fDtupleTree->Branch("Lum", &fLum, "Lum/F");

  // Set branches for composite objects
  fDtupleTree->Branch("lepton", &fLepton);
  fDtupleTree->Branch("jet", &fJet);
  fDtupleTree->Branch("photon", &fPhoton);

  return;
}





//
// Function: SetBranchAddresses
//
// Purpose: to set the branch addresses in the tree
//
// Arguments: none
//
// Return: void
//
void TDtuple::SetBranchAddresses ()
{
  // Creates the TClonesArrays we need and sets the branch addresses for all of the variables.
  // This function should be used when you want to look at an existing dtuple TTree.

  // Get Tclones array for each composite object
  fLepton = new TClonesArray("TLepton");
  fJet = new TClonesArray("TJet");
  fPhoton = new TClonesArray("TPhoton");

  // Set branches for event variables
  fDtupleTree->SetBranchAddress("Run", &fRun);
  fDtupleTree->SetBranchAddress("Event", &fEvent);
  fDtupleTree->SetBranchAddress("RunSection", &fRunSection);
  fDtupleTree->SetBranchAddress("EventFlags", &fEventFlags);
  fDtupleTree->SetBranchAddress("TriggerBits", &fTriggerBits);
  fDtupleTree->SetBranchAddress("MetX", &fMetX);
  fDtupleTree->SetBranchAddress("MetY", &fMetY);
  fDtupleTree->SetBranchAddress("RawMetX", &fRawMetX);
  fDtupleTree->SetBranchAddress("RawMetY", &fRawMetY);
  fDtupleTree->SetBranchAddress("RawSumEt", &fRawSumEt);
  fDtupleTree->SetBranchAddress("SumEt", &fSumEt);
  fDtupleTree->SetBranchAddress("Lum", &fLum);

  // Set branches for composite objects
  fDtupleTree->SetBranchAddress("lepton", &fLepton);
  fDtupleTree->SetBranchAddress("jet", &fJet);
  fDtupleTree->SetBranchAddress("photon", &fPhoton);

  return;
}




//
// Function: GetEntry
//
// Purpose: to get an entry from the dtuple
//
// Arguments: int - the entry to get
//
// Return: int - returns TTree::GetEntry value
// 
int TDtuple::GetEntry(int const ientry)
{
  // Get an entry from the dtuple.  This calls the method TTree::GetEntry() on the dtuple tree.
  // This function also fills the Leptons and Jets vectors and sorts(not yet) them by Pt in decending order:
  // The Highest Pt object will be in Leptons[0]/Jets[0]

  // Value we will return
  int returnValue;

  // Clear to reset all dtuple variables
  Clear();

  // Actually get the entry from the TTree and save the return value
  returnValue = (int) fDtupleTree->GetEntry(ientry);

  // If the return value is zero something is wrong.  Don't try to fill anything just return zero
  if (returnValue == 0) {
    return returnValue;
  }

  // Fill the Leptons vector
  for (int ilepton=0; ilepton != fLepton->GetEntries(); ++ilepton) {
    TLepton *lep = (TLepton*) fLepton->At(ilepton);
    Leptons.push_back((TLepton) *lep);
  }


  // Sort Leptons by Pt - decending order
  SortLeptons(Leptons.begin(), Leptons.end());

  // Fill the Jets vector
  for (int ijet=0; ijet != fJet->GetEntries(); ++ijet) {
    TJet *jet = (TJet*) fJet->At(ijet);
    Jets.push_back(*jet);
  }

  // Sort Jets by Pt - decending order
  SortJets(Jets.begin(), Jets.end());

  // Fill the Photons vector
  for (int iphoton=0; iphoton != fPhoton->GetEntries(); ++iphoton) {
    TPhoton* photon = (TPhoton*) fPhoton->At(iphoton);
    Photons.push_back(*photon);
  }

  // Sort Photons by Pt - decending order
  SortPhotons(Photons.begin(), Photons.end());

  // Return the value from TTree::GetEntry()
  return returnValue;
}




//
// Function: CopyEventVarsFrom
//
// Purpose: to copy the event variables from the input dtuple
//          to this one
//
// Arguments: TDtuple&
//
// Return: void
// 
void TDtuple::CopyEventVarsTo (TDtuple* To)
{
  // This function is used to copy values for the event variables from
  // the input dtuple to this one.  You need to copy the others yourself
  // ie leptons, jets, photons, and so on must be done elsewhere

  To->SetRun(fRun);
  To->SetEvent(fEvent);
  To->SetRunSection(fRunSection);
  To->SetEventFlags(fEventFlags);
  To->SetTriggerBits(fTriggerBits);
  To->SetMetX(fMetX);
  To->SetMetY(fMetY);
  To->SetRawMetX(fRawMetX);
  To->SetRawMetY(fRawMetY);
  To->SetSumEt(fSumEt);
  To->SetRawSumEt(fRawSumEt);
  To->SetLum(fLum);

  return;
}




//
// Function: Fill
//
// Purpose: to save an event to the dtuple TTree
//
// Arguments: none
//
// Return: void
// 
void TDtuple::Fill()
{
  // Save an event to the Dtuple TTree.  This must be called once per event
  // you wish to save

  AddLeptonsToArray();
  // Fill the Dtuple
  fDtupleTree->Fill();

  return;
}





//
// Function: Write
//
// Purpose: To write the output rootfile
//
// Arguments: none
//
// Return: void
// 
void TDtuple::Write()
{
  // Write the dtuple rootfile.  This should be called only once.  If you
  // forget to call it you may in fact lose data.

  // Write the TTree
  //fDtupleTree->Write();
  fDtupleTree->AutoSave();

  return;
}





//
// Function: Clear
//
// Purpose: To clear the TClonesArray objects, reset object counters
//          and to set dummy values for the variables
//
// Arguments: none
//
// Return: void
//
void TDtuple::Clear ()
{
  // Clear the TClonesArray objects, reset object counters
  // and to set dummy values for the variables.  You should
  // call this once per event before you fill any of the variables

  // Clear the lepton and jet TClonesArrays
  fLepton->Clear();
  fJet->Clear();
  fPhoton->Clear();

  // Clear the Lepton and Jet vectors
  Leptons.clear();
  Jets.clear();
  Photons.clear();

  // Reset object counters
  //fNLeptonObjs = 0;
  //fNJetObjs = 0;


  return;
}





//
// Function: DefaultValues
//
// Purpose: To set dummy values to our dtuple variables
//
// Arguments: none
//
// Return: void
//
void TDtuple::DefaultValues()
{
  // Set dummy values to our dtuple variables.  All dtuple event
  // variables will be set to -9999 except for TriggerBits and EventFlags
  // which is set to 0x0

  fRun             = -9999;
  fEvent           = -9999;
  fRunSection      = -9999;
  fMetX            = -9999;
  fMetY            = -9999;
  fRawMetX         = -9999;
  fRawMetY         = -9999;
  fLum             = -9999;
  fTriggerBits     =  0x0;
  fEventFlags      =  0x0;

  Clear();

  return;
}





//
// Function: AddLeptons
//
// Purpose: To add a vector of leptons using iterators
//
// Arguments: iterators - for beginning and end of vector to add
//
// Return: void
//
void TDtuple::AddLeptons (std::vector<TLepton>::iterator it, std::vector<TLepton>::iterator end)
{
  // Add a whole vector of TLepton objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Run through each lepton and add it
  for ( ; it != end; ++it) {
    AddLepton(*it);
  }

  return;
}




//
// Function: AddLeptons
//
// Purpose: To add a vector of leptons all at once to the TClonesArrays
//
// Arguments: vector<Lepton> - a vector of Lepton objects
//
// Return: void
//
void TDtuple::AddLeptons (std::vector<TLepton> const& inleps)
{
  // Add a whole vector of TLepton objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Loop over the vector and add each TLepton
  for (size_t i=0; i != inleps.size(); ++i) {
    AddLepton(inleps[i]);
  }

  return;
}



//
// Function: AddLepton
//
// Purpose: To add one lepton object to the dtuple TLepton vector
//
// Arguments: TLepton - The TLepton you wish to add
//
// Return: void
//
void TDtuple::AddLepton (TLepton const& inlep)
{
  // Add a TLepton the the dtuple lepton 

  // Add the lepton to the vector as well.  This is actually useful.
  Leptons.push_back(inlep);

  // Create a new entry in the TCLonesArray
  //TLepton* newLepton = new ((*fLepton)[fLepton->GetEntries()]) TLepton();

  // Set our new entry with all of the input lepton variables
  //*newLepton = inlep;

  return;
}




void TDtuple::AddLeptonsToArray ()
{
  for (size_t i = 0; i != Leptons.size(); ++i) {
  // Create a new entry in the TCLonesArray
  TLepton* newLepton = new ((*fLepton)[fLepton->GetEntries()]) TLepton();

  // Set our new entry with all of the input lepton variables
  *newLepton = Leptons[i];
  }

  return;
}
//
// Function: AddJets
//
// Purpose: To add a vector of jets all at once to the TClonesArrays via iterators
//
// Arguments: iterators - for beginning and end of vector to add
//
// Return: void
//
void TDtuple::AddJets (std::vector<TJet>::iterator it, std::vector<TJet>::iterator end)
{
  // Add a whole vector of TJet objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Run through each jet and add it
  for ( ; it != end; ++it) {
    AddJet(*it);
  }

  return;
}




//
// Function: AddJets
//
// Purpose: To add a vector of TJet objects to the dtuple TJet TClonesArray
//
// Arguments: vector<TJet> - The vector of TJets you wish to add
//
// Return: void
//
void TDtuple::AddJets (std::vector<TJet> const& injets)
{
  // Add a whole vector of TJet objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Loop over the vector and add each jet
  for (size_t i=0; i != injets.size(); ++i) {
    AddJet(injets[i]);
  }

  return;
}



//
// Function: AddJet
//
// Purpose: To add one TJet object to the dtuple TJet TClonesArray
//
// Arguments: TJet - The TJet you wish to add
//
// Return: void
//
void TDtuple::AddJet (TJet const& injet)
{
  // Add a TJet to the dtuple jet TClonesArray
  // Note: This does not write them, you still need to call Fill()

  // Add the jet to the vector as well.  This is actually useful.
  Jets.push_back(injet);

  // Create a new jet object in the TClonesArray
  TJet* newJet = new ((*fJet)[fJet->GetEntries()]) TJet();

  // Set the new object variables with the input Jet
  *newJet = injet;

  return;
}




//
// Function: AddPhotons
//
// Purpose: To add a vector of jets all at once to the TClonesArrays via iterators
//
// Arguments: iterators - for beginning and end of vector to add
//
// Return: void
//
void TDtuple::AddPhotons (std::vector<TPhoton>::iterator it, std::vector<TPhoton>::iterator end)
{
  // Add a whole vector of Photon objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Run through each jet and add it
  for ( ; it != end; ++it) {
    AddPhoton(*it);
  }

  return;
}




//
// Function: AddPhotons
//
// Purpose: To add a vector of Photon objects to the dtuple Photon TClonesArray
//
// Arguments: vector<Photon> - The vector of Photons you wish to add
//
// Return: void
//
void TDtuple::AddPhotons (std::vector<TPhoton> const& inphotons)
{
  // Add a whole vector of Photon objects to the dtuple TClonesArray.
  // Note: This does not write them, you still need to call Write()

  // Loop over the vector and add each jet
  for (size_t i=0; i != inphotons.size(); ++i) {
    AddPhoton(inphotons[i]);
  }

  return;
}



//
// Function: AddPhoton
//
// Purpose: To add one Photon object to the dtuple Photon TClonesArray
//
// Arguments: Photon - The Photon you wish to add
//
// Return: void
//
void TDtuple::AddPhoton (TPhoton const& inphoton)
{
  // Add a Photon to the dtuple jet TClonesArray
  // Note: This does not write them, you still need to call Fill()

  // Add the jet to the vector as well.  This is actually useful.
  Photons.push_back(inphoton);

  // Create a new jet object in the TClonesArray
  TPhoton* newPhoton = new ((*fPhoton)[fPhoton->GetEntries()]) TPhoton();

  // Set the new object variables with the input Photon
  *newPhoton = inphoton;

  return;
}




//
// Function: SortLeptons
//
// Purpose: Sort a TLepton vector by Pt
//
// Arguments: iter - beginning of vector to sort
//            iter - end of vetor to sort
//
// Return: void
//
void  TDtuple::SortLeptons (std::vector<TLepton>::iterator start, std::vector<TLepton>::iterator end)
{
  // Sort a vector of TLeptons by Et with the highest Et lepton as the first element

  for (std::vector<TLepton>::iterator thispos = start; thispos != end; ++thispos) {
    for (std::vector<TLepton>::iterator opos = thispos; opos != end; ++opos) {
      if (opos->Perp() > thispos->Perp()) {
        std::swap_ranges(thispos, thispos+1, opos);
      }
    }
  }

  return;
}




//
// Function: SortJets
//
// Purpose: Sort a TJet vector by Pt
//
// Arguments: iter - beginning of vector to sort
//            iter - end of vetor to sort
//
// Return: void
//
void TDtuple::SortJets (std::vector<TJet>::iterator start, std::vector<TJet>::iterator end)
{
  // Sort a vector of TJets by Pt with the highest Pt jet as the first element

  for (std::vector<TJet>::iterator thispos = start; thispos != end; ++thispos) {
    for (std::vector<TJet>::iterator opos = thispos; opos != end; ++opos) {
      if (opos->Perp() > thispos->Perp()) {
        std::swap_ranges(thispos, thispos+1, opos);
      }
    }
  }

  return;
}





//
// Function: SortPhotons
//
// Purpose: Sort a Photon vector by Pt
//
// Arguments: iter - beginning of vector to sort
//            iter - end of vetor to sort
//
// Return: void
//
void TDtuple::SortPhotons (std::vector<TPhoton>::iterator start, std::vector<TPhoton>::iterator end)
{
  // Sort a vector of Photons by Pt with the highest Pt jet as the first element

  for (std::vector<TPhoton>::iterator thispos = start; thispos != end; ++thispos) {
    for (std::vector<TPhoton>::iterator opos = thispos; opos != end; ++opos) {
      if (opos->Perp() > thispos->Perp()) {
        std::swap_ranges(thispos, thispos+1, opos);
      }
    }
  }

  return;
}





//
// Function: SetRun
//
// Purpose: To set the internal variable for run number
//
// Arguments: int - run number
//
// Return: void
//
void  TDtuple::SetRun (int const in)
{
  // To set the dtuple variable for run number

  fRun = in;
  return;
}






//
// Function: SetEvent
//
// Purpose: To set the internal variable for event number
//
// Arguments: int - event number
//
// Return: void
//
void  TDtuple::SetEvent (int const in)
{
  // To set the internal variable for event number

  fEvent = in;
  return;
}






//
// Function: SetRunSection
//
// Purpose: To set the internal variable for run section
//
// Arguments: int - run section
//
// Return: void
//
void  TDtuple::SetRunSection (int const in)
{
  // To set the dtuple variable for run section

  fRunSection = in;
  return;
}






//
// Function: SetMetXY
//
// Purpose: To set the internal variables for x and y met components
//
// Arguments: float - x-compoent of missint Et
//            float - y-compoent of missint Et
//
// Return: void
//
void  TDtuple::SetMetXY (float const x, float const y)
{
  // To set the internal variables for x and y components of the missing Et

  SetMetX(x);
  SetMetY(y);

  return;
}






//
// Function: SetRawMetXY
//
// Purpose: To set the internal variables for x and y RawMet components
//
// Arguments: float - x-compoent of raw missint Et
//            float - y-compoent of raw missint Et
//
// Return: void
//
void  TDtuple::SetRawMetXY (float const x, float const y)
{
  // To set the internal variables for x and y components of the raw missing Et

  SetRawMetX(x);
  SetRawMetY(y);

  return;
}






//
// Function: SetMetX
//
// Purpose: To set the internal variable for x component of the missing Et
//
// Arguments: float - x-compoent of missint Et
//
// Return: void
//
void  TDtuple::SetMetX (float const in)
{
  // To set the internal variable for x component of the missing Et

  fMetX = in;
  return;
}






//
// Function: SetMetY
//
// Purpose: To set the internal variable for y component of the missing Et
//
// Arguments: float - y-compoent of missint Et
//
// Return: void
//
void  TDtuple::SetMetY (float const in)
{
  // To set the interal variable for y coponent of missing Et

  fMetY = in;
  return;
}





//
// Function: SetRawMetX
//
// Purpose: To se the internal dtuple variable for Raw Missing Et
//
// Arguments: float - Raw Missing Et
//
// Return: void
//
void  TDtuple::SetRawMetX (float const in)
{
  // To se the internal dtuple variable for x-component of Raw Missing Et

  fRawMetX = in;
  return;
}




//
// Function: SetRawMetY
//
// Purpose: To se the internal dtuple variable for Raw Missing Et
//
// Arguments: float - Raw Missing Et
//
// Return: void
//
void  TDtuple::SetRawMetY (float const in)
{
  // To se the internal dtuple variable for y-component of Raw Missing Et

  fRawMetY = in;
  return;
}




//
// Function: SetSumEt
//
// Purpose: Set the internal dtuple variable for SumEt
//
// Arguments: float - Sum Et
//
// Return: void
void TDtuple::SetSumEt (float const in)
{
  // Set the internal dtuple variable for SumEt

  fSumEt = in;
  return;
}





//
// Function: SetRawSumEt
//
// Purpose: Set the internal dtuple variable for RawSumEt
//
// Arguments: float - Sum Et
//
// Return: void
void TDtuple::SetRawSumEt (float const in)
{
  // Set the internal dtuple variable for RawSumEt

  fRawSumEt = in;
  return;
}





//
// Function: SetLum
//
// Purpose: Set the internal variable for instantaneous luminosity
//
// Arguments: float - luminosity
//
// Return: void
//
void TDtuple::SetLum (float const in)
{
  // Set the internal dtuple variable for instantaneous luminosity

  fLum = in;
  return;
}





//
// Function: GetRun
//
// Purpose: Get the run number
//
// Arguments: none
//
// Return: int - Run number
//
int TDtuple::GetRun ()
{
  // Get the run number

  return fRun;
}






//
// Function: GetEvent
//
// Purpose: Get the run number
//
// Arguments: none
//
// Return: int - Event number
//
int TDtuple::GetEvent ()
{
  // Get the event number

  return fEvent;
}





//
// Function: GetRunSection
//
// Purpose: Get the run section
//
// Arguments: none
//
// Return: int - Run section
//
int TDtuple::GetRunSection ()
{
  // Get the run section

  return fRunSection;
}






//
// Function: GetMetX
//
// Purpose: Get the x-component of the Missing Et
//
// Arguments: none
//
// Return: float - x-component of Missing Et
//
float TDtuple::GetMetX ()
{
  // Get the x-component of the Missing Et

  return fMetX;
}






//
// Function: GetMetY
//
// Purpose: Get the y-component of the Missing Et
//
// Arguments: none
//
// Return: float - y-component of Missing Et
//
float TDtuple::GetMetY ()
{
  // Get the y-component of the Missing Et

  return fMetY;
}




//
// Function: GetMet
//
// Purpose: Get the magnitude of the Missing Et
//
// Arguments: none
//
// Return: float - Missing Et magnitude
//
float TDtuple::GetMet ()
{
  // Get the magnitude of the missing Et.  MetX and MetY
  // must already be set

  return sqrt( pow(fMetX, 2) + pow(fMetY, 2) );
}




//
// Function: GetMetPhi
//
// Purpose: Get the phi of the Met
//
// Arguments: none
//
// Return: float - Missing Et phi direction
//
float TDtuple::GetMetPhi ()
{
  // Get the magnitude of the missing Et.  MetX and MetY
  // must already be set

  return TVector2(fMetX, fMetY).Phi();
}





//
// Function: GetRawMet
//
// Purpose: Get the internal dtuple variable for RawMet
//
// Arguments: none
//
// Return: float - RawMet
//
float TDtuple::GetRawMet ()
{
  // Get the magnitude of the Raw Met

  return sqrt( pow(fRawMetX, 2) + pow(fRawMetY, 2) );
}





//
// Function: GetRawMetX
//
// Purpose: Get the internal dtuple variable for RawMetX
//
// Arguments: none
//
// Return: float - RawMetX
//
float TDtuple::GetRawMetX ()
{
  // Get the internal dtuple variable for RawMetX

  return fRawMetX;
}





//
// Function: GetRawMetY
//
// Purpose: Get the internal dtuple variable for RawMetY
//
// Arguments: none
//
// Return: float - RawMetY
//
float TDtuple::GetRawMetY ()
{
  // Get the internal dtuple variable for RawMetY

  return fRawMetY;
}





//
// Function: GetSumEt
//
// Purpose: Get the internal dtuple variable for SumEt
//
// Arguments: none
//
// Return: float - the SumEt
//
float TDtuple::GetSumEt ()
{
  // Get the internal dtuple variable for SumEt

  return fSumEt;
}





//
// Function: GetRawSumEt
//
// Purpose: Get the internal dtuple variable for RawSumEt
//
// Arguments: none
//
// Return: float - the SumEt
//
float TDtuple::GetRawSumEt ()
{
  // Get the internal dtuple variable for RawSumEt

  return fRawSumEt;
}





//
// Function: GetLum
//
// Purpose: Get the interna dtuple variable for instantaneous luminosity
//
// Arguments: none
//
// Return: float - instantaneous luminosity
//
float TDtuple::GetLum ()
{
  // Get the interna dtuple variable for instantaneous luminosity

  return fLum;
}





//
// Function: GetTriggerBit
//
// Purpose: To see if a bit is set
//
// Arguments: int - which bitto ckeck.
//
// Return: bool - true if bit is set, false otherwise
//
bool TDtuple::GetTriggerBit (std::string const trigger)
{
  // Get a specific bit in the trigger bits dtuple variable

  // Get the trigger map
  //std::map<std::string, int> TriggerMap = GetTriggerBitMap();

  // And it with fTriggerBits
  //return fTriggerBits & TriggerMap[trigger];
  return false;
}





//
// Function: SetTriggerBit
//
// Purpose: To set a bit in the trigger bits int
//
// Arguments: string - which trigger to set
//
// Return: void
//
void TDtuple::SetTriggerBit (std::string const trigger)
{
  // Set a specific bit in the trigger bits dtuple variable

  // Get the trigger map
  //std::map<std::string, int> TriggerMap = GetTriggerBitMap();

  // orequal to set the bit
  //fTriggerBits |= TriggerMap[trigger];
  return;
}





void TDtuple::SetTriggerBits (int const triggers)
{
  fTriggerBits = triggers;
  return;
}





//
// Function: GetTriggerBitMap
//
// Purpose: To get the list of triggers and bits we are saving them to
//
// Arguments: none
//
// Return map - strings and ints for the trigger names and which bits they correspond to
//

//
// Function: SetEventFlag
//
// Purpose: To set a specific event flag
//
// Arguments: string - the flag to set
//
// Return: void
//
void TDtuple::SetEventFlag (std::string const flag)
{
  // Set an event flag in the Dtuple.  Flags are defined in GetEventFlagMap()

  //static std::map<std::string, int> FlagMap = GetEventFlagMap();
  //fEventFlags |= FlagMap[flag];

  return;
}






void TDtuple::SetEventFlags (int const flags)
{
  fEventFlags = flags;
  return;
}






//
// Function: GetEventFlag
//
// Purpose: To get a specific event flag
//
// Arguments: string - the flag to get
//
// Return: bool - true if flag is set
//
bool TDtuple::GetEventFlag (std::string const flag)
{
  // Test if a flag has been set or not.  Flags are defined in GetEventFlagMap()

  // Get the Event Flag Map
  //static std::map<std::string, int> FlagMap = GetEventFlagMap();

  // True if the flag is set in fEventFlags, false otherwise
  //return (fEventFlags & FlagMap[flag]) == 0 ? false : true;
  return false;
}





//
// Function: GetEventFlag
//
// Purpose: To get the event flag int
//
// Arguments:
//
// Return: int - the event flag int
//
int TDtuple::GetEventFlags ()
{
  // Return the value of the internal variable used to store event flags

  return (fEventFlags);
}





//
// Function: GetTriggerBits
//
// Purpose: To get the trigger bit int
//
// Arguments:
//
// Return: int - the trigger bit int
//
int TDtuple::GetTriggerBits ()
{
  // Return the value of the internal variable used to store trigger bits

  return (fTriggerBits);
}





std::vector<TLepton>* TDtuple::GetLeptons ()
{
  // Get the vector of leptons

  return &Leptons;
}




std::vector<TJet>* TDtuple::GetJets ()
{
  // Get the vector of jets

  return &Jets;
}




std::vector<TPhoton>* TDtuple::GetPhotons ()
{
  // Get the vector of photons

  return &Photons;
}








