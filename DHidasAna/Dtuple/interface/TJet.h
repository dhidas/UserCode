#ifndef GUARD_TJet_h
#define GUARD_TJet_h


#include <iostream>
#include <string>
#include <ostream>

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"



//
// Class: TJet
//
// Purpose: This is a class to hold information about TJets.  It
//          does so by storing some information to private variables
//          and then has methods to access those variables as well as
//          methods which combine these basic variables into other quantaties
//          which one might be interested in.  The basic variables which should
//          be filled are listed below.  They are private so as to be nondestructive.
//          In other words, it's harder to mess them up...
//
//          Of you are compiling this internal to root you may want to uncomment the line
//          in this .h file with the ClassDef (no, do not put a semicolon(;) after it,
//          it's a root macro, not a function.
//
//          To otherwise compile this class with gcc I would reccomend the following
//          (on cdf machines):
//
//          source ~cdfsoft/cdf2.cshrc
//          setup gcc v3_4_3
//          setup root v4_02_00a -q GCC_3_4_3
//          g++ `root-config --cflags` -c TJet.cc -o TJet.o
//          
class TJet: public TLorentzVector
{
  // Constructors and destructor.
  // More of these will be added in time.
  public:
    TJet ();
    ~TJet ();


  // The private variables which need to be
  // set by the various methods
  private:
    float DetEta;       // Detector Eta
    float EScaleFactor; // Jet Energy Scale Factor
    float EmF;          // Em Fraction
    float HadF;         // Em Fraction
    float SysPos;       // Positive Systematic uncertainty
    float SysNeg;       // Negative Systematic uncertainty
    float Z0;           // Z0 of the jet
    float dxy;          // Signed displacement of vertex
    float Charge;       // The "charge" of the jet, whatever that is



  // Public methods for accessing the variables in
  // this class
  public:

    // Get quantaties
    float GetDetEta ();
    float GetEmF ();
    float GetHadF ();
    float GetEScaleFactor ();
    float GetSysPos ();
    float GetSysNeg ();
    float GetZ0 ();
    float Getdxy ();
    int   GetCharge();


    // Set quantaties
    void  SetDetEta (float const);
    void  SetEScaleFactor (float const);
    void  SetEmF (float const);
    void  SetHadF (float const);
    void  SetSysPos (float const);
    void  SetSysNeg (float const);
    void  SetZ0 (float const);
    void  Setdxy (float const);
    void  SetCharge (int const);

    void  ReScale (float const);

    // operator def for TJet
    int operator<(TJet&);

  private:
    void  DefaultValues ();


  public:
    // Root Likes ClassDef and ClassImp.
    // Comment them out if you don't need them.
    // There should NOT be a ; since this is a root macro and not a function
    // ClassDef must be the last line of the class before the };
    ClassDef(TJet,3) // TJet class
};


#endif
