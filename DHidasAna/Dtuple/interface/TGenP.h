#ifndef GUARD_TGenP_h
#define GUARD_TGenP_h

#include "TLorentzVector.h"

class TGenP : public TLorentzVector
{
  public:
    TGenP ();
    ~TGenP ();

  private:
    int Id;
    int MotherId;

  public:
    int GetId ();
    int GetMotherId ();

    void SetId (int const);
    void SetMotherId (int const);

    bool IsId (int const);
    bool IsMotherId (int const);


  public:
    // Root Likes ClassDef and ClassImp.
    // Comment them out if you don't need them.
    // There should NOT be a ; since this is a root macro and not a function
    // ClassDef must be the last line of the class before the };
    ClassDef(TGenP,1) // TGenP class
};














#endif
