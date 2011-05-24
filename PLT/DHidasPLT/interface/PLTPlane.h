#ifndef GUARD_PLTPlane_h
#define GUARD_PLTPlane_h

#include "PLTCluster.h"

#include "TH2F.h"


class PLTPlane
{
  public:
    PLTPlane ();
    ~PLTPlane ();

    void   AddHit (PLTHit*);
    float  Charge ();
    int    Channel ();
    bool   AddClusterFromSeed (PLTHit*);
    bool   IsBiggestHitIn3x3 (PLTHit*);
    void   Clusterize ();
    TH2F*  DrawHits2D ();
    size_t NHits ();
    int    ROC ();

  private:
    int fChannel;
    int fROC;
    std::vector<PLTHit*> fHits;
    std::vector<PLTCluster*> fClusters;

};










#endif
