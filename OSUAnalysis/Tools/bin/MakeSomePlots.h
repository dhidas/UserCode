#include "TTree.h"

    double weight;
double PileUpWeight;
    // general event information: 
    int    type,runNumber, numberOfJets, numberOfBJets, eSEL, muSEL, leptoCharge;
unsigned long eventNumber;
double ST, MET, HT;
int nPileUpVtx;
    // leading non-b-jet:
    int    leadingJetPdgId, leadingJetIndGen;
    double leadingJetPtGen, leadingJetEtaGen, leadingJetPhiGen;
    double leadingJetPtRec, leadingJetEtaRec, leadingJetPhiRec;

    double JetPx[100];
    double JetPy[100];
    double JetPz[100];
double JetE[100];
int    JetBTag[100];
double BJetPx[100];
double BJetPy[100];
double BJetPz[100];
double BJetE[100];

    // highest pT leftover jet (gen-level: radiation for the background and d-jet in signals)
    int    freeJetPdgId,    freeJetIndGen;
    double freeJetPtGen,    freeJetEtaGen,    freeJetPhiGen;
    double freeJetPtRec,    freeJetEtaRec,    freeJetPhiRec;
    // a candidate b-jet from the leptonic top:
    double bJetTlPtGen, bJetTlEtaGen, bJetTlPhiGen;
    double bJetTlPtRec, bJetTlEtaRec, bJetTlPhiRec;
    // a candidate b-jet from the hadronic top:
    double bJetThPtGen, bJetThEtaGen, bJetThPhiGen;
    double bJetThPtRec, bJetThEtaRec, bJetThPhiRec;
    // first and second (ordered in pT) jets from the hadronic W
    int    jet1WhPdgId, jet2WhPdgId,  jet1WhIndGen, jet2WhIndGen;
    double jet1WhPtGen, jet1WhEtaGen, jet1WhPhiGen;
    double jet1WhPtRec, jet1WhEtaRec, jet1WhPhiRec;
    double jet2WhPtGen, jet2WhEtaGen, jet2WhPhiGen;
    double jet2WhPtRec, jet2WhEtaRec, jet2WhPhiRec;
    // neutrino and missing ET
    double metPtGen,   metPhiGen,   metEtaGen;
    double metPtRec,   metPhiRec,   metEtaRec;
    // lepton
    int    leptonPdgId;
    double leptonPtGen, leptonEtaGen, leptonPhiGen;
    double leptonPtRec, leptonEtaRec, leptonPhiRec;
    // composite objects 
    //  leptonic and hadronic W and top
    int    tlPdgId,    thPdgId,     wlPdgId,   whPdgId;
    double wlPtGen,    wlEtaGen,    wlPhiGen,  wlMgen;
    double wlPtRec,    wlEtaRec,    wlPhiRec,  wlMrec;
    double tlPtGen,    tlEtaGen,    tlPhiGen,  tlMgen;
    double tlPtRec,    tlEtaRec,    tlPhiRec,  tlMrec;
    double whPtGen,    whEtaGen,    whPhiGen,  whMgen;
    double whPtRec,    whEtaRec,    whPhiRec,  whMrec;
    double thPtGen,    thEtaGen,    thPhiGen,  thMgen;
    double thPtRec,    thEtaRec,    thPhiRec,  thMrec;
    //  W' candidate
    int    wpPdgId;
    double wpPtGen,    wpEtaGen,    wpPhiGen;
    double tPosJetMassGen, tPosJetPtGen, tPosJetEtaGen, tPosJetPhiGen;
    double tNegJetMassGen, tNegJetPtGen, tNegJetEtaGen, tNegJetPhiGen;
    double tPosJetMassRec, tPosJetPtRec, tPosJetEtaRec, tPosJetPhiRec;
    double tNegJetMassRec, tNegJetPtRec, tNegJetEtaRec, tNegJetPhiRec;

    // rec-gen matching for the jets used in the top reco 
    double dRtlBjet, dRthBjet, dRwhJet1, dRwhJet2;
    // matching of gen jets to any jets around
    //double matchRtlJet, matchRthJet, matchRwhJet1, matchRwhJet2;


void initMicroNtuple(TTree* microTuple) {
    microTuple->SetBranchAddress("type",   &type);
    microTuple->SetBranchAddress("weight", &weight);
    microTuple->SetBranchAddress("PileUpWeight", &PileUpWeight);
    microTuple->SetBranchAddress("eventNumber",   &eventNumber);
    microTuple->SetBranchAddress("runNumber",   &runNumber);
    microTuple->SetBranchAddress("numberOfJets",  &numberOfJets);
    microTuple->SetBranchAddress("numberOfBJets", &numberOfBJets);
    microTuple->SetBranchAddress("leptoCharge",   &leptoCharge);
    microTuple->SetBranchAddress("eSEL",          &eSEL);
    microTuple->SetBranchAddress("muSEL",         &muSEL);
    microTuple->SetBranchAddress("ST",            &ST);
    microTuple->SetBranchAddress("HT",            &HT);
    microTuple->SetBranchAddress("MET",           &MET);
    microTuple->SetBranchAddress("nPileUpVtx",            &nPileUpVtx);
    microTuple->SetBranchAddress("leadingJetPdgId", &leadingJetPdgId);
    microTuple->SetBranchAddress("leadingJetIndGen",&leadingJetIndGen);
    microTuple->SetBranchAddress("leadingJetPtGen", &leadingJetPtGen);
    microTuple->SetBranchAddress("leadingJetEtaGen",&leadingJetEtaGen);
    microTuple->SetBranchAddress("leadingJetPhiGen",&leadingJetPhiGen);
    microTuple->SetBranchAddress("leadingJetPtRec", &leadingJetPtRec);
    microTuple->SetBranchAddress("leadingJetEtaRec",&leadingJetEtaRec);
    microTuple->SetBranchAddress("leadingJetPhiRec",&leadingJetPhiRec);

    microTuple->SetBranchAddress("JetPx[numberOfJets]",JetPx);
    microTuple->SetBranchAddress("JetPy[numberOfJets]",JetPy);
    microTuple->SetBranchAddress("JetPz[numberOfJets]",JetPz);
    microTuple->SetBranchAddress("JetE[numberOfJets]",JetE);
    microTuple->SetBranchAddress("JetBTag[numberOfJets]",JetBTag);

    microTuple->SetBranchAddress("BJetPx[numberOfBJets]",BJetPx);
    microTuple->SetBranchAddress("BJetPy[numberOfBJets]",BJetPy);
    microTuple->SetBranchAddress("BJetPz[numberOfBJets]",BJetPz);
    microTuple->SetBranchAddress("BJetE[numberOfBJets]",BJetE);


    microTuple->SetBranchAddress("freeJetPdgId",    &freeJetPdgId);
    microTuple->SetBranchAddress("freeJetIndGen",   &freeJetIndGen);
    microTuple->SetBranchAddress("freeJetPtGen",    &freeJetPtGen);
    microTuple->SetBranchAddress("freeJetEtaGen",   &freeJetEtaGen);
    microTuple->SetBranchAddress("freeJetPhiGen",   &freeJetPhiGen);
    microTuple->SetBranchAddress("freeJetPtRec",    &freeJetPtRec);
    microTuple->SetBranchAddress("freeJetEtaRec",   &freeJetEtaRec);
    microTuple->SetBranchAddress("freeJetPhiRec",   &freeJetPhiRec);

    microTuple->SetBranchAddress("bJetTlPtGen",    &bJetTlPtGen);
    microTuple->SetBranchAddress("bJetTlEtaGen",   &bJetTlEtaGen);
    microTuple->SetBranchAddress("bJetTlPhiGen",   &bJetTlPhiGen);
    microTuple->SetBranchAddress("bJetTlPtRec",    &bJetTlPtRec);
    microTuple->SetBranchAddress("bJetTlEtaRec",   &bJetTlEtaRec);
    microTuple->SetBranchAddress("bJetTlPhiRec",   &bJetTlPhiRec);

    microTuple->SetBranchAddress("bJetThPtGen",    &bJetThPtGen);
    microTuple->SetBranchAddress("bJetThEtaGen",   &bJetThEtaGen);
    microTuple->SetBranchAddress("bJetThPhiGen",   &bJetThPhiGen);
    microTuple->SetBranchAddress("bJetThPtRec",    &bJetThPtRec);
    microTuple->SetBranchAddress("bJetThEtaRec",   &bJetThEtaRec);
    microTuple->SetBranchAddress("bJetThPhiRec",   &bJetThPhiRec);

    microTuple->SetBranchAddress("jet1WhPdgId",    &jet1WhPdgId);
    microTuple->SetBranchAddress("jet2WhPdgId",    &jet2WhPdgId);
    microTuple->SetBranchAddress("jet1WhIndGen",   &jet1WhIndGen);
    microTuple->SetBranchAddress("jet2WhIndGen",   &jet2WhIndGen);
    microTuple->SetBranchAddress("jet1WhPtGen",    &jet1WhPtGen);
    microTuple->SetBranchAddress("jet1WhEtaGen",   &jet1WhEtaGen);
    microTuple->SetBranchAddress("jet1WhPhiGen",   &jet1WhPhiGen);
    microTuple->SetBranchAddress("jet1WhPtRec",    &jet1WhPtRec);
    microTuple->SetBranchAddress("jet1WhEtaRec",   &jet1WhEtaRec);
    microTuple->SetBranchAddress("jet1WhPhiRec",   &jet1WhPhiRec);
    microTuple->SetBranchAddress("jet2WhPtGen",    &jet2WhPtGen);
    microTuple->SetBranchAddress("jet2WhEtaGen",   &jet2WhEtaGen);
    microTuple->SetBranchAddress("jet2WhPhiGen",   &jet2WhPhiGen);
    microTuple->SetBranchAddress("jet2WhPtRec",    &jet2WhPtRec);
    microTuple->SetBranchAddress("jet2WhEtaRec",   &jet2WhEtaRec);
    microTuple->SetBranchAddress("jet2WhPhiRec",   &jet2WhPhiRec);

    microTuple->SetBranchAddress("metPtGen",       &metPtGen);
    microTuple->SetBranchAddress("metPhiGen",      &metPhiGen);
    microTuple->SetBranchAddress("metEtaGen",      &metEtaGen);
    microTuple->SetBranchAddress("metPtRec",       &metPtRec);
    microTuple->SetBranchAddress("metPhiRec",      &metPhiRec);
    microTuple->SetBranchAddress("metEtaRec",      &metEtaRec);

    microTuple->SetBranchAddress("leptonPdgId",    &leptonPdgId);
    microTuple->SetBranchAddress("leptonPtGen",    &leptonPtGen);
    microTuple->SetBranchAddress("leptonEtaGen",   &leptonEtaGen);
    microTuple->SetBranchAddress("leptonPhiGen",   &leptonPhiGen);
    microTuple->SetBranchAddress("leptonPtRec",    &leptonPtRec);
    microTuple->SetBranchAddress("leptonEtaRec",   &leptonEtaRec);
    microTuple->SetBranchAddress("leptonPhiRec",   &leptonPhiRec);

    microTuple->SetBranchAddress("wlPdgId",        &wlPdgId);
    microTuple->SetBranchAddress("wlPtGen",        &wlPtGen);
    microTuple->SetBranchAddress("wlEtaGen",       &wlEtaGen);
    microTuple->SetBranchAddress("wlPhiGen",       &wlPhiGen);
    microTuple->SetBranchAddress("wlMgen",         &wlMgen);
    microTuple->SetBranchAddress("wlPtRec",        &wlPtRec);
    microTuple->SetBranchAddress("wlEtaRec",       &wlEtaRec);
    microTuple->SetBranchAddress("wlPhiRec",       &wlPhiRec);
    microTuple->SetBranchAddress("wlMrec",         &wlMrec);

    microTuple->SetBranchAddress("tlPdgId",        &tlPdgId);
    microTuple->SetBranchAddress("tlPtGen",        &tlPtGen);
    microTuple->SetBranchAddress("tlEtaGen",       &tlEtaGen);
    microTuple->SetBranchAddress("tlPhiGen",       &tlPhiGen);
    microTuple->SetBranchAddress("tlMgen",         &tlMgen);
    microTuple->SetBranchAddress("tlPtRec",        &tlPtRec);
    microTuple->SetBranchAddress("tlEtaRec",       &tlEtaRec);
    microTuple->SetBranchAddress("tlPhiRec",       &tlPhiRec);
    microTuple->SetBranchAddress("tlMrec",         &tlMrec);

    microTuple->SetBranchAddress("whPdgId",        &whPdgId);
    microTuple->SetBranchAddress("whPtGen",        &whPtGen);
    microTuple->SetBranchAddress("whEtaGen",       &whEtaGen);
    microTuple->SetBranchAddress("whPhiGen",       &whPhiGen);
    microTuple->SetBranchAddress("whMgen",         &whMgen);
    microTuple->SetBranchAddress("whPtRec",        &whPtRec);
    microTuple->SetBranchAddress("whEtaRec",       &whEtaRec);
    microTuple->SetBranchAddress("whPhiRec",       &whPhiRec);
    microTuple->SetBranchAddress("whMrec",         &whMrec);

    microTuple->SetBranchAddress("thPdgId",        &thPdgId);
    microTuple->SetBranchAddress("thPtGen",        &thPtGen);
    microTuple->SetBranchAddress("thEtaGen",       &thEtaGen);
    microTuple->SetBranchAddress("thPhiGen",       &thPhiGen);
    microTuple->SetBranchAddress("thMgen",         &thMgen);
    microTuple->SetBranchAddress("thPtRec",        &thPtRec);
    microTuple->SetBranchAddress("thEtaRec",       &thEtaRec);
    microTuple->SetBranchAddress("thPhiRec",       &thPhiRec);
    microTuple->SetBranchAddress("thMrec",         &thMrec);

    microTuple->SetBranchAddress("wpPdgId",        &wpPdgId);
    microTuple->SetBranchAddress("wpPtGen",        &wpPtGen);
    microTuple->SetBranchAddress("wpEtaGen",       &wpEtaGen);
    microTuple->SetBranchAddress("wpPhiGen",       &wpPhiGen);

    microTuple->SetBranchAddress("dRtlBjet",&dRtlBjet);
    microTuple->SetBranchAddress("dRthBjet",&dRthBjet);
    microTuple->SetBranchAddress("dRwhJet1",&dRwhJet1);
    microTuple->SetBranchAddress("dRwhJet2",&dRwhJet2);

    microTuple->SetBranchAddress("tPosJetMassGen", &tPosJetMassGen);
    microTuple->SetBranchAddress("tPosJetPtGen",   &tPosJetPtGen);
    microTuple->SetBranchAddress("tPosJetEtaGen",  &tPosJetEtaGen);
    microTuple->SetBranchAddress("tPosJetPhiGen",  &tPosJetPhiGen);

    microTuple->SetBranchAddress("tNegJetMassGen", &tNegJetMassGen);
    microTuple->SetBranchAddress("tNegJetPtGen",   &tNegJetPtGen);
    microTuple->SetBranchAddress("tNegJetEtaGen",  &tNegJetEtaGen);
    microTuple->SetBranchAddress("tNegJetPhiGen",  &tNegJetPhiGen);

    microTuple->SetBranchAddress("tPosJetMassRec", &tPosJetMassRec);
    microTuple->SetBranchAddress("tPosJetPtRec",   &tPosJetPtRec);
    microTuple->SetBranchAddress("tPosJetEtaRec",  &tPosJetEtaRec);
    microTuple->SetBranchAddress("tPosJetPhiRec",  &tPosJetPhiRec);

    microTuple->SetBranchAddress("tNegJetMassRec", &tNegJetMassRec);
    microTuple->SetBranchAddress("tNegJetPtRec",   &tNegJetPtRec);
    microTuple->SetBranchAddress("tNegJetEtaRec",  &tNegJetEtaRec);
    microTuple->SetBranchAddress("tNegJetPhiRec",  &tNegJetPhiRec);

//    microTuple->Branch("matchRtlJet",  &matchRtlJet);
//    microTuple->Branch("matchRthJet",  &matchRthJet);
//    microTuple->Branch("matchRwhJet1", &matchRwhJet1);
//    microTuple->Branch("matchRwhJet2", &matchRwhJet2);
}
