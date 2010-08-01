#include "VgammaAna/TaiwanAna/interface/VgNtuple.h"

VgNtuple::VgNtuple ()
{
  tree_ = new TChain("VgAnalyzerKit/EventTree", "VgAnalyzerKit/EventTree");
}

VgNtuple::~VgNtuple ()
{
}


void VgNtuple::AddInFiles (std::vector<TString> const& InFiles)
{
  for (size_t i = 0; i != InFiles.size(); ++i) {
    tree_->Add(InFiles[i]);
  }

  SetBranchAddresses();

  return;
}


void VgNtuple::SetBranchAddresses ()
{
  bool doGenParticles_ = false;

  tree_->SetBranchAddress("run", &run_);
  tree_->SetBranchAddress("event", &event_);
  tree_->SetBranchAddress("orbit", &orbit_);
  tree_->SetBranchAddress("bx", &bx_);
  tree_->SetBranchAddress("lumis", &lumis_);
  tree_->SetBranchAddress("isData", &isData_);
  tree_->SetBranchAddress("ttbit", ttbit_);
  tree_->SetBranchAddress("nHLT", &nHLT_);
  tree_->SetBranchAddress("HLT", HLT_);
  tree_->SetBranchAddress("nHFTowersP", &nHFTowersP_);
  tree_->SetBranchAddress("nHFTowersN", &nHFTowersN_);
  tree_->SetBranchAddress("nVtx", &nVtx_);
  tree_->SetBranchAddress("vtx", vtx_);
  tree_->SetBranchAddress("vtxNTrk", vtxNTrk_);
  tree_->SetBranchAddress("vtxNDF", vtxNDF_);
  tree_->SetBranchAddress("vtxD0", vtxD0_);
  tree_->SetBranchAddress("IsVtxGood", &IsVtxGood_);
  tree_->SetBranchAddress("nTrk", &nTrk_);
  tree_->SetBranchAddress("nGoodTrk", &nGoodTrk_);
  tree_->SetBranchAddress("IsTracksGood", &IsTracksGood_);
  if (doGenParticles_) {
    tree_->SetBranchAddress("pdf", pdf_);
    tree_->SetBranchAddress("processID", &processID_);
    // genParticle
    tree_->SetBranchAddress("nMC", &nMC_);
    tree_->SetBranchAddress("mcPID", mcPID);
    tree_->SetBranchAddress("mcPt", mcPt);
    tree_->SetBranchAddress("mcMass", mcMass);
    tree_->SetBranchAddress("mcEta", mcEta);
    tree_->SetBranchAddress("mcPhi", mcPhi);
    tree_->SetBranchAddress("mcGMomPID", mcGMomPID);
    tree_->SetBranchAddress("mcMomPID", mcMomPID);
    tree_->SetBranchAddress("mcMomPt", mcMomPt);
    tree_->SetBranchAddress("mcMomMass", mcMomMass);
    tree_->SetBranchAddress("mcMomEta", mcMomEta);
    tree_->SetBranchAddress("mcMomPhi", mcMomPhi);
    tree_->SetBranchAddress("mcIndex", mcIndex);
    tree_->SetBranchAddress("mcDecayType", mcDecayType); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
    // Gen MET
    tree_->SetBranchAddress("genMET", &genMET_);
    tree_->SetBranchAddress("genMETx", &genMETx_);
    tree_->SetBranchAddress("genMETy", &genMETy_);
    tree_->SetBranchAddress("genMETPhi", &genMETPhi_);
  }
  // Calo MET
  tree_->SetBranchAddress("MET", &MET_);
  tree_->SetBranchAddress("METx", &METx_);
  tree_->SetBranchAddress("METy", &METy_);
  tree_->SetBranchAddress("METPhi", &METPhi_);
  tree_->SetBranchAddress("METsumEt", &METsumEt_);
  tree_->SetBranchAddress("uncorrMET", uncorrMET_); // [0]: uncorrALL, [1]: uncorrJES, [2]: uncorrMUON
  tree_->SetBranchAddress("uncorrMETPhi", uncorrMETPhi_);
  tree_->SetBranchAddress("uncorrMETSumEt", uncorrMETSumEt_);
  // tcMET
  tree_->SetBranchAddress("tcMET", &tcMET_);
  tree_->SetBranchAddress("tcMETx", &tcMETx_);
  tree_->SetBranchAddress("tcMETy", &tcMETy_);
  tree_->SetBranchAddress("tcMETPhi", &tcMETPhi_);
  tree_->SetBranchAddress("tcMETsumEt", &tcMETsumEt_);
  tree_->SetBranchAddress("tcMETmEtSig", &tcMETmEtSig_);
  tree_->SetBranchAddress("tcMETSig", &tcMETSig_);
  // pfMET
  tree_->SetBranchAddress("pfMET", &pfMET_);
  tree_->SetBranchAddress("pfMETx", &pfMETx_);
  tree_->SetBranchAddress("pfMETy", &pfMETy_);
  tree_->SetBranchAddress("pfMETPhi", &pfMETPhi_);
  tree_->SetBranchAddress("pfMETsumEt", &pfMETsumEt_);
  tree_->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig_);
  tree_->SetBranchAddress("pfMETSig", &pfMETSig_);
  // Electron
  tree_->SetBranchAddress("nEle", &nEle_);
  tree_->SetBranchAddress("eleID", eleID_); // [0]: eidRobustLoose, [1]: eidRobustTight, [2]: eidLoose, [3]: eidTight, [4]: eidRobustHighEnergy
  tree_->SetBranchAddress("eleClass", eleClass_);
  tree_->SetBranchAddress("eleCharge", eleCharge_);
  tree_->SetBranchAddress("eleEn", eleEn_);
  tree_->SetBranchAddress("eleSCRawEn", eleSCRawEn_);
  tree_->SetBranchAddress("eleESEn", eleESEn_);
  tree_->SetBranchAddress("eleSCEn", eleSCEn_);
  tree_->SetBranchAddress("elePt", elePt_);
  tree_->SetBranchAddress("elePz", elePz_);
  tree_->SetBranchAddress("eleEta", eleEta_);
  tree_->SetBranchAddress("elePhi", elePhi_);
  tree_->SetBranchAddress("eleSCEta", eleSCEta_);
  tree_->SetBranchAddress("eleSCPhi", eleSCPhi_);
  tree_->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth_);
  tree_->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth_);
  tree_->SetBranchAddress("eleVtx", eleVtx_);
  tree_->SetBranchAddress("eleCaloPos", eleCaloPos_);
  tree_->SetBranchAddress("eleSCPos", eleSCPos_);
  tree_->SetBranchAddress("eleHoverE", eleHoverE_);
  tree_->SetBranchAddress("eleHoverE1", eleHoverE1_);
  tree_->SetBranchAddress("eleHoverE2", eleHoverE2_);
  tree_->SetBranchAddress("eleEoverP", eleEoverP_);
  tree_->SetBranchAddress("elePin", elePin_);
  tree_->SetBranchAddress("elePout", elePout_);
  tree_->SetBranchAddress("eleBrem", eleBrem_);
  tree_->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx_);
  tree_->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx_);
  tree_->SetBranchAddress("eleSigmaEtaEta", eleSigmaEtaEta_);
  tree_->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta_);
  tree_->SetBranchAddress("eleEMax", eleEMax_);
  tree_->SetBranchAddress("eleE2nd", eleE2nd_);
  tree_->SetBranchAddress("eleE2x2", eleE2x2_);
  tree_->SetBranchAddress("eleE3x2", eleE3x2_);
  tree_->SetBranchAddress("eleE3x3", eleE3x3_);
  tree_->SetBranchAddress("eleE4x4", eleE4x4_);
  tree_->SetBranchAddress("eleE5x5", eleE5x5_);
  tree_->SetBranchAddress("eleE2x5Right", eleE2x5Right_);
  tree_->SetBranchAddress("eleE2x5Left", eleE2x5Left_);
  tree_->SetBranchAddress("eleE2x5Top", eleE2x5Top_);
  tree_->SetBranchAddress("eleE2x5Bottom", eleE2x5Bottom_);
  tree_->SetBranchAddress("eleERight", eleERight_);
  tree_->SetBranchAddress("eleELeft", eleELeft_);
  tree_->SetBranchAddress("eleETop", eleETop_);
  tree_->SetBranchAddress("eleEBottom", eleEBottom_);
  if (doGenParticles_) {
    tree_->SetBranchAddress("eleGenIndex", eleGenIndex_);
    tree_->SetBranchAddress("eleGenGMomPID", eleGenGMomPID_);
    tree_->SetBranchAddress("eleGenMomPID", eleGenMomPID_);
    tree_->SetBranchAddress("eleGenMomPt", eleGenMomPt_);
  }
  tree_->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03_);
  tree_->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03_);
  tree_->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03_);
  tree_->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04_);
  tree_->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04_);
  tree_->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04_);
  tree_->SetBranchAddress("eleChi2NDF", eleChi2NDF_);
  tree_->SetBranchAddress("eleD0", eleD0_);
  tree_->SetBranchAddress("eleNumberOfValidHits", eleNumberOfValidHits_);
  // Photon
  tree_->SetBranchAddress("nPho", &nPho_);
  tree_->SetBranchAddress("phoIsPhoton", phoIsPhoton_);
  tree_->SetBranchAddress("phoE", phoE_);
  tree_->SetBranchAddress("phoEt", phoEt_);
  tree_->SetBranchAddress("phoPz", phoPz_);
  tree_->SetBranchAddress("phoEta", phoEta_);
  tree_->SetBranchAddress("phoPhi", phoPhi_);
  tree_->SetBranchAddress("phoR9", phoR9_);
  tree_->SetBranchAddress("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03_);
  tree_->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03_);
  tree_->SetBranchAddress("phoNTrkSolidDR03", phoNTrkSolidDR03_);
  tree_->SetBranchAddress("phoNTrkHollowDR03", phoNTrkHollowDR03_);
  tree_->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03_);
  tree_->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03_);
  tree_->SetBranchAddress("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04_);
  tree_->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04_);
  tree_->SetBranchAddress("phoNTrkSolidDR04", phoNTrkSolidDR04_);
  tree_->SetBranchAddress("phoNTrkHollowDR04", phoNTrkHollowDR04_);
  tree_->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04_);
  tree_->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04_);
  tree_->SetBranchAddress("phoHoverE", phoHoverE_);
  tree_->SetBranchAddress("phoHoverE1", phoHoverE1_);
  tree_->SetBranchAddress("phoHoverE2", phoHoverE2_);
  tree_->SetBranchAddress("phoSigmaEtaEta", phoSigmaEtaEta_);
  tree_->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta_);
  tree_->SetBranchAddress("phoSeedTime", phoSeedTime_);
  tree_->SetBranchAddress("phoPos", phoPos_);
  tree_->SetBranchAddress("phoEMax", phoEMax_);
  tree_->SetBranchAddress("phoE2nd", phoE2nd_);
  tree_->SetBranchAddress("phoE2x2", phoE2x2_);
  tree_->SetBranchAddress("phoE3x2", phoE3x2_);
  tree_->SetBranchAddress("phoE3x3", phoE3x3_);
  tree_->SetBranchAddress("phoE4x4", phoE4x4_);
  tree_->SetBranchAddress("phoE5x5", phoE5x5_);
  tree_->SetBranchAddress("phoE2x5Right", phoE2x5Right_);
  tree_->SetBranchAddress("phoE2x5Left", phoE2x5Left_);
  tree_->SetBranchAddress("phoE2x5Top", phoE2x5Top_);
  tree_->SetBranchAddress("phoE2x5Bottom", phoE2x5Bottom_);
  tree_->SetBranchAddress("phoERight", phoERight_);
  tree_->SetBranchAddress("phoELeft", phoELeft_);
  tree_->SetBranchAddress("phoETop", phoETop_);
  tree_->SetBranchAddress("phoEBottom", phoEBottom_);
  if (doGenParticles_) {
    tree_->SetBranchAddress("phoGenIndex", phoGenIndex_);
    tree_->SetBranchAddress("phoGenGMomPID", phoGenGMomPID);
    tree_->SetBranchAddress("phoGenMomPID", phoGenMomPID);
    tree_->SetBranchAddress("phoGenMomPt", phoGenMomPt);
  }
  tree_->SetBranchAddress("phoSCE", phoSCE_);
  tree_->SetBranchAddress("phoSCEt", phoSCEt_);
  tree_->SetBranchAddress("phoSCEta", phoSCEta_);
  tree_->SetBranchAddress("phoSCPhi", phoSCPhi_);
  tree_->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth_);
  tree_->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth_);
  tree_->SetBranchAddress("phoOverlap", phoOverlap_);
  // Muon
  tree_->SetBranchAddress("nMu", &nMu_);
  tree_->SetBranchAddress("muEta", muEta_);
  tree_->SetBranchAddress("muPhi", muPhi_);
  tree_->SetBranchAddress("muCharge", muCharge_);
  tree_->SetBranchAddress("muPt", muPt_);
  tree_->SetBranchAddress("muPz", muPz_);
  if (doGenParticles_)
    tree_->SetBranchAddress("muGenIndex", muGenIndex_);
  tree_->SetBranchAddress("muIsoTrk", muIsoTrk_);
  tree_->SetBranchAddress("muIsoCalo", muIsoCalo_);
  tree_->SetBranchAddress("muIsoEcal", muIsoEcal_);
  tree_->SetBranchAddress("muIsoHcal", muIsoHcal_);
  tree_->SetBranchAddress("muEmVeto", muEmVeto_);
  tree_->SetBranchAddress("muHadVeto", muHadVeto_);
  tree_->SetBranchAddress("muType", muType_);
  tree_->SetBranchAddress("muID", muID_);
  // [0]: AllArbitrated, [1]: GlobalMuonPromptTight, [2]: TMLSLoose, [3]: TMLSTight, [4]: TM2DCompatLoose, [5]: TM2DCompatTight
  tree_->SetBranchAddress("muD0", muD0_);
  tree_->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits_);
  // Jet
  tree_->SetBranchAddress("nJet", &nJet_);
  tree_->SetBranchAddress("jetEn", jetEn_);
  tree_->SetBranchAddress("jetPt", jetPt_);
  tree_->SetBranchAddress("jetEta", jetEta_);
  tree_->SetBranchAddress("jetPhi", jetPhi_);
  tree_->SetBranchAddress("jetMass", jetMass_);
  tree_->SetBranchAddress("jetEt", jetEt_);
  tree_->SetBranchAddress("jetmaxEInEmTowers", jetmaxEInEmTowers_);
  tree_->SetBranchAddress("jetmaxEInHadTowers", jetmaxEInHadTowers_);
  tree_->SetBranchAddress("jetenergyFractionHadronic", jetenergyFractionHadronic_);
  tree_->SetBranchAddress("jetemEnergyFraction", jetemEnergyFraction_);
  if (doGenParticles_) {
    tree_->SetBranchAddress("jetGenIndex", jetGenIndex_);
    tree_->SetBranchAddress("jetGenJetIndex", jetGenJetIndex_);
    tree_->SetBranchAddress("jetGenJetEn", jetGenJetEn_);
    tree_->SetBranchAddress("jetGenJetPt", jetGenJetPt_);
    tree_->SetBranchAddress("jetGenJetEta", jetGenJetEta_);
    tree_->SetBranchAddress("jetGenJetPhi", jetGenJetPhi_);
    tree_->SetBranchAddress("jetGenJetMass", jetGenJetMass_);
    tree_->SetBranchAddress("jetGenPartonID", jetGenPartonID_);
    tree_->SetBranchAddress("jetGenPartonMomID", jetGenPartonMomID_);
    tree_->SetBranchAddress("jetGenEn", jetGenEn_);
    tree_->SetBranchAddress("jetGenPt", jetGenPt_);
    tree_->SetBranchAddress("jetGenEta", jetGenEta_);
    tree_->SetBranchAddress("jetGenPhi", jetGenPhi_);
  }
  tree_->SetBranchAddress("jetpartonFlavour", jetpartonFlavour_);
  tree_->SetBranchAddress("jetRawPt", jetRawPt_);
  tree_->SetBranchAddress("jetRawEn", jetRawEn_);
  tree_->SetBranchAddress("jetCharge", jetCharge_);
  tree_->SetBranchAddress("jetCombinedSVBJetTags", jetCombinedSVBJetTags_);
  tree_->SetBranchAddress("jetCombinedSVMVABJetTags", jetCombinedSVMVABJetTags_);
  tree_->SetBranchAddress("jetConeIsoTauJetTags", jetConeIsoTauJetTags_);
  tree_->SetBranchAddress("jetImpactParaMVABJetTags", jetImpactParaMVABJetTags_);
  tree_->SetBranchAddress("jetJetBProbBJetTags", jetJetBProbBJetTags_);
  tree_->SetBranchAddress("jetJetProbBJetTags", jetJetProbBJetTags_);
  tree_->SetBranchAddress("jetSimpleSVBJetTags", jetSimpleSVBJetTags_);
  tree_->SetBranchAddress("jetSoftElecBJetTags", jetSoftElecBJetTags_);
  tree_->SetBranchAddress("jetSoftMuonBJetTags", jetSoftMuonBJetTags_);
  tree_->SetBranchAddress("jetSoftMuonNoIPBJetTags", jetSoftMuonNoIPBJetTags_);
  tree_->SetBranchAddress("jetTrackCountHiEffBJetTags", jetTrackCountHiEffBJetTags_);
  tree_->SetBranchAddress("jetTrackCountHiPurBJetTags", jetTrackCountHiPurBJetTags_);
  tree_->SetBranchAddress("jetJetLRval", jetJetLRval_);
  tree_->SetBranchAddress("jetJetProb", jetJetProb_);
  // Zee candiate
  tree_->SetBranchAddress("nZee", &nZee_);
  tree_->SetBranchAddress("ZeeMass", ZeeMass_);
  tree_->SetBranchAddress("ZeePt", ZeePt_);
  tree_->SetBranchAddress("ZeeEta", ZeeEta_);
  tree_->SetBranchAddress("ZeePhi", ZeePhi_);
  tree_->SetBranchAddress("ZeeLeg1Index", ZeeLeg1Index_);
  tree_->SetBranchAddress("ZeeLeg2Index", ZeeLeg2Index_);
  // Zmumu candiate
  tree_->SetBranchAddress("nZmumu", &nZmumu_);
  tree_->SetBranchAddress("ZmumuMass", ZmumuMass_);
  tree_->SetBranchAddress("ZmumuPt", ZmumuPt_);
  tree_->SetBranchAddress("ZmumuEta", ZmumuEta_);
  tree_->SetBranchAddress("ZmumuPhi", ZmumuPhi_);
  tree_->SetBranchAddress("ZmumuLeg1Index", ZmumuLeg1Index_);
  tree_->SetBranchAddress("ZmumuLeg2Index", ZmumuLeg2Index_);
  // Wenu candidate
  tree_->SetBranchAddress("nWenu", &nWenu_);
  tree_->SetBranchAddress("WenuMassTCaloMET", WenuMassTCaloMET_);
  tree_->SetBranchAddress("WenuEtCaloMET", WenuEtCaloMET_);
  tree_->SetBranchAddress("WenuACopCaloMET", WenuACopCaloMET_);
  tree_->SetBranchAddress("WenuMassTTcMET", WenuMassTTcMET_);
  tree_->SetBranchAddress("WenuEtTcMET", WenuEtTcMET_);
  tree_->SetBranchAddress("WenuACopTcMET", WenuACopTcMET_);
  tree_->SetBranchAddress("WenuMassTPfMET", WenuMassTPfMET_);
  tree_->SetBranchAddress("WenuEtPfMET", WenuEtPfMET_);
  tree_->SetBranchAddress("WenuACopPfMET", WenuACopPfMET_);
  // Wmunu candidate
  tree_->SetBranchAddress("nWmunu", &nWmunu_);
  tree_->SetBranchAddress("WmunuMassTCaloMET", WmunuMassTCaloMET_);
  tree_->SetBranchAddress("WmunuEtCaloMET", WmunuEtCaloMET_);
  tree_->SetBranchAddress("WmunuACopCaloMET", WmunuACopCaloMET_);
  tree_->SetBranchAddress("WmunuMassTTcMET", WmunuMassTTcMET_);
  tree_->SetBranchAddress("WmunuEtTcMET", WmunuEtTcMET_);
  tree_->SetBranchAddress("WmunuACopTcMET", WmunuACopTcMET_);
  tree_->SetBranchAddress("WmunuMassTPfMET", WmunuMassTPfMET_);
  tree_->SetBranchAddress("WmunuEtPfMET", WmunuEtPfMET_);
  tree_->SetBranchAddress("WmunuACopPfMET", WmunuACopPfMET_);



  return;
}
