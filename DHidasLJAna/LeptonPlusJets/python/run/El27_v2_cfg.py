import FWCore.ParameterSet.Config as cms

process = cms.Process("data")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  #input = cms.untracked.int32(5000)
  input = cms.untracked.int32(-1) 
)
process.load("DHidasTopBSM.LeptonPlusJets.El27_v2_cfi");
##process.source = cms.Source("PoolSource",
##  skipEvents = cms.untracked.uint32(0), 
##  fileNames = cms.untracked.vstring(
##    '',
##    ''
##    ), 
##  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
##)

process.patana = cms.EDAnalyzer('DHidasPatAna',
  debug = cms.untracked.bool(False), 
  OutFileName = cms.untracked.string('El27_v2.root'),
  IsData = cms.untracked.bool(False)
)                               

process.p = cms.Path(process.patana)
