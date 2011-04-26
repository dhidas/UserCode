import FWCore.ParameterSet.Config as cms

process = cms.Process("dataana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  #input = cms.untracked.int32(5000)
  input = cms.untracked.int32(-1) 
)
process.source = cms.Source("PoolSource",
  skipEvents = cms.untracked.uint32(0), 
  fileNames = cms.untracked.vstring(
    'file:/cms/data21/dhidas/TopBSMSkims/El27_vX_Skim.root'
    ), 
  duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
)

process.patana = cms.EDAnalyzer('DHidasPatAna',
  debug = cms.untracked.bool(False), 
  OutFileName = cms.untracked.string('RunOnSkim_El27_vX.root'),
  IsData = cms.untracked.bool(False),
  MakeDtuple = cms.untracked.bool(True)
)                               

process.p = cms.Path(process.patana)
