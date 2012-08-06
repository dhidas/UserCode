import FWCore.ParameterSet.Config as cms

process = cms.Process("ttbar")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  #input = cms.untracked.int32(5000)
  input = cms.untracked.int32(-1) 
)
process.load("MyAnalysis.MyPatLayer2s.DeanTest_cfi");
##process.source = cms.Source("PoolSource",
##  skipEvents = cms.untracked.uint32(0), 
##  fileNames = cms.untracked.vstring(
##    '',
##    ''
##    ), 
##  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
##)

process.ttbar = cms.EDAnalyzer('myPatAna',
  sumPtNJetsMin = cms.untracked.double(200.0), 
  debug = cms.untracked.bool(False), 
  outputFilename = cms.untracked.string('File1.root'),
  outputFilename2= cms.untracked.string('File2.root'),
  PatJetType     = cms.untracked.string('selectedPatJetsAK5PF')
)                               

process.p = cms.Path(process.ttbar)
