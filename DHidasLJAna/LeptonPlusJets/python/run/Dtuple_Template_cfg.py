import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


options = VarParsing.VarParsing("analysis")
options.inputFiles = cms.untracked.vstring('')
options.outputFile = cms.untracked.string('')
options.parseArguments()


process = cms.Process("DtupleMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  input = cms.untracked.int32(-1) 
)
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    options.inputFiles
    ), 
)

process.patana = cms.EDAnalyzer('DHidasPatAna',
  debug = cms.untracked.bool(False), 
  #OutFileName = cms.untracked.string('Test_El27_vX_Dtuple.root'),
  OutFileName = cms.untracked.string(options.outputFile),
  IsData = cms.untracked.bool(False),
  MakeDtuple = cms.untracked.bool(True)
)                               



process.p = cms.Path(process.patana)
