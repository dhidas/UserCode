import FWCore.ParameterSet.Config as cms

process = cms.Process("data")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  input = cms.untracked.int32(-1) 
)
process.load("DHidasLJAna.LeptonPlusJets.Mu15_v2_cfi");
##process.source = cms.Source("PoolSource",
##  skipEvents = cms.untracked.uint32(0), 
##  fileNames = cms.untracked.vstring(
##    '',
##    ''
##    ), 
##  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
##)

process.patfilter = cms.EDFilter('Filter',
  debug = cms.untracked.bool(False), 
#  OutFileName = cms.untracked.string(''),
  IsData = cms.untracked.bool(False)
)

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string('/cms/data21/dhidas/TopBSMSkims/Mu15_v2_Skim.root'), 
    # save only events passing the full path 
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ), 
    # save PAT Layer 1 output 
    outputCommands = cms.untracked.vstring('keep *' ) # you need a '*' to unpack the list of commands 'patEventContent' 
    ) 
process.outpath = cms.EndPath(process.out)


process.p = cms.Path(process.patfilter)
