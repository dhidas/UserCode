import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


options = VarParsing.VarParsing("analysis")
options.inputFiles = cms.untracked.vstring('')
options.outputFile = cms.untracked.string('')
options.register ('IsData',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "IsData: True or False [1,0]")
options.parseArguments()


process = cms.Process("Skim")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  input = cms.untracked.int32(-1) 
)
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    options.inputFiles
    ), 
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

IsData = False
if options.IsData == 1:
  print "This IS Data"
  IsData = True

process.patfilter = cms.EDFilter('Filter',
  debug = cms.untracked.bool(False), 
#  OutFileName = cms.untracked.string(''),
  IsData = cms.untracked.bool(IsData),
  JSONFilename = cms.untracked.string('json/Cert_160404-166502_7TeV_PromptReco_Collisions11_JSON.txt'),
  TriggerNames = cms.untracked.vstring(
    'HLT_Mu12_v1',
    'HLT_IsoMu12_v1',
    'HLT_Mu15_v2',
    'HLT_IsoMu12_v1',
    'HLT_IsoMu12_v1',
    'HLT_Mu24_v2',
    'HLT_IsoMu17_v6',
    'HLT_IsoMu17_v8',
    'HLT_IsoMu24_v5',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1'
    )
)

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string( options.outputFile ), 
    # save only events passing the full path 
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ), 
    # save PAT Layer 1 output 
    outputCommands = cms.untracked.vstring('keep *' ) # you need a '*' to unpack the list of commands 'patEventContent' 
    ) 
process.outpath = cms.EndPath(process.out)


process.p = cms.Path(process.patfilter)
