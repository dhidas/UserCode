import FWCore.ParameterSet.Config as cms

process = cms.Process("WhatsInANameAnyway")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# The input file name, which won't matter for crab jobs anyway
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      "/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root"
    )
)


# Fill the Dtuple
process.FillDtuple = cms.EDAnalyzer('FillDtuple',
    )

# This will provide the FillDtuple analyzer the output file
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("Dtuple.root")
    )

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('PATLayer1_Output.Dean.root'),
                               # save only events passing the full path
                               #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('') 
                               )
#process.outpath = cms.EndPath(process.out)



# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'])



# require physics declared
process.physDecl = cms.EDFilter("PhysDecl",
    applyfilter = cms.untracked.bool(True)
)

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )


# configure HLT
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )



## global tag for data
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR10_P_V2::All')


process.selectedPatMuons.cut = "pt > 3"
process.selectedPatElectrons.cut = "pt > 3"
process.selectedPatJets.cut = "pt > 3"
process.selectedPatPhotons.cut = "pt > 3"

# sometimes you might need to restrict to only AOD input
#from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD
#restrictInputToAOD(process, ['All'])


process.p = cms.Path(
#    process.hltLevel1GTSeed*
#    process.scrapingVeto*
#    process.physDecl*
    process.primaryVertexFilter*
    process.patDefaultSequence*
    process.FillDtuple
    )

