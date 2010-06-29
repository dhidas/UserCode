import FWCore.ParameterSet.Config as cms
#Runs pat and saves the ntuple on 354_patch1 for data
process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



# Input your file or list of files (separate by ",")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/tmp/dhidas/RECO_3X/Wgamma/F295BB2E-0647-DF11-881C-0030487D858D.root'
        'file:/cms/data16/dhidas/RECO_3X/Wgamma/F295BB2E-0647-DF11-881C-0030487D858D.root'
    )
)


#Maximum number of events, -1=all
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### Geometry and Detector Conditions

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#Give the global tag. For MC I think it is MX_3XY_V20, and data is (I think) START3X_V20:All
process.GlobalTag.globaltag = cms.string('MC_3XY_V20::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")




from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
removeMCMatching(process,['All'])


MyEventContent = [
  'keep *_selectedPatPhotons_*_*',
  'keep *_selectedPatElectrons_*_*',
  'keep *_selectedPatMuons_*_*',
  'keep *_selectedPatJets*_*_*',
  'keep *_patMETs*_*_*',
  'keep *_selectedPatTrackCands_*_*'
]

process.patElectrons.embedTrack = True
process.patElectrons.embedGsfTrack = True
process.patElectrons.embedSuperCluster =  True
process.patElectrons.embedPFCandidate = False

#embedCaloTowers =  True


# Output module configuration
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
              fileName = cms.untracked.string('PATNtuple.root'),
              # save only events passing the full path
              #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
              # save PAT Layer 1 output
              # you need a '*' to unpack the list of commands 'patEventContent'
              outputCommands = cms.untracked.vstring('drop *', *MyEventContent )
)



process.outpath = cms.EndPath(process.patDefaultSequence+process.out)







