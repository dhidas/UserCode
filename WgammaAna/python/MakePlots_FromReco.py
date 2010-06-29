import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/tmp/dhidas/RECO_3X/Wgamma/F295BB2E-0647-DF11-881C-0030487D858D.root'
        'file:/cms/data16/dhidas/RECO_3X/Wgamma/F295BB2E-0647-DF11-881C-0030487D858D.root'
    )
)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#Give the global tag. For MC I think it is MX_3XY_V20, and data is (I think) START3X_V20:All
process.GlobalTag.globaltag = cms.string('MC_3XY_V20::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent




from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
removeMCMatching(process,['All'])




# This will provide the FillDtuple analyzer the output file
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("Plots.root")
)



process.demo = cms.EDAnalyzer('PlotElectronVariables',
    electronSrc = cms.untracked.InputTag('selectedPatElectrons')
)

process.PlotPhotonVariables = cms.EDAnalyzer('PlotPhotonVariables',
    photonSrc = cms.untracked.InputTag('selectedPatPhotons')
)


process.p = cms.Path(process.patDefaultSequence+process.demo+process.PlotPhotonVariables)
