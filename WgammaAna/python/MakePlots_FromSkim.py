import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:PATNtuple.root'
    )
)



process.PlotElectronVariables = cms.EDAnalyzer('PlotElectronVariables',
    electronSrc = cms.untracked.InputTag('selectedPatElectrons')
)

process.PlotPhotonVariables = cms.EDAnalyzer('PlotPhotonVariables',
    photonSrc = cms.untracked.InputTag('selectedPatPhotons')
)

# This will provide the FillDtuple analyzer the output file
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("Plots.root")
)



process.p = cms.Path(process.PlotElectronVariables+process.PlotPhotonVariables)
