import FWCore.ParameterSet.Config as cms

process = cms.Process("WhatsInANameAnyway")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# The input file name, which won't matter for crab jobs anyway
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      "file://INFILE"
    )
)


# Fill the Dtuple
process.FillDtuple = cms.EDAnalyzer('FillDtuple',
    )

# This will provide the FillDtuple analyzer the output file
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("Dtuple.root")
    )




# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.selectedLayer1Muons.cut = "pt > 8"
process.selectedLayer1Electrons.cut = "pt > 8"
process.selectedLayer1Jets.cut = "pt > 15"
process.selectedLayer1Photons.cut = "pt > 8"

# sometimes you might need to restrict to only AOD input
#from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD
#restrictInputToAOD(process, ['All'])


process.p = cms.Path(
    process.patDefaultSequence
    *process.FillDtuple
    )

