import FWCore.ParameterSet.Config as cms

process = cms.Process("WhatsInANameAnyway")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# The input file name, which won't matter for crab jobs anyway
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/data/dhidas/B633FC04-5FB2-DE11-BED2-001E0B5FC57A.root'
      )
)


# Fill the Dtuple
process.FillDtuple = cms.EDAnalyzer('FillDtuple',
    OutFileName = cms.untracked.string("")
    )

# This will provide the FillDtuple analyzer the output file
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("Dtuple.root")
    )




# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from PhysicsTools.PatAlgos.selectionLayer1.selectedLayer1Objects_cff import *
#from PhysicsTools.PatAlgos.cleaningLayer1.cleanLayer1Objects_cff import *
#cleanLayer1Objects = cms.Sequence(
#    cleanLayer1Muons *        # NOW WE MUST USE '*' AS THE ORDER MATTERS
#    cleanLayer1Electrons *
#    cleanLayer1Photons *
#    #cleanLayer1Taus *
#    cleanLayer1Jets *
#    #cleanLayer1Hemispheres *
#    cleanLayer1Summary
#    )
#process.cleanLayer1Electrons.preselection = "(electronID('eidRobustHighEnergy') > 0)" # && (trackIso < 3) && pt > 20"
#process.cleanLayer1Muons.preselection = "isGood('GlobalMuonPromptTight')" # && pt > 20"
#process.cleanLayer1Muons.preselection = "pt > 10"
#process.cleanLayer1Electrons.preselection = "pt > 10"
#process.cleanLayer1Jets.preselection = "pt > 10"
#process.cleanLayer1Photons.preselection = "pt > 10"
##process.cleanLayer1Photons.photonIDSource = cms.InputTag("PhotonIDProd","PhotonAssociatedID")
##process.cleanLayer1Photons.photonIDSource = cms.InputTag("PhotonCutBasedIDTight","PhotonAssociatedID")
##process.cleanLayer1Photons.preselection = "photonID('PhotonCutBasedIDTight') > 0"
#
#process.cleanLayer1Electrons.checkOverlaps.muons.requireNoOverlaps = True
##process.cleanLayer1Photons.checkOverlaps.muons.requireNoOverlaps = True
#process.cleanLayer1Photons.checkOverlaps.electrons.requireNoOverlaps = True
#process.cleanLayer1Jets.checkOverlaps.muons.requireNoOverlaps = True
#process.cleanLayer1Jets.checkOverlaps.electrons.requireNoOverlaps = True
#process.cleanLayer1Jets.checkOverlaps.photons.requireNoOverlaps = Truea

process.selectedLayer1Muons.cut = "pt > 8"
process.selectedLayer1Electrons.cut = "pt > 8"
process.selectedLayer1Jets.cut = "pt > 15"
process.selectedLayer1Photons.cut = "pt > 8"


process.p = cms.Path(
    process.patDefaultSequence
    *process.FillDtuple
    )


