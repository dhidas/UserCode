import FWCore.ParameterSet.Config as cms

process = cms.Process("ttbar")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
  input = cms.untracked.int32(-1)
#input = cms.untracked.int32(3000) 
)
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring( 
<<<<<<< Test.py
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v2/trigVal_patTuple_3_2_PaP.root'
#                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleMu/Mu15_v2/trigVal_patTuple_48_1_UW8.root'
=======
#                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleMu/Mu15_v2/trigVal_patTuple_38_1_7bm.root'
                             'file:/cms/dan/Data/Collisions11/v2/DCSOnly/SingleMu/trigVal_patTuple_230_2_HEZ.root'
>>>>>>> 1.5
                                                             )
                              )
##process.source = cms.Source("PoolSource",
##  skipEvents = cms.untracked.uint32(0), 
##  fileNames = cms.untracked.vstring(
##    '',
##    ''
##    ), 
##  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
##)



process.patana = cms.EDAnalyzer('DHidasPatAna',
  debug = cms.untracked.bool(False), 
  OutFileName = cms.untracked.string('Test.root'),
  MakeDtuple = cms.untracked.bool(True),
  IsData = cms.untracked.bool(False)
)                               

process.p = cms.Path(process.patana)
