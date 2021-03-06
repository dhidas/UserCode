import FWCore.ParameterSet.Config as cms



source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_10_1_WxU.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_11_2_0b7.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_12_1_qMT.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_13_1_nng.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_14_1_iNI.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_15_1_HKS.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_16_1_nZX.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_17_1_bbH.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_18_1_42o.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_19_2_CrL.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_1_1_q2j.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_20_1_hga.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_21_2_eKs.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_22_1_An2.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_23_3_jPi.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_24_1_La0.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_24_3_1NF.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_25_1_ES7.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_26_1_ngx.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_27_1_JeZ.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_28_1_kF1.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_29_1_jHL.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_2_1_G99.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_30_0_0PP.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_31_0_tAF.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_32_0_9mT.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_3_1_uNd.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_4_1_VSr.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_5_1_fP2.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_6_1_XZT.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_7_1_Br7.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_8_1_GhH.root',
                             'file:/cms/se/store/user/duggan/Collisions11/PreTechStop/SemiOfficial_v2/SingleElectron/El27_v1/trigVal_patTuple_9_1_FVe.root'
                                                             )                              )