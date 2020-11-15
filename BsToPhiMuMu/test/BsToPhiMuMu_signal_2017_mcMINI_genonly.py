import FWCore.ParameterSet.Config as cms
#import os


#FILE1 = os.environ.get('FILE1')
#FILE2 = os.environ.get('FILE2')


process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

RunOnMC = True
KeepGen = True

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = cms.string('94X_dataRun2_v10') ## for 2016
# process.GlobalTag.globaltag = cms.string('94X_dataRun2_v11') ## for 2017
#process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v20')
process.GlobalTag.globaltag = cms.string('102X_mcRun2_asymptotic_v7')

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring(FILE1))
                            fileNames = cms.untracked.vstring('file:/uscms/home/ckar/nobackup/PYTHIA8-BsToPhiMuMu_2017-XXXX_nofilter_step4_100.root'))
#   '/store/mc/RunIIAutumn18MiniAOD/BsToMuMuPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/40000/C8C2B108-2643-A448-97E7-7EB799346E5C.root'))
        #'/store/data/Run2016G/Charmonium/MINIAOD/17Jul018-v1/40000/FE986973-338D-E811-A968-0242AC130002.root',
        #'file:/eos/home-c/ckar/BSTOPHIMUMU/files/bsTophimumu_mc_2016.root', 
        #'file:/eos/home-c/ckar/BSTOPHIMUMU/files/bsTjpsiophi_mc_2016.root',
        #'file:/eos/home-c/ckar/BSTOPHIMUMU/files/charmonium_2018.root',
        #'file:/eos/home-c/ckar/BSTOPHIMUMU/files/charmonium_2017.root',
#     ),
# )

taskB0 = cms.Task()

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
taskB0.add(process.unpackedTracksAndVertices)

g_TriggerNames_LastFilterNames = [
    ('HLT_DoubleMu4_JpsiTrk_Displaced',  'hltDisplacedmumuFilterDoubleMu4Jpsi'),
    ('HLT_DoubleMu4_PsiPrimeTrk_Displaced', 'hltDisplacedmumuFilterDoubleMu4PsiPrime'),
    ('HLT_DoubleMu4_LowMassNonResonantTrk_Displaced',  'hltLowMassNonResonantTkVertexFilter')
]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]


#process.load("myanalyzer.BsToPhiMuMu.slimmedMuonsTriggerMatcher_cfi")  
process.load("slimmedMuonsTriggerMatcher_cfi")  
process.ntuple = cms.EDAnalyzer('BsToPhiMuMu',
                                OutputFileName = cms.string("BsToPhiMuMu_2018_SignalMC_Mini_unfilter.root"),
                                #OutputFileName = cms.string(FILE2),
                                BuildBsToPhiMuMu = cms.untracked.bool(True),
                                
                                MuonMass = cms.untracked.double(0.10565837),
                                MuonMassErr = cms.untracked.double(3.5e-9),
                                KaonMass = cms.untracked.double(0.493677),
                                KaonMassErr = cms.untracked.double(1.6e-5),
                                BsMass = cms.untracked.double(5.36677),          ## put the Bs Mass (pdg value)

                                pruned = cms.InputTag("prunedGenParticles"),
                                packed = cms.InputTag("packedGenParticles"),
                                
                                TriggerResultsLabel = cms.InputTag("TriggerResults","", "HLT"),
                                prescales     = cms.InputTag("patTrigger"),
                                objects       = cms.InputTag("slimmedPatTrigger"),
                                TriggerNames  = cms.vstring("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v", 
                                                            "HLT_DoubleMu4_JpsiTrk_Displaced_v",
                                                            "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"
                                                        ),
                                
                                BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                VertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                #MuonLabel = cms.InputTag("slimmedMuons"),
                                MuonLabel = cms.InputTag('slimmedMuonsWithTrigger'),
                                TrackLabel = cms.InputTag("unpackedTracksAndVertices"),
                                
                                
                                PuInfoTag        = cms.InputTag("slimmedAddPileupInfo"),
                                #TriggerNames = cms.vstring([]),
                                LastFilterNames = cms.vstring([]),
                                

                                IsMonteCarlo = cms.untracked.bool(RunOnMC),
                                KeepGENOnly  = cms.untracked.bool(KeepGen),

                                TruthMatchMuonMaxR = cms.untracked.double(0.004), # [eta-phi]
                                TruthMatchKaonMaxR = cms.untracked.double(0.3), # [eta-phi]

                                MuonMinPt = cms.untracked.double(4.0), # 3.0 [GeV]
                                MuonMaxEta = cms.untracked.double(2.4),
                                MuonMaxDcaBs = cms.untracked.double(2.0), # [cm]
                                MuMuMinPt = cms.untracked.double(6.9),      # [GeV/c]
                                MuMuMinInvMass = cms.untracked.double(1.0), # [GeV/c2]
                                MuMuMaxInvMass = cms.untracked.double(4.8), # [GeV/c2]
                                MuMuMinVtxCl = cms.untracked.double(0.10), # 0.05
                                MuMuMinLxySigmaBs = cms.untracked.double(3.0),
                                MuMuMaxDca = cms.untracked.double(0.5), # [cm]
                                MuMuMinCosAlphaBs = cms.untracked.double(0.9),

                                TrkMinPt = cms.untracked.double(0.8), # 0.4 [GeV/c]
                                TrkMinDcaSigBs = cms.untracked.double(0.8), # 0.8 hadron DCA/sigma w/respect to BS (=>changed Max to Min)
                                TrkMaxR = cms.untracked.double(110.0), # [cm] ==> size of tracker volume in radial direction
                                TrkMaxZ = cms.untracked.double(280.0), # [cm] ==> size of tracker volume in Z direction
                                PhiMinMass = cms.untracked.double(1.00), # [GeV/c2]  - 3 sigma of the width(~5MeV)
                                PhiMaxMass = cms.untracked.double(1.04), # [GeV/c2]  + 3 sigma of the width

                                BsMinVtxCl = cms.untracked.double(0.01),
                                BsMinMass = cms.untracked.double(4.7), # [GeV/c2]
                                BsMaxMass = cms.untracked.double(6.0), # [GeV/c2]
                                
                                #printMsg         = cms.untracked.bool(False)
                               )

# process.TFileService = cms.Service('TFileService', 
#                                    fileName = cms.string(
#                                        'BsToPhiMuMu_miniaod.root'
#                                    ), 
#                                    closeFileFast = cms.untracked.bool(True)
#                                )
#process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)

#process.ntupPath = cms.Path(process.ntuple)
#process.ntupPath.associate(taskB0)


process.p = cms.Path(process.slimmedMuonsWithTriggerSequence * process.unpackedTracksAndVertices * process.ntuple)

