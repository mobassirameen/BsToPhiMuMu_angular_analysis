from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'job_crab_MCOfficial_Bs2JpsiPhi_finaljob_MINIAODSIM_CMSSW10218'
config.General.workArea = 'crab_MCOfficial_Bs2JpsiPhi_finaljob_MINIAODSIM_CMSSW10218'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'BsToPhiMuMu_Jpsi_2016_mcMINI.py'
config.JobType.outputFiles = ['BsToPhiMuMu_Mini2016_JpsiMc.root']

config.JobType.inputFiles = ['PileupMC_2016.root','DataPileupHistogram2016_rereco.root']
config.Data.inputDataset = '/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 5#10

config.JobType.allowUndistributedCMSSW = True
config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_CH_*', 'T2_UK_*', 'T2_IT_*', 'T2_US_*']


config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_BMM4_crabtest_Charmonium_Run2015B'

config.Site.storageSite = 'T3_US_FNALLPC'  # Your output destination. Useful to us: T2_CH_CSCS, T3_CH_PSI, T2_US_Nebraska

