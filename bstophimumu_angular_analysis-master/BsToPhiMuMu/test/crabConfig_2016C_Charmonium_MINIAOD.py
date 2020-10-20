from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'job_crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_16C_v2'
config.General.workArea = 'crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_16C_v2'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Charmonium_2016B_dataMINI.py'
config.JobType.outputFiles = ['BsToPhiMuMu_2016_Data.root']

config.JobType.inputFiles = ['PileupMC_2016.root','DataPileupHistogram2016_rereco.root']
config.Data.inputDataset = '/Charmonium/Run2016C-17Jul2018-v1/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10#10

config.JobType.allowUndistributedCMSSW = True
config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_CH_*', 'T2_UK_*', 'T2_IT_*', 'T2_US_*']

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON_MuonPhys.txt'

#config.Data.runRange = '272007-275376'# Era B                                                                                                                               
config.Data.runRange = '275657-276283'# Era C                                                                                                                                
#config.Data.runRange = '276315-276811' #Era D                                                                                                                                
#config.Data.runRange = '276831-277420' #Era E                                                                                                                                
#config.Data.runRange  = '277772-278808' #Era F                                                                                                                               
#config.Data.runRange = '278820-280385' #Era G                                                                                                                                
#config.Data.runRange = '280919-284044' #Era H                                                                                                                               

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_BMM4_crabtest_Charmonium_Run2015B'

config.Site.storageSite = 'T2_IN_TIFR'  # Your output destination. Useful to us: T2_CH_CSCS, T3_CH_PSI, T2_US_Nebraska

