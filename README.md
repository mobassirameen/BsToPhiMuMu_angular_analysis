# Run2-Bs2PhiMuMu
Run2 Bs --> Phi(KK) mu+ mu- analysis

# How to build:
-------------
CMSSW release For 2016:
For 2018:

$ cmsrel CMSSW_10_2_18


$ cd CMSSW_X_X_XX/src/

$ cmsenv

$ mkdir Run2_analysis

$ cd Run2_analysis/

$ git clone https://gitlab.cern.ch/rraturi/bstophimumu_angular_analysis.git

$ scram b -j8


# How to run:
-----------
$ cd Run2_analysis/BsToPhiMuMu/python/

Recommendation: 
Before running over MC or Data for any period, please do check the GT is correct or not
inside the config file.
Latest GT should be taken.
https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable

JSON file 2016: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON_MuonPhys.txt

JSON file 2017: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt

JSON file 2018: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt


For 2016:

For MC Signal,

$ cmsRun bstophimumu_2016Official_signalMC.py

For JpsiPhi MC,

$ cmsRun bstophimumu_2016Official_JpsiPhiMC.py

For PsiPrimePhi MC,

$ cmsRun bstophimumu_2016Official_PsiprimeMC.py

For 2016 Charmonium dataset

$ cmsRun bstophimumu_2016Data_Charmonium.py

For 2016 DoubleMuonLowmass dataset,

$ cmsRun bstophimumu_2016Data_DoubleMuon.py

For 2017:

For MC Signal,

$ cmsRun bstophimumu_official2017_MC_cfg.py

For JpsiPhi MC,

$ cmsRun bstophimumu_official2017_JpsiPhiMC_cfg.py

For PsiPrimePhi MC,

$ cmsRun bstophimumu_official2017_PsiprimeMC_cfg.py

For 2017 Charmonium dataset

$ cmsRun bstophimumu_2017Data_Charmonium_cfg.py

For 2017 DoubleMuonLowmass dataset,

$ cmsRun bstophimumu_2017Data_DoubleMuon_cfg.py

For 2018:

For MC Signal,

$ cmsRun bstophimumu_official2018_MC_cfg.py

For JpsiPhi MC,

$ cmsRun bstophimumu_official2018_JpsiPhiMC_cfg.py

For PsiPrimePhi MC,

$ cmsRun bstophimumu_official2018_PsiprimeiMC_cfg.py

For 2018 Charmonium dataset

$ cmsRun bstophimumu_2018Data_Charmonium_cfg.py

For 2018 DoubleMuonLowmass dataset,

$ cmsRun bstophimumu_2018Data_DoubleMuon_cfg.py

Example of crab job:

$ crab submit -c crabConfig_bstophimumu_2017B_data_Charmonium.py

For submitting crab job, you need pileup histogram and these are hardcoded
in side the src/BsToPhiMuMu.cc as well crabConfig_bstophimumu_2017B_data_Charmonium.py.

For 2916, 2017, and 2018 pileup files, 

DataPileupHistogram2016_rereco.root
DataPileupHistogram2017_rereco.root
DataPileupHistogram2018_rereco.root

PileupMC_2016.root
PileupMC_2017.root
PileupMC_2018.root


