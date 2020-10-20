# Run2-Bs2PhiMuMu
Run2 Bs --> Phi(KK) mu+ mu- analysis

We run Selector code over the ntuple produced from MINIAOD/MINIAODSIM samples
Which are present now at root://cmseos.fnal.gov//eos/uscms/store/user/ckar/

Inside selector code, we chose the single candidate depending upon the best Bs selection.

BDT is implemented inside the MVAnalysis class. The xml file need to be change inside that.

To run over the data, 
$ sh process_singlecandtree_data.sh
To run over the MC,
$ sh process_singlecandtree_phimm.sh