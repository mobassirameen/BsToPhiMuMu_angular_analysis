#!/bin/bash                                                                                                                                                                       

echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
echo -e "\n=> processing singlecand ntuple for Bs2phimm sample."

# 2000000


#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016B_MINI_Correction_Iso.root  DoubleMuonLowMass_2016B_Mini_Iso -s 0  -n 2000000
#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016C_MINI_Correction_Iso.root  DoubleMuonLowMass_2016C_Mini_Iso -s 0  -n 2000000
#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016D_MINI_Correction_Iso.root  DoubleMuonLowMass_2016D_Mini_Iso -s 0  -n 2000000
#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016E_MINI_Correction_Iso.root  DoubleMuonLowMass_2016E_Mini_Iso -s 0  -n 2000000
#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016F_MINI_Correction_Iso.root  DoubleMuonLowMass_2016F_Mini_Iso -s 0  -n 2000000
#./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016G_MINI_Correction_Iso.root  DoubleMuonLowMass_2016G_Mini_Iso -s 0  -n 2000000
./sel data DoubleMuonLowMass cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/DoubleMuonLowMass_2016H_MINI_Correction_Iso.root  DoubleMuonLowMass_2016H_Mini_Iso -s 0  -n 2000000

#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016B_MINI_Correction_Iso.root Charmonium_2016B_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016C_MINI_Correction_Iso.root Charmonium_2016C_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016D_MINI_Correction_Iso.root Charmonium_2016D_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016E_MINI_Correction_Iso.root Charmonium_2016E_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016F_MINI_Correction_Iso.root Charmonium_2016F_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016G_MINI_Correction_Iso.root Charmonium_2016G_Mini_Iso -s 0  -n 100000
#./sel data Charmonium cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/crab_processed/Charmonium_2016H_MINI_Correction_Iso.root Charmonium_2016H_Mini_Iso -s 0  -n 100000


# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016B_MINI.root  DoubleMuonLowMass_2016B_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016C_MINI.root  DoubleMuonLowMass_2016C_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016D_MINI.root  DoubleMuonLowMass_2016D_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016E_MINI.root  DoubleMuonLowMass_2016E_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016F_MINI.root  DoubleMuonLowMass_2016F_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016G_MINI.root  DoubleMuonLowMass_2016G_Mini_Presel -s 0  -n 2000000
# ./sel data DoubleMuonLowMass cut0 /uscms/home/ckar/nobackup/0000/DoubleMuonLowMass_2016H_MINI.root  DoubleMuonLowMass_2016H_Mini_Presel -s 0  -n 2000000

# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016B_MINI.root Charmonium_2016B_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016C_MINI.root Charmonium_2016C_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016D_MINI.root Charmonium_2016D_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016E_MINI.root Charmonium_2016E_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016F_MINI.root Charmonium_2016F_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016G_MINI.root Charmonium_2016G_Mini_Presel -s 0  -n 2000000
# ./sel data Charmonium cut0 /uscms/home/ckar/nobackup/0000/Charmonium_2016H_MINI.root Charmonium_2016H_Mini_Presel -s 0  -n 2000000

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#done
