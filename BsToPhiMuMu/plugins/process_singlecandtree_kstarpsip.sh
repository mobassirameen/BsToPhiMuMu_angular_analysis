#!/bin/bash                                                                                                                                                                       
declare -i nfiles
nfiles=(7471553/2000000)+1
echo "nfiles = $nfiles"


echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
echo -e "\n=> processing singlecand ntuple for BdToKstarmumu sample."

for ((num=4;num<=$nfiles;num++))
do

    ./sel mc.lite mc cutopt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToPhiMuMu_2018_psikstarMC_Mini_Correction_Iso_v2.root BsToPhiMuMu_2018_psikstarMC_Mini_Correction_Iso_v2 -s $((2000000*$((num-1))))  -n 2000000
    #./sel mc.lite cut_bdt /uscms/home/ckar/nobackup/BsToPhiMuMu_Signal_2016MINI_MC_v2.root BsToPhiMuMu_OfficialMC_signal_2016Mini_Presel -s $((2000000*$((num-1))))  -n 2000000
    echo -e "\n==> file$num is done."

    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done

