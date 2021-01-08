#!/bin/bash                                                                                                                                                                       
declare -i nfiles
nfiles=(86750317/2000000)+1
echo "nfiles = $nfiles"


echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
#root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToMuMuPhi-SignalMC_2016MINI_v2.rootBsToMuMuPhi-SignalMC_2016MINI_v1.root
echo -e "\n=> processing singlecand ntuple for Bs2phimm sample."

for ((num=1;num<=$nfiles;num++))
do

    ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToMuMuPhi-SignalMC_2016MINI_v1_Correction_Iso.root BsToPhiMuMu_OfficialMC_signal_2016Mini_Presel_v1 -s $((2000000*$((num-1))))  -n 2000000
    #./sel mc.lite cut_bdt /uscms/home/ckar/nobackup/BsToPhiMuMu_Signal_2016MINI_MC_v2.root BsToPhiMuMu_OfficialMC_signal_2016Mini_Presel -s $((2000000*$((num-1))))  -n 2000000
    echo -e "\n==> file$num is done."

    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done

declare -i nfiles2
nfiles2=(86805846/2000000)+1
echo "nfiles v2 = $nfiles2"

echo -e "\n=> processing singlecand ntuple for Bs2phimm sample."

for ((num=1;num<=$nfiles2;num++))
do
    ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToMuMuPhi-SignalMC_2016MINI_v2_Correction_Iso.root BsToPhiMuMu_OfficialMC_signal_2016Mini_Presel_v2 -s $((2000000*$((num-1))))  -n 2000000
    echo -e "\n==> file$num is done."
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done

#declare -i nfiles2
#files2=(128147817/2000000)+1
#echo "nfiles2 =$nfiles2"

#echo -e "\n=> processing singlecand ntuple for kstarmumu sample."

#for ((num=1;num<=$nfiles2;num++))
#do
#    ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/  KstarMuMu_OfficialMC_signal_2016Mini_Presel -s $((2000000*$((num-1))))  -n 2000000
#    echo -e "\n==> file$num is done."
#    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#done
