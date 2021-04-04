#!/bin/bash                                                                                                                                                                       


declare -i nfiles
nfiles=(43141393/2000000)+1
echo "nfiles = $nfiles"


echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
echo -e "\n=> processing singlecand ntuple for Bs2psipphi sample."


for ((num=1;num<=$nfiles;num++))
do
 #echo -e "\n./sel mc.lite cut_bdt /eos/user/r/rraturi/BsToPhiMuMu/DataBase/IsoVar/2016/MC/uPhi-PsiMC_2016MINI_Iso.root sel_BsToPsiPhi_OfficialMC_signal_2016 -s $((2000000*$((num-1))))  -n 2000000"
 ./sel mc.lite mc cutopt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToPhiMuMu_2018_PsiPhiMC_Mini_Correction_Iso.root BsToPhiMuMu_2018_PsiPhiMC_Mini_Correction_Iso -s $((2000000*$((num-1))))  -n 2000000
 echo -e "\n==> file$num is done."
 echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done




# for psikstar

#declare -i nfiles2
#nfiles2=(9237624/2000000)+1
#echo "nfiles2 =$nfiles2"

#echo -e "\n=> processing singlecand ntuple for kstarmumu sample."

#for ((num=1;num<=$nfiles2;num++))
#do
#    ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToPhiMuMu_2017_psikstarMC_Mini_Iso.root PsiKstarMuMu_OfficialMC_signal_2017Mini_Presel -s $((2000000*$((num-1))))  -n 2000000
#    echo -e "\n==> file$num is done."
#    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#done

