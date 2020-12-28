#!/bin/bash                                                                                                                                                                       


declare -i nfiles
nfiles=(93599410/2000000)+1
echo "nfiles = $nfiles"


echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
echo -e "\n=> processing singlecand ntuple for Bs2jpsiphi sample."


for ((num=1;num<=$nfiles;num++))
do
 #echo -e "\n./sel mc.lite mc cut_bdt /eos/uscms/store/user/ckar/BsToMuMuPhi-JpsiMC_2016MINI_Correction_Iso.root -s $((2000000*$((num-1))))  -n 2000000"
 ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToMuMuPhi-JpsiMC_2016MINI_Correction_Iso.root BsToJpsiPhi_2016MC_Official  -s $((2000000*$((num-1))))  -n 2000000

 echo -e "\n==> file$num is done."
 echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done

#declare -i nfiles2
#nfiles2=(67440041/2000000)+1
#echo "nfiles2 =$nfiles2"

#echo -e "\n=> processing singlecand ntuple for kstarmumu sample."

#for ((num=1;num<=$nfiles2;num++))
#do
 #   ./sel mc.lite mc cut_bdt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToPhiMuMu_2017_JpsikstarMC_Mini_Iso.root  JpsiKstarMuMu_OfficialMC_signal_2017Mini_Presel -s $((2000000*$((num-1))))  -n 2000000
 #   echo -e "\n==> file$num is done."
 #   echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#done

