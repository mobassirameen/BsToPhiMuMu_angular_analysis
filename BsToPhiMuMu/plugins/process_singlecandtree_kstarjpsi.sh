#!/bin/bash                                                                                                                                                                       
declare -i nfiles
nfiles=(72244124/2000000)+1
echo "nfiles = $nfiles"


echo -e "\n@@@@  shell script is going to produce $nfiles files  @@@@"
echo -e "\n=> processing singlecand ntuple for BdToJpsikstar sample." 
for ((num=1;num<=$nfiles;num++))
do

    ./sel mc.lite mc cutopt root://cmseos.fnal.gov//eos/uscms/store/user/ckar/BsToMuMuPhi-JpsiKstar_2016MINI_Correction_Iso_v2.root BsToMuMuPhi-JpsiKstar_2016MINI_Correction_Iso_v2 -s $((2000000*$((num-1))))  -n 2000000

    echo -e "\n==> file$num is done."

    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
done

