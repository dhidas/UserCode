#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -

OutDir=$1
InFileName=data6jetpt70_160offset_1089pb3fits10GeVBins.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_lpcfarm.sh $OutDir
cp -v RunRooStats3Param $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_lpcfarm.jdl << +EOF
OutDir = ./
Executable = run_lpcfarm.sh

universe = vanilla
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $InFileName, RunRooStats3Param
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = dhidas@FNAL.GOV
Arguments = \$(Process)
Queue 1000
+EOF


cd $OutDir
condor_submit submit_lpcfarm.jdl
cd $ThisDir
