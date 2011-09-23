#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -

OutDir=$1
InFileName=dataGE6jets_2177pb_10GeVBins_Hist.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_lpc.sh $OutDir
cp -v bin/RunPValue.cc $OutDir
cp -v RunPValue $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_lpc.jdl << +EOF
OutDir = ./
Executable = run_lpc.sh

universe = vanilla
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $InFileName, RunPValue
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = dhidas@FNAL.GOV
Arguments = \$(Process)
Queue 1
+EOF


cd $OutDir
condor_submit submit_lpc.jdl
cd $ThisDir
