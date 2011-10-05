#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_5_0_0_pre1/src
eval `scramv1 runtime -sh`
cd -

OutDir=$1
InFileName=dataGE6jets_2177pb_10GeVBins_Hist.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_lpcfarm_cls.sh $OutDir
cp -v RunMultJetsCLs $OutDir
cp -v bin/RunMultJetsCLs.cc $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_lpcfarm_cls.jdl << +EOF
OutDir = ./
Executable = run_lpcfarm_cls.sh

universe = vanilla
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $InFileName, RunMultJetsCLs
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = dhidas@FNAL.GOV
Arguments = \$(Process) \$(OutDir)
Queue 63
+EOF

# Add this for tests
#  +LENGTH="SHORT"

cd $OutDir
condor_submit submit_lpcfarm_cls.jdl
cd $ThisDir
