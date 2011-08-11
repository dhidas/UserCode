#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh

source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}

cd /users/h2/dhidas/CMSSW_4_2_5/src
eval `scramv1 runtime -sh`
cd -



OutDir=$1
InFileName=data6jetpt70_160offset_1089pb3fits10GeVBins.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_hexfarm_data.sh $OutDir
cp -v RunRooStats3Param $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_hexfarm_data.jdl << +EOF
OutDir = ./
Executable = run_hexfarm_data.sh

universe = vanilla
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = NO
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = dhidas@FNAL.GOV
Arguments = \$(Process)
Queue 26
+EOF



cd $OutDir
condor_submit submit_hexfarm_data.jdl
cd $ThisDir
