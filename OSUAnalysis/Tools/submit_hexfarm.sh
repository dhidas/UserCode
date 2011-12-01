#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh



source /condor/current/setup/condor_setup.sh



OutDir=$1
ListName=$2

mkdir -pv $OutDir/Log $OutDir
cp $ListName $OutDir
cp run_hexfarm.sh $OutDir
cp pileup_160404-165970.root  pileup_160404-166502.root $OutDir

FileName=`basename $ListName`

cat > $OutDir/submit_hexfarm.jcl << +EOF
OutDir = ./
Executable = run_hexfarm.sh

universe = vanilla
Should_Transfer_Files = NO
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = dhidas@physics.rutgers.edu
Arguments = \$(Process) $FileName
Queue 20
+EOF



cd $OutDir
condor_submit submit_hexfarm.jcl
cd $ThisDir
