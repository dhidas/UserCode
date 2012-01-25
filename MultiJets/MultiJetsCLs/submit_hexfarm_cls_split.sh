#!/bin/sh

if [ $# -ne 2 ]
then
echo "Usage: submit_hexfarm_cls_split.sh OutDir SeedOffset"
exit
else
echo 'Submitting to hexfarm!'
fi


ThisDir=$PWD

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh

source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}




OutDir=$1
SeedOffset=$2

InFileName=Data6Jet_4632pb-1_hist.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_hexfarm_cls_split.sh $OutDir
cp -v RunMultJetsCLsSplit $OutDir
cp -v bin/RunMultJetsCLsSplit.cc $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_hexfarm_cls_split.jdl << +EOF
OutDir = ./
Executable = run_hexfarm_cls_split.sh

universe = vanilla
Should_Transfer_Files = NO
Output = \$(OutDir)/Log/Log_\$(Process)_$SeedOffset.out
Error =  \$(OutDir)/Log/Log_\$(Process)_$SeedOffset.err
Log =    \$(OutDir)/Log/Log_\$(Process)_$SeedOffset.log
notify_user = dhidas@physics.rutgers.edu
Arguments = \$(Process) $SeedOffset
Queue 63
+EOF



cd $OutDir
condor_submit submit_hexfarm_cls_split.jdl
cd $ThisDir
