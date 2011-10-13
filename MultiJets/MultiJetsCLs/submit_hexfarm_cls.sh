#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh

source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}




OutDir=$1
InFileName=dataGE6jets_2177pb_10GeVBins_Hist.root

mkdir -pv $OutDir/Log $OutDir

cp -v run_hexfarm_cls.sh $OutDir
cp -v RunMultJetsCLs $OutDir
cp -v bin/RunMultJetsCLs.cc $OutDir
cp -v Data/$InFileName $OutDir


cat > $OutDir/submit_hexfarm_cls.jdl << +EOF
OutDir = ./
Executable = run_hexfarm_cls.sh

universe = vanilla
Should_Transfer_Files = NO
Output = \$(OutDir)/Log/Log_\$(Process).out
Error =  \$(OutDir)/Log/Log_\$(Process).err
Log =    \$(OutDir)/Log/Log_\$(Process).log
notify_user = clseitz@physics.rutgers.edu
Arguments = \$(Process)
Queue 63
+EOF



cd $OutDir
condor_submit submit_hexfarm_cls.jdl
cd $ThisDir
