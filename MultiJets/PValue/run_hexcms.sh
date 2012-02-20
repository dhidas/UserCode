#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1

source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh


cd /users/h2/dhidas/CMSSW_4_2_5/src
eval `scramv1 runtime -sh`
cd -


EXE=./RunPValue

DAT=dataGE6jets_2177pb_10GeVBins_Hist.root


echo "$EXE $DAT $1"
$EXE $DAT $1

