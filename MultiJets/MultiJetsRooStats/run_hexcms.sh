#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2


source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh


cd /users/h2/dhidas/CMSSW_4_2_5/src
eval `scramv1 runtime -sh`
cd -


EXE=/users/h2/dhidas/UserCode/dhidas/MultiJets/MultiJetsRooStats/RunRooStats3Param

mkdir -pv $2/Log
cd $2

echo "$EXE -1 $1"
$EXE -1 $1

