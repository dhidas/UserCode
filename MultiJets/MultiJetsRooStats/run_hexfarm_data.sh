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


EXE=$PWD/RunRooStats3Param
DAT=$PWD/data6jetpt70_160offset_1089pb3fits10GeVBins.root

echo "HOSTNAME " $HOSTNAME

echo "$EXE $DAT -1 $1"
$EXE $DAT -1 $1


