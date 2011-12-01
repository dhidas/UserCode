#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'List:    ' $2



source /condor/current/setup/condor_setup.sh


export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh



cd ~/CMSSW_4_2_4/src
cmsenv
cd -



OATexpert $1 $2


