#!/bin/sh

echo 'Hello World!'
echo 'Section:    ' $1
echo 'SeedOffset: ' $2


source /osg/apps/osg/setup.sh

export PATH=/home/cdfcaf/condor/dist/bin:${PATH}

export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh


export ROOTSYS=/cms/data24/dhidas/CernRoot/root_2011.09.17_HEAD/root
export LD_LIBRARY_PATH=$ROOTSYS/lib
PATH=$PATH:$ROOTSYS/bin




EXE=$PWD/RunMultJetsCLsSplit
DAT=$PWD/Data6Jet_4632pb-1_hist.root

echo "HOSTNAME " $HOSTNAME

echo "$EXE $DAT $1 $2"
$EXE $DAT $1 $2


