#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2


source /uscmst1/prod/sw/cms/shrc prod
source /uscmst1/prod/grid/gLite_SL5.sh

cd ~dhidas/CMSSW_3_10_0/src
cmsenv
cd -


export CMS_PATH=/uscmst1/prod/sw/cms
export PYTHONDIR=/uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8

export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.27.06/bin:$PATH
export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.27.06
export PYTHONPATH=${ROOTSYS}/lib

export LD_LIBRARY_PATH=${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export ROOT_INCLUDE=${ROOTSYS}/include

export PATH=${PWD}:${PATH}




INPUTFILE=/uscms/home/dhidas/Data35pb/LandauFit_data_35pb-1_6jets_and_scaled_4jets_pt45.root
EXE=/uscms/home/dhidas/UserCode/dhidas/MultiJets/MultiJetsRooStats/RunMultiJetsRooStats

mkdir -pv $2
cd $2

echo "$EXE $1"
$EXE $1

