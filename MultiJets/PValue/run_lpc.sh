#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2

InFileName=dataGE6jets_2177pb_10GeVBins_Hist.root

source /uscmst1/prod/sw/cms/shrc prod
#source /uscmst1/prod/grid/gLite_SL5.sh

cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -


#export CMS_PATH=/uscmst1/prod/sw/cms
#export PYTHONDIR=/uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8

#export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.27.06/bin:$PATH
#export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.27.06
#export PYTHONPATH=${ROOTSYS}/lib

export LD_LIBRARY_PATH=${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export ROOT_INCLUDE=${ROOTSYS}/include
export PATH=${PWD}:${PATH}




EXE=./RunPValue


echo "$EXE $InFileName $1"
$EXE $InFileName $1

