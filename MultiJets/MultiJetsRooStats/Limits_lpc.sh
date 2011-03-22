#!/bin/sh

echo 'Hello World!'
DIR=$1

source /uscmst1/prod/sw/cms/shrc prod
#source /uscmst1/prod/grid/gLite_SL5.sh

cd /uscms/home/dhidas/CMSSW_3_10_0/src
eval `scramv1 runtime -sh`
cd -


export CMS_PATH=/uscmst1/prod/sw/cms
export PYTHONDIR=/uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8

export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.27.06/bin:$PATH
export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.27.06
export PYTHONPATH=${ROOTSYS}/lib

export LD_LIBRARY_PATH=${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export ROOT_INCLUDE=${ROOTSYS}/include

export PATH=${PWD}:${PATH}



./CompileLimits $DIR/Limits_Data.dat $DIR/PELimits_*.dat
./PlotLimits Limits.dat
cp -v Limits.dat $DIR
cp -v Limits95.dat $DIR
cp -v GraphLimits.eps $DIR
cp -v LimitsPlot_*.eps NOW
mv -v LimitsPlot_*.eps $DIR
mv -v LimitsPlotN95_*.eps $DIR
