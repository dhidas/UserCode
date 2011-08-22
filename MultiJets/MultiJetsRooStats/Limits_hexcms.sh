#!/bin/sh

echo 'Hello World!'
DIR=$1

source /osg/apps/osg/setup.sh
export PATH=/home/cdfcaf/condor/dist/bin:${PATH}
cd /users/h2/dhidas/CMSSW_4_2_5/src
eval `scramv1 runtime -sh`
cd -



./CompileLimits $DIR/Limits_Data.dat $DIR/PELimits_*.dat
./PlotLimits Limits.dat
cp -v Limits.dat $DIR
cp -v Limits95.dat $DIR
cp -v GraphLimits.eps $DIR
cp -v LimitsPlot_*.eps NOW
mv -v LimitsPlot_*.eps $DIR
mv -v LimitsPlotN95_*.eps $DIR

