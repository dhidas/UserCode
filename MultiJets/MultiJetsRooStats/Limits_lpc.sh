#!/bin/sh

echo 'Hello World!'
DIR=$1

source /uscmst1/prod/sw/cms/shrc prod
#source /uscmst1/prod/grid/gLite_SL5.sh

cd /uscms/home/dhidas/CMSSW_4_2_3/src
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

