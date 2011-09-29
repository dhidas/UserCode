#!/bin/sh

echo 'Hello World!'
DIR=$1

source /uscmst1/prod/sw/cms/shrc prod
<<<<<<< Limits_lpc.sh
cd /uscms/home/dhidas/CMSSW_4_2_3/src
=======
#source /uscmst1/prod/grid/gLite_SL5.sh

cd /uscms/home/dhidas/CMSSW_4_2_3/src
>>>>>>> 1.2
eval `scramv1 runtime -sh`
cd -


<<<<<<< Limits_lpc.sh

OUTDIR=$DIR/LimitPlots
mkdir -pv $OUTDIR
=======
>>>>>>> 1.2

./CompileLimits $DIR/Limits_Data.dat $DIR/PELimits_*.dat
./PlotLimits Limits.dat
cp -v Limits.dat $OUTDIR
cp -v Limits95.dat $OUTDIR
cp -v GraphLimits.eps $OUTDIR
mv -v LimitsPlot_*.eps $OUTDIR
mv -v LimitsPlotN95_*.eps $OUTDIR

