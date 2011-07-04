#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2


source /uscmst1/prod/sw/cms/shrc prod
#source /uscmst1/prod/grid/gLite_SL5.sh

cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -



EXE=/uscms/home/dhidas/UserCode/dhidas/MultiJets/MultiJetsRooStats/RunMultiJetsRooStats

mkdir -pv $2
cd $2

echo "$EXE $1"
$EXE $1

