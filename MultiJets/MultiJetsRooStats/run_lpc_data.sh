#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2



source /uscmst1/prod/sw/cms/cshrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -csh`
cd -



EXE=/uscms/home/dhidas/UserCode/dhidas/MultiJets/MultiJetsRooStats/RunMultiJetsRooStats

mkdir -pv $2
cd $2

echo "$EXE -1 $1"
$EXE -1 $1

