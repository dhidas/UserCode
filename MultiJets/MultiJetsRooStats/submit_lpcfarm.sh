#!/bin/sh

echo 'Submitting to hexfarm!'

ThisDir=$PWD

source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -

OutDir=$1

mkdir -pv $OutDir/Log $OutDir

cp -v run_lpcfarm.sh $OutDir
cp -v RunRooStats3Param $OutDir
cp -v Data/DijetMassFit_data_881pb-1_6jets_pt70.root $OutDir
cp -v submit_lpcfarm.jdl $OutDir

cd $OutDir
condor_submit submit_lpcfarm.jdl
cd $ThisDir
