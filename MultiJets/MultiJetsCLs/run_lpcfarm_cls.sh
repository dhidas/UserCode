echo 'Hello World!'
echo 'Section: ' $1



source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_5_0_0_pre1/src
eval `scramv1 runtime -sh`
cd -



EXE=$_CONDOR_SCRATCH_DIR/RunMultJetsCLs
DAT=$_CONDOR_SCRATCH_DIR/dataGE6jets_2177pb_10GeVBins_Hist.root


echo "$EXE $DAT $1"
$EXE $DAT $1

