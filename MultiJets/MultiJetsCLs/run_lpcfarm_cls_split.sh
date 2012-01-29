echo 'Hello World!'
echo 'Section:    ' $1
echo 'SeedOffset: ' $2



source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_5_0_0_pre1/src
eval `scramv1 runtime -sh`
cd -



EXE=$_CONDOR_SCRATCH_DIR/RunMultJetsCLsSplit
DAT=$_CONDOR_SCRATCH_DIR/Data6Jet_4632pb-1_hist.root


echo "$EXE $DAT $1 $2"
$EXE $DAT $1 $2

