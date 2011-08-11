echo 'Hello World!'
echo 'Section: ' $1



source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -



EXE=$_CONDOR_SCRATCH_DIR/RunRooStats3Param
DAT=$_CONDOR_SCRATCH_DIR/data6jetpt70_160offset_1089pb3fits10GeVBins.root


echo "$EXE $DAT $1"
$EXE $DAT $1

