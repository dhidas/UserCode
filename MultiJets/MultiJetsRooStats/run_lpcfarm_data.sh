echo 'Hello World!'
echo 'Section: ' $1



source /uscmst1/prod/sw/cms/shrc prod
cd /uscms/home/dhidas/CMSSW_4_2_3/src
eval `scramv1 runtime -sh`
cd -



EXE=$_CONDOR_SCRATCH_DIR/RunRooStats3Param
DAT=$_CONDOR_SCRATCH_DIR/dataGE6jets_2177pb_10GeVBins_Hist.root
#DAT=$_CONDOR_SCRATCH_DIR/DijetMassFit_data_881pb-1_6jets_pt70.root


echo "$EXE $DAT -1 $1"
$EXE $DAT -1 $1

