#!/bin/sh

echo 'Hello World!'

SECTION=$1
OUTDIR=$2
INPUTLIST=$3
RELEASEDIR=$4

PYFILE=$OUTDIR/Dtuple_Template_cfg.py

echo 'Section:     ' $SECTION
echo 'OUTDIR:      ' $OUTDIR
echo 'INPUTLIST:   ' $INPUTLIST
echo 'RELEASEDIR:  ' $RELEASEDIR
echo ''
echo 'PYFILE:      ' $PYFILE

# Setup release
cd $RELEASEDIR
export VO_CMS_SW_DIR=/cms/base/cmssoft
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh
eval `scramv1 ru -sh`


#Go to the out dir
cd $OUTDIR

# get the infile
let LINENUMBER=$SECTION+1;
echo 'LINENUMBER: ' $LINENUMBER
INPUTFILE=file:`head -$LINENUMBER $INPUTLIST | tail -1`
echo "File is: $INPUTFILE"

OUTFILENAME="Dtuple_`basename $INPUTFILE`"

# run the cmsRun job
echo cmsRun $PYFILE print inputFiles=$INPUTFILE outputFile=$OUTDIR/$OUTFILENAME
cmsRun $PYFILE print inputFiles=$INPUTFILE outputFile=$OUTDIR/$OUTFILENAME





