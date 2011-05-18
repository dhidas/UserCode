#!/bin/sh

echo 'Hello World!'

SECTION=$1
OUTDIR=$2
INPUTLIST=$3
NSECTIONS=$4
RELEASEDIR=$5

PYFILE=$OUTDIR/Skim_Template_cfg.py

echo 'Section:     ' $SECTION
echo 'OUTDIR:      ' $OUTDIR
echo 'INPUTLIST:   ' $INPUTLIST
echo 'NSECTIONS:   ' $NSECTIONS
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

LINES=`wc -l $INPUTLIST  | awk '{print \$1}'`
let NPERSEC=$LINES/$NSECTIONS;
let REMAINDER=$LINES%$NSECTIONS;
echo 'NPERSEC:     ' $NPERSEC
echo 'REMAINDER:   ' $REMAINDER

# if section < remainder - 1 add another one =)

if [[ $SECTION -lt $REMAINDER ]]
then let BEGIN=$SECTION*$NPERSEC+$SECTION
else let BEGIN=$SECTION*$NPERSEC+$REMAINDER
fi

if [[ $SECTION -lt $REMAINDER ]]
then let NUM=$NPERSEC+1
else let NUM=$NPERSEC
fi



echo 'BEGIN:   ' $BEGIN
echo 'NUM:     ' $NUM


# get the infile
let LINENUMBER=$BEGIN+$NUM;
echo 'LINENUMBER: ' $LINENUMBER
echo "head -$LINENUMBER $INPUTLIST | tail -$NUM"
echo `head -$LINENUMBER $INPUTLIST | tail -$NUM`
declare -a INPUTFILE
INPUTFILE=`head -$LINENUMBER $INPUTLIST | tail -$NUM  | sed 's/^/file:/g' | tr '[:space:]' ','`
INPUTFILE=${INPUTFILE%?}
#echo "File is: $INPUTFILE"
OUTFILENAME="Skim_$SECTION.root"

# run the cmsRun job
echo cmsRun $PYFILE print inputFiles=$INPUTFILE outputFile=$OUTDIR/$OUTFILENAME
cmsRun $PYFILE print inputFiles=$INPUTFILE outputFile=$OUTDIR/$OUTFILENAME





