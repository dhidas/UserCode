#!/bin/bash

export WORKINGDIR=/cms/data15/dhidas/SUSYTrilepton/CMSSW_3_1_4/src/

echo 'WorkingDir: ' $WORKINGDIR
echo 'Process: ' $1
echo 'Name:    ' $2
echo 'OutDir:  ' $3

let FILENUMBER=${1}+1
export INFILE=`head -$FILENUMBER $2 | tail -1`
export OUTDIR=$3/$2
mkdir -vp $OUTDIR
export OUTFILE="$3/$2/$2_$1.root"
echo $INFILE $OUTFILE

rm -f ${2}_${1}.py


export VO_CMS_SW_DIR="/cms/base/cmssoft"
export COIN_FULL_INDIRECT_RENDERING=1
source $VO_CMS_SW_DIR/cmsset_default.sh

cd $WORKINGDIR
eval `scramv1 runtime -sh`

export TEMPLATE="$WORKINGDIR/DHidasAna/FillDtuple/python/FillDtuple_Template.py"

cd $OUTDIR

export MODINFILE=`echo $INFILE | sed s%\/%\\\\\\\\\/%g`
echo 'MODINFILE: ' $MODINFILE
export MODOUTFILE=`echo $OUTFILE | sed s%\/%\\\\\\\\\/%g`
echo 'MODOUTFILE: ' $MODOUTFILE

echo "cat $TEMPLATE | sed s/INFILE/$MODINFILE/g | sed s/Dtuple.root/$MODOUTFILE/g"
cat $TEMPLATE | sed s/INFILE/$MODINFILE/g | sed "s/Dtuple.root/$MODOUTFILE/g"  >> ${2}_${1}.py
echo "cmsRun $2_$1.py >& $2_$1.log"
cmsRun $2_$1.py

