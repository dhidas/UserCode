#!/bin/sh

echo 'Hello World!'
echo 'Section: ' $1
echo 'OUTDIR:  ' $2

export PATH=/cms/data15/dhidas/CernRoot/root_v5.28.00/root/bin:${PATH}
export LD_LIBRARY_PATH=/cms/data15/dhidas/CernRoot/root_v5.28.00/root/lib
. /cms/data15/dhidas/CernRoot/root_v5.28.00/root/bin/thisroot.sh


EXE=/users/h2/dhidas/Expo/RunMultiJetsRooStats

mkdir -pv $2
cd $2

echo "$EXE -1 $1"
$EXE -1 $1

