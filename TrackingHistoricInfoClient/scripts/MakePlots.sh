#!/bin/sh

TagName=$1
Password=$2
RunStart=$3
RunEnd=$4

PlotDir="CurrentPlots"

rm -rf $PlotDir

if [ $4 ]; then
TrackingHDQMInspector $TagName $Password $RunStart $RunEnd
else
TrackingHDQMInspector $TagName $Password $RunStart
fi
mkdir -pv $PlotDir
mv *.gif $PlotDir
mv *.eps $PlotDir
cp diow.pl $PlotDir
cd $PlotDir
./diow.pl
cd ..

