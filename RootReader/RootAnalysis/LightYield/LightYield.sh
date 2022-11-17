#!/bin/bash

#################
### COMPILE   ###
#################
echo "Compiling LightYield.C"
g++ LightYield.C -L/opt/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/opt/root/lib -pthread -lm -ldl -rdynamic -pthread -std=c++14 -m64 -I/opt/root/include -o pe_dist


##################
### INITIALIZE ###
##################

cwd=`pwd`
hist_dir="ly_histograms"
method_dir="charge"

if [ ! -d "$cwd/"$hist_dir"" ]; then
		mkdir $cwd/"$hist_dir"
fi

if [ ! -d "$cwd/"$hist_dir"/"$method_dir"" ]; then
		mkdir $cwd/"$hist_dir"/"$method_dir"
fi

if [ ! -d "$cwd/"$hist_dir"/"$method_dir"/CH32" ]; then
		mkdir $cwd/"$hist_dir"/"$method_dir"/CH32
fi

if [ ! -d "$cwd/"$hist_dir"/"$method_dir"/CH16" ]; then
		mkdir $cwd/"$hist_dir"/"$method_dir"/CH16
fi

peak_gFit()
{
	./pe_dist $1
}

for entry in "../rootfiles/"/*
do
  echo "$entry"
#../rootfiles//26_pos7_angle0_e52_ch32.root
 if [ -z "$entry" ]; then
    echo "VAR is empty"
	else 
 peak_gFit "$entry"
fi
 
done