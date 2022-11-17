#!/bin/bash


cd "${0%/*}"
cwd=`pwd`


#################
### COMPILE   ###
#################

g++ BaselineAverage.C -L/opt/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/opt/root/lib -pthread -lm -ldl -rdynamic -pthread -std=c++14 -m64 -I/opt/root/include -o BL_average_new

##################
### INITIALIZE ###
##################



hist_dir="BL_histograms"

if [ ! -d "$cwd/"$hist_dir"" ]; then
		mkdir $cwd/"$hist_dir"
fi



#################
### CALIB FIT ###
#################


for file in "../rootfiles/"/*; do
	#without extension

	rootFilePath=$file
	
	prefix="../rootfiles//"
	runName=${file#$prefix} #Remove prefix
	suffix=".root"
	runName=${runName%$suffix} #Remove suffix
	
	runNr=`echo "$runName" | sed -r 's/^([^.]+).*$/\1/; s/^[^0-9]*([0-9]+).*$/\1/'`
	xLimit=1.3

	
	string="$runName $runNr  $xLimit $rootFilePath"
	#printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 1 bash -c "do_calib $string"
	#"9_calib_vb58_tune8530_pcbb 9 6075PE 40-SiPM sw5 3 83 6050PE 8-ArrayTB17 none 2."

	./BL_average_new $runName $runNr  $xLimit $rootFilePath
done

