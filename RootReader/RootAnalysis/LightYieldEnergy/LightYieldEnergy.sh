#!/bin/bash

#################
### COMPILE   ###
#################
cd "${0%/*}"

g++ LightYieldEnergy.C -lstdc++fs -L/opt/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/opt/root/lib -pthread -lm -ldl -rdynamic -pthread -std=c++14 -m64 -I/opt/root/include -o lightyieldenergy -lSpectrum

##################
### INITIALIZE ###
#################

paths14=()
paths26=()
paths52=()
positions=()

for file in "../rootfiles/"/*; do
	#without extension

	if [ ! -z "$1" ]; then

		if [[ $file == *"$1"* ]]; then
			echo $line
		else
			continue
		fi

	fi

	rootFilePath=$file

	prefix="../rootfiles//"
	runName=${file#$prefix} #Remove prefix
	suffix=".root"
	runName=${runName%$suffix} #Remove suffix

	if [[ $runName == *"dc"* ]]; then
		continue
	fi

	runNr=$(echo "$runName" | sed -r 's/^([^.]+).*$/\1/; s/^[^0-9]*([0-9]+).*$/\1/')

	if [[ $runName == *"14"* ]]; then
		paths14+=($runName)
	fi
	if [[ $runName == *"26"* ]]; then
		paths26+=($runName)
	fi
	if [[ $runName == *"52"* ]]; then
		paths52+=($runName)
	fi
	#15_pos4_angle0_e26_ch32

	IFS='_' read -ra runInfos <<<"$runName"

	for i in "${runInfos[@]}"; do

		if [[ $i == *"pos"* ]]; then
			position=$(echo "$i" | sed 's/[^0-9]*//g')
			#echo "Position: $position"
		fi
	done

	positions+=($position)

done

echo "1.4Gev: ${paths14[@]}"
echo "2.6Gev: ${paths26[@]}"
echo "5.2Gev: ${paths52[@]}"

positionsU=($(printf "%s\n" "${positions[@]}" | sort -u | tr '\n' ' '))

echo "RunNumbers: ${positions[@]}"
echo "RunNumbersUnique: ${positionsU[@]}"


IFS=\, eval 'str14="${paths14[*]}"'
IFS=\, eval 'str26="${paths26[*]}"'
IFS=\, eval 'str52="${paths52[*]}"'
IFS=\, eval 'strR="${positionsU[*]}"'

./lightyieldenergy ${str14} ${str26} ${str52} ${strR}
