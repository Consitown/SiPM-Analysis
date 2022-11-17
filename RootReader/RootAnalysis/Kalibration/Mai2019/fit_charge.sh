#!/bin/bash
rm -r ./calib_histograms/charge/

# compile c++ fit routine
g++ fit_charge.C $(root-config --libs --cflags) -o fit_charge

#################
### CALIB FIT ###
#################

do_calib() {
	channel=$0

	f=calib_vb$1

	# f=77_40SiPM_0_sw2_Calib_lOn_Int590_HV60_TestbeamConfig_210319

	./fit_charge $f $channel

}
export -f do_calib

channel_list=(0 1 2 3 4 5 6 7)

# channel_list=(2 3 4 5 6 7 8)
# channel_list=(1)

for file in rootfiles/*; do
	#without extension
	string=$file
	prefix="rootfiles/calib_vb"
	string=${string#$prefix} #Remove prefix
	suffix=".root"
	string=${string%$suffix} #Remove suffix
	printf "DOING:  %s \n" "$string"
	printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 1 bash -c "do_calib $string"
done

#for variable in 60
#  do
##  printf "TEST $variable"
#  printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 1 bash -c "do_calib $variable"

#done
