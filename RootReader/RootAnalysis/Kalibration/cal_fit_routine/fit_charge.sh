#!/bin/bash

# compile c++ fit routine
g++ fit_charge.C `root-config --libs --cflags` -o fit_charge

#################
### CALIB FIT ###
#################

do_calib()
{	
	channel=$0
	
	f=65_8SiPM_0_Calib_lOn_Int550_HV60_TestbeamConfig_210319
	# f=77_40SiPM_0_sw2_Calib_lOn_Int590_HV60_TestbeamConfig_210319	

	./fit_charge $f $channel
	
}
export -f do_calib

channel_list=(1 2 3 4 5 6 7 8)
# channel_list=(2 3 4 5 6 7 8)
# channel_list=(1)

printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 8 bash -c "do_calib"