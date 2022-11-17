#!/bin/bash

#################
### COMPILE   ###
#################

#Run with number for single file: ./FitCharge.sh 12 for run 12, without number of all 



g++ FitCharge.C -lstdc++fs $(root-config --libs --cflags) -o fit_charge_new

##################
### INITIALIZE ###
##################

cwd=$(pwd)

hist_dir="calib_histograms"
method_dir="charge"

if [ ! -d "$cwd/"$hist_dir"" ]; then
	mkdir $cwd/"$hist_dir"
fi

if [ ! -d "$cwd/"$hist_dir"/"$method_dir"" ]; then
	mkdir $cwd/"$hist_dir"/"$method_dir"
fi

#################
### CALIB FIT ###
#################

do_calib() {
	argument=$0

	./fit_charge_new $argument

}

export -f do_calib

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

	runNr=$(echo "$runName" | sed -r 's/^([^.]+).*$/\1/; s/^[^0-9]*([0-9]+).*$/\1/')
	numberOfPeaks=9

	#echo $runName $runNr "6075PE" "40-SiPM" "sw5" $numberOfPeaks $rootFilePath
	#echo $string
	#printf "%s\n" ${channel_list[@]} | xargs -n 1 -P 1 bash -c "do_calib $string"

	#"9_calib_vb58_tune8530_pcbb 9 6075PE 40-SiPM sw5 10"

	./fit_charge_new $runName $runNr "6075PE" "40-SiPM" "sw5" $numberOfPeaks $rootFilePath
done

data_list=(
	# DARK COUNT

	# f=74_40SiPM_0_DarkCount_lOff_HV60_TestbeamConfig_210319
	# f=78_40SiPM_0_sw2_DarkCount_lOff_HV60_TestbeamConfig_210319
	# f=80_40SiPM_0_sw1_DarkCount_lOff_HV60_TestbeamConfig_210319
	# f=83_8SiPM_17_DarkCount_lOff_HV60_TestbeamConfig_210319
	# f=86_40SiPM_1_DarkCount_lOff_HV60_TestbeamConfig_210319
	# f=86_40SiPM_2_sw4_DarkCount_lOff_HV60_TestbeamConfig_210319

	#CALIB RUNS
	# scheme:
	# filename run_nr sipm_id wom_id sw_id n_peaks

	#"\na_I3_HV6000_p331_allChannels_allSum"
	#"\nb_I3_HV6000_p331_allChannels_allSum"
	#"\nc_I3_HV6000_p331_allChannels_allSum"
	#"\nd_I3_HV6000_p331_allChannels_allSum"
	#"\ne_I3_HV6000_p331_allChannels_allSum"
	#"\nf_I3_HV6000_p331_allChannels_allSum"
	#"\ng_I3_HV6000_p331_allChannels_allSum"
	#"\nh_I3_HV6000_p331_allChannels_allSum"

	# "\n38_40SiPM_1_sw4_Calib_lOff_Int554_HV60_TestbeamConfig_240119 38 3075PE WOM-C sw4"

	# "\n44_40SiPM_2_sw4_Calib_lOff_Int554_HV60_TestbeamConfig_240119 44 3075PE WOM-D sw4"

	#"\n51_8SiPM_1_Calib_lOff_Int536_HV60_TestbeamConfig_040219 51 6075PE WOM-A none"
	#"\n52_8SiPM_2_Calib_lOff_Int536_HV60_TestbeamConfig_040219 52 6075PE WOM-B none"
	#"\n64_8SiPM_0_Calib_lOn_Int544_HV60_TestbeamConfig_210319 64 6075PE 8-Array0 none"

	#"\n72_40SiPM_0_sw4_Calib_lOn_Int590_HV60_TestbeamConfig_210319 72 3075PE 40-Array0 sw4"

	#"\n75_40SiPM_0_sw4_Calib_lOn_Int570_HV60_TestbeamConfig_210319 75 3075PE 40-Array0 sw4"
	#"\n76_40SiPM_0_sw2_Calib_lOn_Int570_HV60_TestbeamConfig_210319 76 3075PE 40-Array0 sw2"

	#"\n79_40SiPM_0_sw1_Calib_lOn_Int590_HV60_TestbeamConfig_210319 79 3075PE 40-Array0 sw1"

	#"\n85_40SiPM_1_sw4_Calib_lOn_Int570_HV60_TestbeamConfig_210319 85 3075PE WOM-C sw4"
	#"\n87_40SiPM_2_sw4_Calib_lOn_Int570_HV60_TestbeamConfig_210319 87 3075PE WOM-D sw4"

	# "\n53_8SiPM_17_Calib_lOff_Int530_HV60_TestbeamConfig_040219 53 6050PE 8-ArrayTB17 none"

	## calib data for final analysis, Sep. 2019

	# // WOM-C/-D,40Array0 sw1
	# "41_40SiPM_1_sw2_Calib_lOff_Int590_HV60_TestbeamConfig_240119 41 3075PE WOM-C sw2 6"
	# "\n47_40SiPM_2_sw2_Calib_lOff_Int574_HV60_TestbeamConfig_040219 47 3075PE WOM-D sw2 5"
	# "\n81_40SiPM_0_sw1_Calib_lOn_Int600_HV60_TestbeamConfig_210319 81 3075PE 40-Array0 sw1 6"

	# # // WOM-C/-D,40Array0 sw2
	# "\n42_40SiPM_1_sw1_Calib_lOff_Int600_HV60_TestbeamConfig_240119 42 3075PE WOM-C sw1 5"
	# "\n50_40SiPM_2_sw1_Calib_lOff_Int600_HV60_TestbeamConfig_040219 50 3075PE WOM-D sw1 5"
	# "\n77_40SiPM_0_sw2_Calib_lOn_Int590_HV60_TestbeamConfig_210319 77 3075PE 40-Array0 sw2 6"

	# # // WOM-C/-D,40Array0 sw4
	# "\n84_40SiPM_1_sw4_Calib_lOn_Int590_HV60_TestbeamConfig_210319 84 3075PE WOM-C sw4 7"
	# "\n88_40SiPM_2_sw4_Calib_lOn_Int580_HV60_TestbeamConfig_210319 88 3075PE WOM-D sw4 5"
	# "\n73_40SiPM_0_sw4_Calib_lOn_Int590_HV60_TestbeamConfig_210319 73 3075PE 40-Array0 sw4 7"

	# # // WOM-A/-B,8Array0
	# "\n65_8SiPM_0_Calib_lOn_Int550_HV60_TestbeamConfig_210319 65 6075PE 8-Array0 none 6"
	# "\n68_8SiPM_1_Calib_lOn_Int550_HV60_TestbeamConfig_210319 68 6075PE WOM-A none 6"
	# "\n69_8SiPM_2_Calib_lOn_Int550_HV60_TestbeamConfig_210319 69 6075PE WOM-B none 6"

	# // 8ArrayTB17
	# "\n82_8SiPM_17_Calib_lOn_Int550_HV60_TestbeamConfig_210319 82 6050PE 8-ArrayTB17 none 5"
	"9_calib_vb58_tune8530_pcbb 9 6075PE 40-SiPM sw5 10"

)

#echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "do_calib"

echo

#################
### DEBUGGING ###
#################

# Test functions for fast debugging with xargs
test_fn() {
	echo -n "-> "
	for a in "$0"; do echo -n "\"$a\" "; done
	echo
	# sleep 1 # show xargs parallel mode
}
export -f test_fn

test_fn2() {
	echo $0
	echo $1
	echo $2
	echo $3
}
export -f test_fn2

# echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "test_fn"
