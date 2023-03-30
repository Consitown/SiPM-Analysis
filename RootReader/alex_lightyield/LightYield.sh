#!/bin/bash

#################
### COMPILE   ###
#################
echo "Compiling LightYield.C"
g++ LightYield.C `root-config --libs --cflags` -o pe_dist


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

export -f peak_gFit

data_list=(

	# runName runNr pos WC_version particle energy target_wom

	# ## Ch32
	# ___ POS 0 ––––––
#	"1_pos0_angle0_e26_ch32 1 0 CH32 electron 26 CH32"
	# # ___ POS 2 ––––––
#	 "\n2_pos0_angle0_e14_ch32 2 0 CH32 electron 14 CH32"
#	# # ___ POS 3 ––––––
#	 "\n3_pos0_angle0_e52_ch32 3 0 CH32 electron 52 CH32"
	# # ___ POS 4 ––––––
#	 "24_pos7_angle0_e14_ch32 24 7 CH32 electron 14 CH32"
#	 "\n25_pos7_angle0_e26_ch32 25 7 CH32 electron 26 CH32"
	#"\n26_pos7_angle0_e52_ch32 26 7 CH32 electron 52 CH32"

 	"\n27_pos8_angle0_e52_ch32 27 8 CH32 electron 52 CH32"
	 "\n28_pos8_angle0_e26_ch32 27 8 CH32 electron 26 CH32"
	"\n29_pos8_angle0_e14_ch32 29 8 CH32 electron 14 CH32"
	# # ___ POS 7 ––––––
	# "\n25_pos7_angle0_e26_ch32 25 7 CH32 electron 26 CH32"
	# # ___ SAT TEST w/ offstet w/ CH16  ––––––
	# "\n100_saturation_test 100 666 CH32 electron 26 CH32"
	# # ___ SAT TEST w/o offstet w/o CH16  ––––––
	# "\n101_wcsmalloff 101 666 CH32 electron 26 CH32"
	# # ___ SAT TEST w/ offstet w/o CH16  ––––––
	# "\n102_wcsmallof_offset113 102 666 CH32 electron 26 CH32"
	# # ___ POS 10 angle 30, SAT TEST w/ offstet w/o CH16  ––––––
	# "\n103_pos10_offsetON 103 666 CH32 electron 26 CH32"
	# # ___ POS 10 angle 30, SAT TEST w/o offstet w/o CH16  ––––––
	# "\n104_pos10_offsetOFF 104 666 CH32 electron 26 CH32"
	

	

	)



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



#echo ${arr[@]}




#echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "peak_gFit"
#echo -e ${arr[@]} | xargs -0 -n 1 -P 1 bash -c "peak_gFit"







#################
### DEBUGGING ###
#################

# Test functions for fast debugging with xargs
test_fn()
{
	echo -n "-> "; for a in "$0"; do echo -n "\"$a\" "; done; echo
	# sleep 1 # show xargs parallel mode 
}
export -f test_fn

test_fn2()
{
	echo $0; echo $1; echo $2; echo $3;
}
export -f test_fn2

# echo -e ${data_list[@]} | tr '\n' '\0' | xargs -0 -n 1 -P 1 bash -c "test_fn"