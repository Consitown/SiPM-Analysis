#!/bin/bash

#################
### COMPILE   ###
#################

g++ EnergyAnalysis.C -lstdc++fs `root-config --libs --cflags` -o energyAnalysis

##################
### INITIALIZE ###
#################


 ./energyAnalysis

 



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