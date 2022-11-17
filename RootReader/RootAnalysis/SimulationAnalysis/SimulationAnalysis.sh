#!/bin/bash

#################
### COMPILE   ###
#################

g++ SimulationAnalysis.C -lstdc++fs `root-config --libs --cflags` -o SimulationAnalysis

##################
### INITIALIZE ###
#################


for entry in "../rootfiles/"/*
do
  echo "$entry"
#../rootfiles//26_pos7_angle0_e52_ch32.root
 if [ -z "$entry" ]; then
    echo "VAR is empty"
	else 
 ./SimulationAnalysis "$entry"
fi
 
done

 


