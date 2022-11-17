#!/bin/bash

#################
### COMPILE   ###
#################
cd "${0%/*}"

g++ EfficiencyPos.C -lstdc++fs `root-config --libs --cflags` -o efficiencypos -lSpectrum 

##################
### INITIALIZE ###
#################

./efficiencypos

 


