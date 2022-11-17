#!/bin/bash

#################
### COMPILE   ###
#################

g++ SimulationAnalysisPos.C -lstdc++fs `root-config --libs --cflags` -o SimulationAnalysisPos

##################
### INITIALIZE ###
#################



 ./SimulationAnalysisPos

 


