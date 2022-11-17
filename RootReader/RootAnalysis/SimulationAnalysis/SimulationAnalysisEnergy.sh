#!/bin/bash

#################
### COMPILE   ###
#################

g++ SimulationAnalysisEnergy.C -lstdc++fs `root-config --libs --cflags` -o SimulationAnalysisEnergy

##################
### INITIALIZE ###
#################



 ./SimulationAnalysisEnergy

 
done

 


