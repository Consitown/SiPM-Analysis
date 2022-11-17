#!/bin/bash

#################
### COMPILE   ###
#################

g++ TestGround.C -lstdc++fs `root-config --libs --cflags` -o testGround

##################
### INITIALIZE ###
#################


 ./testGround

 


