#!/bin/bash

g++ SumAnalysis.cpp `root-config --libs --cflags` -o sumfile

./sumfile
