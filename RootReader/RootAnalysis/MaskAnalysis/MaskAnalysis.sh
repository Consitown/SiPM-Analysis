#!/bin/bash

g++ MaskAnalysis.cpp `root-config --libs --cflags` -o maskfile

./maskfile
