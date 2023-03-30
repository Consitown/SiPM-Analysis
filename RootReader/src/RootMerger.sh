#!/bin/bash
#Teilweise bekommt man Dateien mit Unicode Zeichen im Namen, liegt entweder an CRLF oder sonst woran. Einfach diese Datei als Basis nehmen und Copy Pasten, dann gehts.
#rm read
# -rpath option might be necessary with some ROOT installations 
# -lSpectrum option might be necessary with some ROOT installations 
#g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -lSpectrum -o read 
g++ mergeROOTFiles.C `root-config --libs --cflags` -lSpectrum -o mergeROOTFiles

here=`pwd` #Print Working Directory Substitution


#$here/mergeROOTFiles $here/runs/Fast/22_muon6_pos4/
