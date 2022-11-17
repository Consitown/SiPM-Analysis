#!/bin/bash

# Print out list of the .root files to be analyzed:
echo "Analyzing the following files: "
for file in "../data/"*;do
	echo $file
done

# Ask user if they want to proceed (choose N/n e.g. if the file list wasn't correct):
while true; do
	read -p "Proceed? (This will remove the previously compiled version of the program, so be careful!) [Y/n] " proceed
	case $proceed in
		[Yy]* ) rm PhotonCountDistribution;
				echo "Proceeding..."; break;;
		[Nn]* ) echo "Program aborted."; exit;;
		* ) echo "Please type Y/y or N/n.";;
	esac
done

# Compile the .cpp-file with the g++ compiler, including the root libraries:
g++ PhotonCountDistribution.cpp `root-config --libs --cflags` -o PhotonCountDistribution

# Check if the file was successfully compiled:
if [ -e PhotonCountDistribution ]
	then
	echo "Compilation successful."
else
	echo "Compilation error. Check settings and program for errors and try again."
	exit
fi

# Create a list of root files, hand it to the compiled program as an argument and run it:
fileList="";
for file in "../data/"*;do
	fileList+="$file "
done
./PhotonCountDistribution $fileList
