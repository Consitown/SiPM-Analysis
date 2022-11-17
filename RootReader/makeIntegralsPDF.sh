#!/bin/bash

FILES="/mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/finishedRootfiles/*"
for file in $FILES;do
    echo $file
done
saveFolder=mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/integralAnalysis
running=true
# Ask user, which file(s) should get analysed
while $running; do
    echo "Enter run number to analyse (for all enter 'a'). 'q' to quit." # XXX
    read proceed
    case $proceed in
        [a]* ) runNumber=$a;;
        [q]* ) echo "quit!"; exit;;
        *) runNumber=$proceed; echo "Runnumber is ${runNumber}!";;
    esac

    echo "Start analysis? [y/n]"
    read answer
    case $answer in
        [Yy]* ) echo "Proceeding..."; break;;
		[Nn]* ) echo "Program aborted."; running=false; exit;;
		* ) echo "Please type Y/y or N/n.";;
    esac
    
done

# compile makeIntegralsPDF.cpp with libraries
echo "Compiling!"
g++ ./src/integralAnalysisMain.cpp -L/opt/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/opt/root/lib -pthread -lm -ldl -rdynamic -pthread -std=c++14 -m64 -I/opt/root/include -o ./src/integralAnalysisMain
g++ ./src/meanAnglePlots.cpp -L/opt/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/opt/root/lib -pthread -lm -ldl -rdynamic -pthread -std=c++14 -m64 -I/opt/root/include -o ./src/meanAnglePlots

# Check if the file was successfully compiled:
if [ -e ./src/makeIntegralsPDF ]
	then
	echo "Compilation successful."
else
	echo "Compilation error. rip"
	exit
fi

file=$(find /mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/finishedRootfiles -maxdepth 1 -name "${runNumber}*" )
# echo "file here!"
# echo $file

filename=$(basename -- "$file")
filename="${filename%.*}"
# echo $filename
saveFolder=/mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/integralAnalysis/$filename
mkdir $saveFolder
# echo $saveFolder

# echo $file

./src/integralAnalysisMain "$file"
./src/meanAnglePlots "$file"

python3 plotIntegralMeans.py "${saveFolder}/propagatedintegralMeans.txt" "$saveFolder"





















