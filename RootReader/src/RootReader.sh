#!/bin/bash

#rm read
# -rpath option might be necessary with some ROOT installations
# -lSpectrum option might be necessary with some ROOT installations
#g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -lSpectrum -o read
#g++ geometry.C read.C analysis.C main.C `root-config --libs --cflags` -lSpectrum -o read
threadsStarted=0
threadCounter=0
counter=0
stopNewThreads="0"

fastRunNumber=0
fastInFolder=0
fastOutfolder=0
fastHeaderSize=0
fastRunDir=0
fastRunName=0
fastLineArr=0

showInformation() {

    echo "  _____             _   _____                _           "
    echo " |  __ \           | | |  __ \              | |          "
    echo " | |__) |___   ___ | |_| |__) |___  __ _  __| | ___ _ __ "
    echo " |  _  // _ \ / _ \| __|  _  // _ \/ _\` |/ _\` |/ _ \ '__|"
    echo " | | \ \ (_) | (_) | |_| | \ \  __/ (_| | (_| |  __/ |   "
    echo " |_|  \_\___/ \___/ \__|_|  \_\___|\__,_|\__,_|\___|_|   "
    echo "                                                         "
echo "   _   _   _   _   _     _   _   _   _   _   _   _  ";
echo "  / \ / \ / \ / \ / \   / \ / \ / \ / \ / \ / \ / \ ";
echo " ( F | i | n | a | l ) ( V | e | r | s | i | o | n )";
echo "  \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/ \_/ ";
    echo "                                                         "


    echo "This tool helps you to start your prefered ROOT reading analysis."

    echo "Made with love by Jan Zimmermann in 2019/2020 (jan4995@gmail.com)"
    echo "-------------------------------------------------------"

}

checkDependencies() {

    PKG_OK=$(dpkg-query -W --showformat='${Status}\n' zenity | grep "ok installed")
    if [ "" == "$PKG_OK" ]; then
        echo Please enter the sudo password to install Zenity

        sudo apt-get install zenity
    else
        echo "Zenity is installed! "
    fi

    BUS_OK=$(dpkg-query -W --showformat='${Status}\n' dbus-x11 | grep "ok installed")
    if [ "" == "$BUS_OK" ]; then
        echo Please enter the sudo password to install DBUSx11
        sudo apt-get install dbus-x11
    else
        echo "dbus-x11 is installed! "
    fi
}

chooseOutFolder() {
    checkDependencies
    outFolder=$(zenity --file-selection --directory --title "Select output Folder (.bin Files)?")
    #  echo $outFolder
}

chooseInFolder() {
    checkDependencies
    inFolder=$(zenity --file-selection --directory --title "Select input Folder (.bin Files)?")
}

createInAndOutFolder() {
    if [ ! -d "$1" ]; then
        mkdir "$1"
    fi
    if [ ! -d "$2" ]; then
        mkdir "$2"
    fi
}

parseRunList() {
    dir_name=$0
    rl_file=$1
    nameSchema=$2
    # construct array. contains the elements separated by "_" delimiter in $dir_name
    IFS="_" read -r -a fields <<<"$dir_name"

    nfields=${#fields[@]}

    # get run number
    runNr=${fields[0]}
    if [ "$nameSchema" == "2019" ]; then

        # get pdgID, beam energy
        #11_pos3_angle0_e26_ch32
        # pos=$(echo ${fields[1]} | cut -c 4-)
        pos=$(echo ${fields[1]} | cut -c 4-)
        angle=$(echo ${fields[2]} | cut -c 6-)
        energy=$(echo ${fields[3]} | cut -c 2-)
        channel=$(echo ${fields[4]} | cut -c 3-)
        #  echo "Runlist (2019) creating..."
        line_to_runlsit="$runNr $dir_name $pos $angle $energy $channel"

    else

        runParticle=${fields[1]}
        case $runParticle in
        muon[6])
            pdgID="-13"
            energy="6"
            ;;
        pion[1-6])
            pdgID="211"
            energy=$(echo $runParticle | cut -c 5)
            ;;
        e5)
            pdgID="-11"
            energy=5
            ;;
        *)
            echo "UNKNOWN particle description in run $runNr"
            ;;
        esac

        # get measurement position
        runMP=${fields[2]}
        case $runMP in
        pos[0-9] | pos[1][0-9]) # 0°/30° measurements, no xy-coordinates
            frontMP=$(echo $runMP | cut -c 4-)
            ;;
        scanBD) # 90° WOM scans, hardcoded $sideMP
            side_pos_x=${fields[3]}
            side_pos_y=${fields[4]}
            case $side_pos_y in
            0)
                case $side_pos_x in
                0) sideMP=18 ;; 6) sideMP=19 ;; 12) sideMP=20 ;; 18) sideMP=21 ;; 24) sideMP=22 ;;
                esac
                ;;
            1)
                case $side_pos_x in
                0) sideMP=23 ;; 6) sideMP=24 ;; 12) sideMP=25 ;; 18) sideMP=26 ;; 24) sideMP=27 ;;
                esac
                ;;
            2)
                case $side_pos_x in
                0) sideMP=28 ;; 6) sideMP=29 ;; 12) sideMP=30 ;; 18) sideMP=31 ;; 24) sideMP=32 ;;
                esac
                ;;
            3)
                case $side_pos_x in
                0) sideMP=33 ;; 6) sideMP=34 ;; 12) sideMP=35 ;; 18) sideMP=36 ;; 24) sideMP=37 ;;
                esac
                ;;
            esac
            ;;
        *)
            echo "UNKNOWN position description in run $runNr"
            ;;
        esac

        # get angle
        case $nfields in
        3)
            angle="0"
            ;;
        4)
            angle="0"
            ;;
        5)
            angle=$(echo "${fields[3]}" | cut -c 6-)
            ;;
        7)
            angle=$(echo "${fields[5]}" | cut -c 6-)
            ;;
        *)
            echo "UNKNOWN angele description in run $runNr"
            ;;

        esac

        ######################
        ## PRINT TO RUNLIST ##
        ######################
        # line_to_runlsit
        case $nfields in
        [3-5])
            line_to_runlsit="$runNr $dir_name $frontMP $pdgID $energy $angle"
            ;;
        7)
            line_to_runlsit="$runNr $dir_name $sideMP $pdgID $energy $angle $side_pos_x $side_pos_y"
            ;;
        esac
        echo "Runlist (2018) creating..."
    fi

    #echo "$line_to_runlsit"
    echo "$line_to_runlsit" >>"$rl_file"

}

#g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -lSpectrum -o read
#g++ geometry.C read.C analysis.C main.C `root-config --libs --cflags` -lSpectrum -o read
compileRead() {
    echo "Compiling Read.C..."
    g++ ./src/geometry.C ./src/read.cpp ./src/analysis.C ./src/main.C ./src/misc.C `root-config --libs --cflags` -lSpectrum -o ./src/read
    echo "Compiling done!"
}

readFull() {
    if [ "$shouldCompile" = true ]; then
        compileRead
    fi
    runNr=$1
    readAll=false
    inFolder=$2
    outFolder=$3
    headerSize=$4
    saveFolder=$outFolder

    echo "$runNr"
    echo "$2"
    echo "$3"
    echo "$4"
    echo "$saveFolder"

    if [[ " ${runNumber[@]} " =~ "a" ]]; then
        readAll=true
    fi

    mkdir "$saveFolder"

    while read line; do
        [[ $line == \#* ]] && continue
        lineArr=($line)
		echo HIERHER
		echo $line

        doRun="0"
        if [[ " ${runNumber[@]} " =~ " ${lineArr[0]} " ]]; then
            # echo "PRINT 2: " ${runNumber[@]} ${lineArr[0]}
            doRun="1"
        fi

        if [ "$doRun" = "1" ] || [ $readAll = true ]; then
            runName=${lineArr[1]}
            runDir=$saveFolder/$runName
            mkdir "$runDir"
            if [ ! -e runDir/$runName.list ]; then
                ls $inFolder/$runName | grep \.bin >$runDir/$runName.list
            fi

            #time $here/readFull $runDir/$runName.list $inFolder/$runName/ $runDir/out.root ${lineArr[0]} ${lineArr[2]} ${lineArr[3]} ${lineArr[4]} ${lineArr[5]} ${lineArr[6]}

            time ./src/read $runDir/$runName.list $inFolder/$runName/ $runDir/$runName.root $runName $headerSize "$isDC" "$dynamicBL" "$useCalibValues" "${lineArr[0]}" "${lineArr[1]}" "${lineArr[2]}" "${lineArr[3]}" "${lineArr[4]}" "${lineArr[5]}" "$automaticWindow" "$iWForceRun"

        fi
    done <./RootRunlist.txt

}

compileMerger() {
    echo "Compiling Root Merger..."
    g++ ./src/mergeROOTFiles.C $(root-config --libs --cflags) -lSpectrum -o ./src/mergeROOTFiles
    echo "Compiling done!"
}

merger() {
    rootFileList=$1
    ./src/mergeROOTFiles $rootFileList $2 $3

    #merge PDF Waves.pdf
    pdfunite $(find $2 -name "*waveforms.pdf") $2/waveforms.pdf
    pdfunite $(find $2 -name "*waveforms_chSum.pdf") $2/waveforms_chSum.pdf
    pdfunite $(find $2 -name "*waveforms_womSum.pdf") $2/waveforms_womSum.pdf

}

installPDFUnite() {
    PKG_OK=$(dpkg-query -W --showformat='${Status}\n' poppler-utils | grep "ok installed")
    if [ "" == "$PKG_OK" ]; then
        echo Please enter the sudo password to install PDF Unite to create PDF files

        sudo apt-get install poppler-utils
    else
        echo "PDFUnite is installed! "
    fi
}

readFast() {
    if [ "$shouldCompile" = true ]; then
        compileRead
    fi
    installPDFUnite
    compileMerger
    readAll=false
    inFolder=$2
    outFolder=$3
    headerSize=$4
    saveFolder=$outFolder

    if [[ " ${runNumber[@]} " =~ "a" ]]; then
        readAll=true
    fi

    mkdir "$saveFolder"

    rootFileList=""

    # " ${runNr[@]} " =~ " ${lineArr[0]} "

    time (
        while read line; do
            # echo $line
            [[ $line == \#* ]] && continue #Wenn die Zeile leer ist wird "continue" ausgeführt
            lineArr=($line)

            doRun="0"
            if [[ " ${runNumber[@]} " =~ " ${lineArr[0]} " ]]; then
                # echo "PRINT 2: " ${runNumber[@]} ${lineArr[0]}
                doRun="1"
            fi

            if [ "$doRun" = "1" ] || [ $readAll = true ]; then #lineArr ist ein Array aus jeder Zeile, getrennt durch Leerzeichen. linearray[0] ist die Run NUMMER
                rootFileList=""
                runName=${lineArr[1]}
                runDir=$saveFolder/$runName
                mkdir "$runDir"
                rm -rf $runDir/*
                if [ ! -e $runDir/$runName.list ]; then
                    ls $inFolder/$runName | grep \.bin >$runDir/$runName.list #Durchsucht das RunName verzeichnis nach bins files und erstellt eine Runlist
                fi

                counter=0
                threadCounter=0

                #Für jede Bin einzeln durchlaufen
                while read lineTemp; do

                    #    ls *.cfg | xargs -P 4 -n 1 read_cfg.sh
                    #   Create RunList for every file

                    echo $lineTemp >$runDir/$counter.list
                    mkdir $runDir/$counter
                    runDirRelative="${runDir//$here/}"
                    rootTreeFilePath=".$runDirRelative/$counter/out.root/T"
                    rootFileList="${rootFileList}||$rootTreeFilePath"
                    counter=$((counter + 1))

                done \
                    <$runDir/$runName.list

                fastInFolder=$inFolder
                fastOutFolder=$outFolder
                fastHeaderSize=$headerSize
                fastRunDir=$runDir
                fastLineArr=($line)
                fastRunName=${fastLineArr[1]}
                fastRunNumber=$runNr

                remainder=$((counter % threads))
                loopNumber=$((counter / threads))
                #loopNumber=$((loopNumber - 1)) #to have to correct number in the loop
                # if [ "$remainder" != "0" ]; then
                #   loopNumber=$((loopNumber + 1))
                # fi

                echo "Loop Number: $loopNumber Number of Files: $counter Extra Threads: $remainder Threads: $threads"
                for currentBin in $(seq 1 $loopNumber); do
                    echo "Main Thread Loop Started"

                    for ((i = 0; i < "$threads"; i++)); do
                        readFastIteration $threadCounter &
                        threadCounter=$((threadCounter + 1))
                    done

                    wait

                done
                for currentBinRemain in $(seq 1 $remainder); do
                    echo "Remainder Thread Started"
                    readFastIteration $threadCounter &
                    threadCounter=$((threadCounter + 1))

                done

                wait

                merger $rootFileList $runDir $runName
                #rm $runDir/*.list
                rm -rf $runDir/*/
                find $runDir -name "*.list" -type f -delete

            fi

        done \
            <./RootRunlist.txt

        #reads lines of "runslist"

    )

}

readRoot() {

    threads=$1
    runNumber=$2 #a=ALL
    inFolder=$3
    outFolder=$4
    headerSize=$5
    if [ "$threads" -gt "1" ]; then
        readFast $runNumber $inFolder $outFolder $headerSize
    else
        readFull $runNumber $inFolder $outFolder $headerSize
    fi

}

readFastIteration() {
	echo "readFastITeration "
    ./src/read $fastRunDir/$1.list $fastInFolder/$fastRunName/ $fastRunDir/$1/out.root $fastRunName $fastHeaderSize "$isDC" "$dynamicBL" "$useCalibValues" "${fastLineArr[0]}" "${fastLineArr[1]}" "${fastLineArr[2]}" "${fastLineArr[3]}" "${fastLineArr[4]}" "${fastLineArr[5]}" "$automaticWindow" "$iWForceRun"

}

saveConfig() {
    rm config.txt 2>/dev/null
    destdir=config.txt
    echo "$inFolder" >>"$destdir"
    echo "$outFolder" >>"$destdir"
    echo "$threads" >>"$destdir"
    echo "${runNumber[@]}" >>"$destdir"
    echo "$headerSize" >>"$destdir"
    echo "$isDC" >>"$destdir"
    echo "$dynamicBL" >>"$destdir"
    echo "$useCalibValues" >>"$destdir"
    echo "$useExistingRunList" >>"$destdir"
    echo "$automaticConfig" >>"$destdir"
    echo "$automaticWindow" >>"$destdir"

    echo "config saved!"

}

saveAnalysisPath() {
    rm analysisPath.txt 2>/dev/null
    destdir=analysisPath.txt
    echo "$analysisPath" >>"$destdir"
}
loadAnalysisPath() {
    if test -f "analysisPath.txt"; then
        if ! test -z "$(awk 'NR == 1' analysisPath.txt)"; then
            analysisPath=$(awk 'NR == 1' analysisPath.txt)
        fi

    fi
}

checkAutomaticConfig() {
    if test -f "config.txt"; then
        tempConfig=$(awk 'NR == 10' config.txt)
    fi
    if [ "$tempConfig" = true ]; then
        loadConfig

    fi

}

loadConfig() {
    #destdir=config.txt
    #$inFolder=$( sed -n '1 p' config.txt )

    if test -f "config.txt"; then
        inFolder=$(awk 'NR == 1' config.txt)
        outFolder=$(awk 'NR == 2' config.txt)
        threads=$(awk 'NR == 3' config.txt)
        runNumberTemp=$(awk 'NR == 4' config.txt)
        runNumber=($runNumberTemp)
        headerSize=$(awk 'NR == 5' config.txt)
        isDC=$(awk 'NR == 6' config.txt)
        dynamicBL=$(awk 'NR == 7' config.txt)
        useCalibValues=$(awk 'NR == 8' config.txt)
        useExistingRunList=$(awk 'NR == 9' config.txt)
        automaticConfig=$(awk 'NR == 10' config.txt)
        automaticWindow=$(awk 'NR == 11' config.txt)
        echo "config loaded!"
    else
        echo "config file not found! Save a config first!"

    fi
}

startAutomaticEffRun() {
    echo "Start? (Run:  ${runNumber[*]})"
    select yn in "Yes" "No"; do
        case $yn in

        \
            "No")
            break 2
            ;;

        "Yes")
            echo "Started for runs: (${runNumber[*]}) "
            echo "Select DC Run:"
            read runNumberDC
            iwScriptDir=$(find $analysisPath -name 'IntegrationWindowAnalysis.sh' -printf "%h\n")
            rootfileFolderDir=$analysisPath/rootfiles
            rootfileFolderDirFinished=$analysisPath/finishedRootfiles
            
            mkdir -p $rootfileFolderDirFinished;

            echo "Location of the Integration Scripts $iwScriptDir"
            echo "Location of the Rootfile Folder  $rootfileFolderDir"

            echo "Make sure there are already the correct rootfiles in the analysis/rootfiles folder or press copy to move them from the runfolder to /rootfiles!"
            echo "Yes -> starts automatic Baseline File creation"
            select yn in "Create" "Copy" "Continue"; do
                case $yn in
                "Create")
                    echo "Start creating Rootfiles with a constant Baseline? (Run: ${runNumber[*]})"
                    startAutomaticBaseline
                    break
                    ;;
                "Continue")
                    break
                    ;;
                "Copy")
                    for n in ${runNumber[@]}; do
                        find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
                    done
                    break
                    ;;
                esac
            done

            echo "Calculating the Integration windows by the sum Histograms (Check the histograms)"
            ($iwScriptDir/IntegrationWindowAnalysis.sh)
            echo "Done Integration Window Script! "
            echo "Moving IntegrationWindows.txt to RunHelper/src"
            cwd=$(pwd)
            find $iwScriptDir -type f -name "IntegrationWindows.txt" -exec cp {} $cwd/src/ \;

            automaticWindow=true #make sure this is on
            echo "Run with new Integration Windows... "
            start
            echo "Moving files..."
            for n in ${runNumber[@]}; do
                find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
            done

            echo "Calculating the Correction Factor..."
            ($iwScriptDir/IntegrationWindowAnalysis.sh)
            echo "Done Integration Window Script! "
            echo "Moving CorrectionValues.txt to RunHelper/src"
            cwd=$(pwd)
            find $iwScriptDir -type f -name "CorrectionValues.txt" -exec cp {} $cwd/src/ \;

            echo "Run with Correction Factor... "
            start
            echo "Moving IW files..."
            for n in ${runNumber[@]}; do
                find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} "$rootfileFolderDirFinished" \;
                find $rootfileFolderDir -type f -name "${n}_*" -a -name '*.root' -exec rm -rf {} \;
            done

            echo "Automatic Integration Done! Moved files to $rootfileFolderDirFinished"

            dcScriptDir=$(find $analysisPath -name 'DCProbability.py' -printf "%h\n")

            rootfileFolderDir=$analysisPath/rootfiles
            savedRunNumber=("${runNumber[@]}")

            unset runNumber
            runNumber=("${runNumberDC[@]}")

            for wRun in "${savedRunNumber[@]}"; do
                #find runname
                aRunName=$(grep -I "IW_${wRun}_" ./src/IntegrationWindows.txt | tail -1 | cut -f1 -d"=")
                echo "Using $aRunName"
                iWForceRun=$aRunName

                automaticWindow=false
                #start DC reading with IW from aRunName
                start
                echo "Moving DC files..."

                suffix="iw$wRun"
                path=$(find . -type f -name "${runNumber}_dc*" -a -name '*.root')
                echo "origin path: $path, destination: $rootfileFolderDir/${runNumber}_dc_${suffix}.root"
                if cp $path "$rootfileFolderDir/${runNumber}_dc_${suffix}.root"; then
                    echo "Moving successfull!"
                else
                    echo "Moving not successfull! :("
                    exit 222
                fi

            done
            runNumber=("${savedRunNumber[@]}")
            echo "Done with the automatic efficiency run!"

            break
            ;;
        esac
    done

}

startAutomaticDCRun() {
    #Goal: Calculate 2 DC Limits, with 1 calib IW and measurement IW -> 2 DC Analysis, with different IW

    dcScriptDir=$(find $analysisPath -name 'DCProbability.py' -printf "%h\n")

    rootfileFolderDir=$analysisPath/rootfiles

    echo "Using runs: ${runNumber[*]}"
    savedRunNumber=("${runNumber[@]}")
    unset runNumber
    echo "Select DC Run:"
    read runNumber

    useCalibValues="1"

    for wRun in "${savedRunNumber[@]}"; do
        #find runname
        aRunName=$(grep -I "IW_${wRun}_" ./src/IntegrationWindows.txt | tail -1 | cut -f1 -d"=")
        echo "Using $aRunName"
        iWForceRun=$aRunName
        automaticWindow=false
        #start DC reading with IW from aRunName
        echo "DEBUG: $runNumber    ${savedRunNumber[*]}"
        start
        echo "Moving files..."
        shouldCompile=false
        suffix="iw$wRun"
        path=$(find . -type f -name "${runNumber}_dc*" -a -name '*.root')
        echo "path: $path"
        if cp $path "$rootfileFolderDir/${runNumber}_dc_${suffix}.root"; then
            echo Moving successfull!
        else
            exit 222
        fi

    done
    runNumber=("${savedUunNumber[@]}")
    automaticWindow=true

    #python $dcScriptDir/DCProbability.py

}

startAutomaticBaseline() {
    while true; do
        select yn in "Yes" "No"; do
            case $yn in

            \
                "No")
                break 2
                ;;

            "Yes")
                echo "Started for runs: (${runNumber[*]}) "

                baselineScriptDir=$(find $analysisPath -name 'BaselineAverage.sh' -printf "%h\n")
                rootfileFolderDir=$analysisPath/rootfiles

                echo "Location of the Baseline Script $baselineScriptDir"
                echo "Location of the Rootfile Folder  $rootfileFolderDir"

                echo "Starting Runs with dynamic Baseline..."
                dynamicBL="1"
                start

                echo "Moving files..."

                for n in ${runNumber[@]}; do
                    find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
                done

                echo "Starting Baseline Script"
                ($baselineScriptDir/BaselineAverage.sh)
                echo "Done Baseline Script! "
                echo "Moving Baselines.txt to RunHelper/src"
                cwd=$(pwd)
                find $baselineScriptDir -type f -name "Baselines.txt" -exec cp {} $cwd/src/ \;

                echo "Starting Runs with constant Baseline..."
                dynamicBL="0"
                start

                echo "Moving files..."

                for n in ${runNumber[@]}; do
                    find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
                done
                echo "Automatic Baseline done!"

                break 2
                ;;
            esac
        done
    done

}

startAutomaticIntegration() {

    echo "Start? (Run:  ${runNumber[*]})"
    select yn in "Yes" "No"; do
        case $yn in

        \
            "No")
            break 2
            ;;

        "Yes")
            echo "Started for runs: (${runNumber[*]}) "

            iwScriptDir=$(find $analysisPath -name 'IntegrationWindowAnalysis.sh' -printf "%h\n")
            rootfileFolderDir=$analysisPath/rootfiles
            rootfileFolderDirFinished=$analysisPath/finishedRootfiles

            mkdir -p $rootfileFolderDirFinished;
            echo "Make sure there are already the correct rootfiles in the analysis/rootfiles folder or press copy to move them from the runfolder to /rootfiles!"
            echo "Yes -> starts automatic Baseline File creation"
            select yn in "Create" "Copy" "Continue"; do
                case $yn in
                "Create")
                    echo "Start creating Rootfiles with a constant Baseline? (Run: ${runNumber[*]})"
                    startAutomaticBaseline
                    break
                    ;;
                "Continue")
                    break
                    ;;
                "Copy")
                    for n in ${runNumber[@]}; do
                        find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
                    done
                    break
                    ;;
                esac
            done

            echo "Location of the Integration Scripts $iwScriptDir"
            echo "Location of the Rootfile Folder  $rootfileFolderDir"

            echo "Calculating the Integration windows by the sum Histograms (Check the histograms)"
            ($iwScriptDir/IntegrationWindowAnalysis.sh)
            echo "Done Integration Window Script! "
            echo "Moving IntegrationWindows.txt to RunHelper/src"
            cwd=$(pwd)
            find $iwScriptDir -type f -name "IntegrationWindows.txt" -exec cp {} $cwd/src/ \;

            automaticWindow=true #make sure this is on
            echo "Run with new Integration Windows... "
            start
            echo "Moving files..."
            for n in ${runNumber[@]}; do
                find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} $rootfileFolderDir \;
            done

            echo "Calculating the Correction Factor..."
            ($iwScriptDir/IntegrationWindowAnalysis.sh)
            echo "Done Integration Window Script! "
            echo "Moving CorrectionValues.txt to RunHelper/src"
            cwd=$(pwd)
            find $iwScriptDir -type f -name "CorrectionValues.txt" -exec cp {} $cwd/src/ \;

            echo "Run with Correction Factor... "
            start
            echo "Moving files..."
            for n in ${runNumber[@]}; do
                find . -type f -name "${n}_*" -a -name '*.root' -exec cp {} "$rootfileFolderDirFinished" \;
                find $rootfileFolderDir -type f -name "${n}_*" -a -name '*.root' -exec rm -rf {} \;
            done


            echo "Automatic Integration Done! Moved to: $rootfileFolderDirFinished"

            break
            ;;
        esac
    done

}

start() {
    echo "   _____ _             _   "
    echo "  / ____| |           | |  "
    echo " | (___ | |_ __ _ _ __| |_ "
    echo "  \___ \| __/ _\` | '__| __|"
    echo "  ____) | || (_| | |  | |_ "
    echo " |_____/ \__\__,_|_|   \__|"
    echo "                           "
    echo "                           "
    echo "--------------------------------PARAMETER--------------------------------"

    echo "Input data folder: $inFolder"
    echo "Output data folder: $outFolder"
    echo "Compiles the scripts: $shouldCompile"
    echo "Threads: $threads"
    echo "RunNumber: ${runNumber[*]}"
    echo "Use Existing RunList: $useExistingRunList"
    echo "Header Size: $headerSize"
    echo "Automatic Config: $automaticConfig"
    echo "---------------------------------------------"

    if [ "$automaticConfig" = true ]; then
        saveConfig
    fi

    createInAndOutFolder $inFolder $outFolder

    if [ "$useExistingRunList" = false ]; then
        echo "Creating a runlist for all input folder Files... (Name: RootRunlist.txt)"
        export -f parseRunList
        runlist="RootRunlist.txt"
        rm "$runlist"
        ls $inFolder | sort -n | xargs -n 1 -P 1 bash -c "parseRunList $runlist $nameSchema"

        #Check if readList exists
        if [ -f $runList ]; then
            echo "Runlist created successfully!"
        else
            echo "Runlist creation failed!"
        fi
    fi

    readRoot $threads $runNumber $inFolder $outFolder $headerSize
    echo "---------------------------------------------"
    echo "  ____    ___   _   _  _____ "
    echo " |  _ \  / _ \ | \ | || ____|"
    echo " | | | || | | ||  \| ||  _|  "
    echo " | |_| || |_| || |\  || |___ "
    echo " |____/  \___/ |_| \_||_____|"
    echo "                             "
    echo "---------------------------------------------"
}

# ██████╗ ██████╗ ██████╗ ███████╗
#██╔════╝██╔═══██╗██╔══██╗██╔════╝
#██║     ██║   ██║██║  ██║█████╗
#██║     ██║   ██║██║  ██║██╔══╝
#╚██████╗╚██████╔╝██████╔╝███████╗
# ╚═════╝ ╚═════╝ ╚═════╝ ╚══════╝

showInformation

here=$(pwd)
src=$here/src/
inFolder=$here/data
outFolder=$here/runs
shouldCompile=true
automaticConfig=true
threads=6
runNumber=a
headerSize=a
useExistingRunList=false
automaticWindow=true
useCalibValues="0"
dynamicBL="1"
isDC="0"
nameSchema="2019"
checkAutomaticConfig
loadAnalysisPath

echo "-------------------------------------------------------"

echo "Select an option!"

while true; do
    echo "-------------------------------------------------------"
    echo "-------------------------------------------------------"
    echo "-------------------------------------------------------"

    # select yn in "Start" "RunMode ($threads)" "RunNumber (${runNumber[*]})" "Change Input Folder ($inFolder)" "Change Output Folder ($outFolder)" "Compile Readscripts ($shouldCompile)" "Binary Headersize ($headerSize)" "Calibration/Baseline" "Use existing RunList ($useExistingRunList)" "RunList Naming Schema ($nameSchema)" "Load Config" "Save Config" "Informations"; do
    select yn in "Start" "RunNumber (${runNumber[*]})" "Settings" "Tools" "Load Config" "Save Config" "Informations"; do
        case $yn in

        \
            \
            "Settings")
            while true; do
                select yn in "Threads ($threads)" "Automatic Integration Window ($automaticWindow)" "Automatic Config Saving ($automaticConfig)" "Compile Readscripts ($shouldCompile)" "Binary Headersize ($headerSize)" "Calibration/Baseline" "Use existing RunList ($useExistingRunList)" "RunList Naming Schema ($nameSchema)" "Change Input Folder ($inFolder)" "Change Output Folder ($outFolder)" "<- Back"; do
                    case $yn in

                    \
                        "<- Back")
                        break 2
                        ;;
                    "Threads ($threads)")
                        echo "How many threads should run parallel? (1 Thread -> Slow, but no ROOT File merging), 8 CORE CPU -> Max 8 Threads, 1 Thread takes about 1GB RAM"
                        read threads
                        break
                        ;;
                    "Automatic Integration Window ($automaticWindow)")
                        echo "Uses the values stored in the IntegrationWindow.txt File, generated by the analysis tool 'IntegrationWindowAnalysis' , if false -> hardcoded values in read.C"
                        if [ "$automaticWindow" = false ]; then
                            automaticWindow=true
                        else
                            automaticWindow=false
                        fi
                        break
                        ;;

                    "Automatic Config Saving ($automaticConfig)")
                        echo "Saves and restores the latest started run config! (Changing this automatically saves the config)"
                        if [ "$automaticConfig" = false ]; then
                            automaticConfig=true
                        else
                            automaticConfig=false
                        fi
                        saveConfig
                        break
                        ;;

                    \
                        "Change Input Folder ($inFolder)")
                        chooseInFolder
                        break
                        ;;
                    "Binary Headersize ($headerSize)")
                        echo "Enter the headersize of the binary files (Depends on the WC version)"
                        echo "Version: <2.8.14: 327,  ==2.8.14: 328,  >2.8.14: 403   -> a=automatic"

                        read headerSize
                        break
                        ;;

                    "Change Output Folder ($outFolder)")
                        chooseOutFolder
                        break
                        ;;
                    "Compile Readscripts ($shouldCompile)")
                        if [ "$shouldCompile" = false ]; then
                            shouldCompile=true
                        else
                            shouldCompile=false
                        fi

                        break
                        ;;
                    "Is Dark Count Measurement ($isDC)")
                        if [ "$isDC" = "0" ]; then
                            isDC="1"
                        else
                            isDC="0"
                        fi
                        break
                        ;;
                    "Use existing RunList ($useExistingRunList)")
                        if [ "$useExistingRunList" = false ]; then
                            useExistingRunList=true

                        else
                            useExistingRunList=false
                        fi
                        break
                        ;;

                    "RunList Naming Schema ($nameSchema)")
                        echo "Enter Name Schema ($nameSchema)"
                        read nameSchema
                        break
                        ;;
                    "Dynamic Baseline ($dynamicBL)")
                        if [ "$dynamicBL" = "0" ]; then
                            dynamicBL="1"
                        else
                            dynamicBL="0"
                        fi
                        break
                        ;;
                    "Calibrate Measurement ($useCalibValues)")
                        if [ "$useCalibValues" = "0" ]; then
                            useCalibValues="1"
                        else
                            useCalibValues="0"
                        fi
                        break
                        ;;

                    "Calibration/Baseline")
                        while true; do
                            select yn in "Is Dark Count Measurement ($isDC)" "Dynamic Baseline ($dynamicBL)" "Calibrate Measurement ($useCalibValues)" "<- Back"; do
                                case $yn in

                                \
                                    "<- Back")
                                    break 2
                                    ;;

                                "Is Dark Count Measurement ($isDC)")
                                    if [ "$isDC" = "0" ]; then
                                        isDC="1"
                                    else
                                        isDC="0"
                                    fi
                                    break
                                    ;;
                                "Dynamic Baseline ($dynamicBL)")
                                    if [ "$dynamicBL" = "0" ]; then
                                        dynamicBL="1"
                                    else
                                        dynamicBL="0"
                                    fi
                                    break
                                    ;;
                                "Calibrate Measurement ($useCalibValues)")
                                    if [ "$useCalibValues" = "0" ]; then
                                        useCalibValues="1"
                                    else
                                        useCalibValues="0"
                                    fi
                                    break
                                    ;;
                                esac
                            done
                        done

                        break
                        ;;
                    esac
                done
            done
            break
            ;;

        \
            "Start")
            printf "\033c" #Clean Screen
            start

            break 2
            ;;

        \
            "RunNumber (${runNumber[*]})")
            echo "Enter a runNumber (ALL=a, Example: 22 OR multiple: 31,32,40 OR range: 10-15 OR combined)"
            read runNumberRaw
            runNumberRaw=($(echo "$runNumberRaw" | tr ',' '\n'))
            # echo "${runNumber[*]}"

            for i in "${runNumberRaw[@]}"; do
                :
                if [[ $i == *"-"* ]]; then
                    range=($(echo $i | tr "-" "\n"))
                    startValue=${range[0]}
                    endValue=${range[1]}

                    while [ $startValue -le $endValue ]; do
                        runNumberRaw=("$startValue" "${runNumberRaw[@]}")
                        startValue=$(($startValue + 1))
                    done
                fi

                #Remove Range
                unset runNumber
                for z in "${runNumberRaw[@]}"; do
                    :
                    if [[ $z != *"-"* ]]; then
                        runNumber=("$z" "${runNumber[@]}")
                    fi
                done
            done
            break
            ;;

        "Tools")
            while true; do
                select yn in "Automatic Baseline Run" "Automatic Integration Run" "Automatic Dark Count Runs" "Automatic Efficiency Run" "<- Back"; do
                    case $yn in

                    \
                        "<- Back")
                        break 2
                        ;;

                    "Automatic Baseline Run")
                        echo "This will start an automatic run, creating a root file with a constant baseline inside the run folder"
                        if [ -z ${analysisPath+x} ]; then
                            echo "Your analysis script path is not set. Do it now:"
                            analysisPath=$(zenity --file-selection --directory --title "Select Analysis Script path")
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                saveAnalysisPath
                            fi

                            echo "Your analysis script path is: $analysisPath"
                            echo "If you want to change it, just delete or edit the file: analysisPath.txt"
                            echo "Start? (Run: ${runNumber[*]})"
                            startAutomaticBaseline

                        else
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                echo "Your analysis script path is: $analysisPath"
                                echo "If you want to change it, just delete or edit the file: analysisPath.txt"
                                echo "Start? (Run: ${runNumber[*]})"
                                startAutomaticBaseline
                            fi

                        fi

                        # echo "Your analysis script path is set to: $analysisPath"

                        #outFolder=$(zenity --file-selection --directory --title "Select output Folder (.bin Files)?")

                        break 2
                        ;;
                    "Automatic Integration Run")
                        echo "This will start an automatic run, creating a root file with correct integrals"
                        if [ -z ${analysisPath+x} ]; then
                            echo "Your analysis script path is not set. Do it now:"
                            analysisPath=$(zenity --file-selection --directory --title "Select Analysis Script path")
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                saveAnalysisPath
                            fi

                            echo "Your analysis script path is: $analysisPath"
                            echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                            startAutomaticIntegration

                        else
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                echo "Your analysis script path is: $analysisPath"
                                echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                                startAutomaticIntegration
                            fi

                        fi

                        # echo "Your analysis script path is set to: $analysisPath"

                        #outFolder=$(zenity --file-selection --directory --title "Select output Folder (.bin Files)?")

                        break 2
                        ;;
                    "Automatic Dark Count Runs")
                        echo "This will start an automatic run, calculating Dark Count Limits for multiple integration windows!"
                        if [ -z ${analysisPath+x} ]; then
                            echo "Your analysis script path is not set. Do it now:"
                            analysisPath=$(zenity --file-selection --directory --title "Select Analysis Script path")
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                saveAnalysisPath
                            fi

                            echo "Your analysis script path is: $analysisPath"
                            echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                            startAutomaticDCRun

                        else
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                echo "Your analysis script path is: $analysisPath"
                                echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                                startAutomaticDCRun
                            fi

                        fi

                        # echo "Your analysis script path is set to: $analysisPath"

                        #outFolder=$(zenity --file-selection --directory --title "Select output Folder (.bin Files)?")

                        break 2
                        ;;
                    "Automatic Efficiency Run")
                        echo "This will start an automatic run, calculating DC Limit and IW."
                        if [ -z ${analysisPath+x} ]; then
                            echo "Your analysis script path is not set. Do it now:"
                            analysisPath=$(zenity --file-selection --directory --title "Select Analysis Script path")
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                saveAnalysisPath
                            fi

                            echo "Your analysis script path is: $analysisPath"
                            echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                            startAutomaticEffRun

                        else
                            if test -z "$analysisPath"; then
                                echo "Analysis Path not set!"
                            else
                                echo "Your analysis script path is: $analysisPath"
                                echo "If you want to change it, just delete or edit the file: analysisPath.txt"

                                startAutomaticEffRun
                            fi

                        fi

                        # echo "Your analysis script path is set to: $analysisPath"

                        #outFolder=$(zenity --file-selection --directory --title "Select output Folder (.bin Files)?")

                        break 2
                        ;;
                    esac

                done
            done

            break
            ;;

        "Load Config")
            loadConfig
            break
            ;;
        "Save Config")
            saveConfig
            break
            ;;
        "Informations")
            echo "-------------------------------------------------------"
            echo "Thread Informations: "
            echo "1 Thread:  High memory usage, very slow. But there is no merging of individual rootfiles, which could probably be an error source."
            echo "More than 1 Thread: Creats multiple Rootfiles at the same time and merges them together with TChain, Histograms are missing if not included in MergeRootFiles.C"
            echo "RunList: Has to be in the same directory as this script, named: RootRunlist, Content like: 22 22_muon6_pos4 4 -13 6 0 AB per Line"
            echo "-------------------------------------------------------"
            echo "Header: The header size (open binary file with editor and count chars at the top)"
            echo "changes from version to version. And also from different WC channel sizes"
            echo "There are some changes that have to be made in the read.C scripts to use a new header size"
            echo "Changes like reading new data in the event structure or reading a dummy byte!"
            echo "Just put headersize on 'a' to let the system detect the headersize by reading the printable chars at the .bin beginning"
            echo "Maybe helping: https://owncloud.lal.in2p3.fr/index.php/s/3c1804f473fc1b2d3de454a23329b4a7"
            echo "If there are still questions, feel free to contact me: jan4995@gmail.com"
            echo "-------------------------------------------------------"

            break
            ;;

        esac

    done

done

# sudo apt-get install zenity
#FILE=$(dialog --title "Delete a file" --stdout --title "Please choose a file to delete" --fselect /tmp/ 14 48)
