# SiPM-Analysis
SiPM = silicon photomultiplier

## RootReader
the original RootReader: https://github.com/Uni2K/RootReader \
modified RootReader by Andrea: https://github.com/andi-matter/bachelor/tree/andrea \
modified RootReader-scripts by Alex: https://github.com/AlexVagts/wavecatcher-analysis (he did not really tidy up his repository) \
I kept all important Macros and changes from Alex/Andrea and I will try to document them here.

## wavecatcher-analysis (analysis_programm)
the original wavecatcher-analysis: https://github.com/cscharf-hub/wavecatcher-analysis \
Since the wavecatcher-analysis's main Script, the "ReadRun.cc", gets updated constantly, I will try to keep it up to date here. More important are the Macros I used to analyse the data.

## General infos
More details of the two pieces of software will be given after this section. \
(I am copying bits and pieces from before mentioned repositories for this, but also give some infos myself) \
If I speak of a Macro, I mean a file that one executes in the terminal. \
If I speak of a Script, I mean a file that is used by a macro (all Scripts are .cpp or .cc or .C - files).

### Some background for context
I used the here presented software to analyse data generated by SiPM-arrays which are placed on PCB's (Printed Circuit Board). These SiPM's have the task to collect light giuded by polymethylmethacryate (PMMA) tubes/cylinders. The light is produced by liquid-szintillator (LS) and the wavelength of the produced light is shifted by the coating on the PMMA tubes/cylinders (so called WOM-tubes; wavelength shifting optical modules). The shifted wavelength yields the highest quantum effeciency for the used SiPM's. The LS (consists of LAB+PPO) is confined inside a steel box and the PMMA tubes/cylinders are emerged into the LS within a PMMA vessel (it's like a pot). These steel boxes are part of the planned SHiP (search for hidden particles) project. In particular the boxes are to be used as the surround background tagger (SBT, basically a massive veto detector). See this pdf for some more infos on the SHiP-projects: https://www.physik.uni-hamburg.de/iexp/gruppe-hagner/pdf/bick-dpg-2017-online.pdf (the SBT is the large, gray, increasing in thickness when going to the right thing on page 21)

### Taking data
The eMUSIC Miniboard v2 by Scientifica (http://siub.ub.edu/down/music/emusic-miniboard.user_guide.r2.5a.pdf) is in charge of collecting the electric signals from the SiPM's. These signals are then digitalized using a WaveCatcher from Caen (https://www.caen.it/products/wavecatcher/) and recorded via the piece of software that accompanies the WaveCatcher. I usually save only the waveforms of the recorded events (enough for all analysis purposes). These waveforms are stored in .bin-files.

### What are .bin files?
Binary files are just holding raw informations, e.g. the voltage values of a measurement. In case of the SHiP HU group, they are created by an analog to digital converter called WaveCatcher. It digitizes the raw SiPM signals and records them using a Windows client. The way the raw data is stored in the .bin files is set in the Windows client. The used configuration puts the 1024 voltage measurements of 1000  events in a single .bin file. The files of the Windows client can be found in the WaveCatcher folder in this repository. The calibration/configuration files are also in the folders (use them to find the exact saving configuration of the bin files). The installation of the WaveCatcher client may requires to restart the PC (it needs to find drivers, which may only be available after a restart).

### What are .root files?
Root files are files created by the CERN ROOT framework (https://root.cern.ch/). They store informations in a hierachical tree-like structure called ROOT tree. Each piece of information is then called leaf, while each type of information is called branch.
Example: You measured the voltage of a photodiode over a timespan of 10ns. A possible branch would then be the amplitude of the voltage, while a leaf holds for example 5mV.
Each Root File contains the information that you want it to contain. In this workgroup (SHiP HU Berlin), we are interested in the time integral of a SiPM signal. We are also interested in several timing obersables. All those branches need to be calculated and stored in the .root file.

### Whats the reason?
The size of the .bin files is usually large, because all the raw data is included. For the analysis, however, you only need specific informations like the amplitude, integral, start point etc. So, the idea is to convert the large amounts of data into single condensed .root files. In this way, the analysis is faster and only relies on some "small" .root files. This method is also more user friendly as it allows a central debugging of the reading process.
