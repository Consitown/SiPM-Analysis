# wavecatcher-analysis

To get started download repository and navigate to examples folder. Open CERN ROOT in examples folder and execute:
.x read_exampledata_notcompiled.cc



You can speed up the analysis by compiling library. To compile it on linux or mac do: 
make -f makefile

Then add gSystem->Load("ReadRunLib.sl"); to macro or to rootlogon.C and then execute .x read_exampledata.cc


Needs CERN ROOT Release 6.24/02 or higher compiled with c++17 (otherwise need to change line 15 in makefile --std=c++17 -> --std=c++14)
