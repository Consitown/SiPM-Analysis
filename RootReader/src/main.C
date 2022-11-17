//root
#include <TString.h>
#include <TH1F.h>

//C, C++
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <map>

//local
#include "geometry.h"
#include "analysis.h"
#include "read.h"

using namespace std;

int main(int argc, char *argv[])
{
  TString inFileList;
  TString inDataFolder;
  TString outFile;

  inFileList = argv[1];
  inDataFolder = argv[2];
  outFile = argv[3];
  string runName = argv[4];
  string headerSize = argv[5];
  if (headerSize == "a")
    headerSize = "-1";
  string dynamicBL = argv[7];
  string isDC = argv[6];
  string useConstCalibValues = argv[8];

  string iWForceRun = argv[16];
  std::string::size_type i = iWForceRun.find("IW_");
  if (i != std::string::npos)
   iWForceRun.erase(i, 3);







  map<string, string> readParameters;
  readParameters.insert(make_pair("inFileList", inFileList));
  readParameters.insert(make_pair("inDataFolder", inDataFolder));
  readParameters.insert(make_pair("outFile", outFile));
  readParameters.insert(make_pair("runName", runName));
  readParameters.insert(make_pair("headerSize", headerSize));
  readParameters.insert(make_pair("dynamicBL", dynamicBL));
  readParameters.insert(make_pair("isDC", isDC));
  readParameters.insert(make_pair("useConstCalibValues", useConstCalibValues));
  if (argc >= 10)
    readParameters.insert(make_pair("runNumber", string(argv[9])));
  if (argc >= 12)
    readParameters.insert(make_pair("runPosition", string(argv[11])));
    // std::cout<<"RUNPOSITION"<<string(argv[11]);
  if (argc >= 13)
    readParameters.insert(make_pair("runAngle", string(argv[12])));
  if (argc >= 14)
    readParameters.insert(make_pair("runEnergy", string(argv[13])));
  if (argc >= 15)
    readParameters.insert(make_pair("runChannelNumberWC", string(argv[14])));
 if (argc >= 16)
    readParameters.insert(make_pair("automaticWindow", string(argv[15])));
  if (argc >= 17)
    readParameters.insert(make_pair("iWForceRun", iWForceRun));

  read(readParameters);

  return 0;
}