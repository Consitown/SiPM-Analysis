// reads finished Rootfile from (input) folder. Prints Integral histograms for each channel to pdf. Collects Mean Integral values.
// Calculates normalised integral men = integral mean/sum of all means & errors.
// lists all these in .txt file propagatedIntegralMeans in Order of angles!


//Including root functionalities:
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TNtuple.h>
#include <TImage.h>
#include <TAttImage.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <numeric>
#include <tuple>
#include <map>

//andrea
#include "plastScintTimes.h"
#include "plastScintAmps.h"
#include "integralsSiPM.h"
#include "timeResDependencies.h"

// git test

using namespace std;

int main(int argc, char *argv[]) {
  string rootfileStr = argv[1]; // takes in format e.g. ../RootAnalysis/finishedRootfiles/18_cosmics_vb58_50k_2808_PARTS.root
  
  // string directory = "" + (string) argv[2];
  // find out run file name, assuming it ends in .root
  // std::cout << "2 " << directory << endl;

  string runNameRaw;
  // cout << rootfileStr.size()-6 << endl;
  for (int i = rootfileStr.size()-6; i >= 0; i--) {
    if (rootfileStr.at(i) == '/' || rootfileStr.at(i) == '\\') {
      if (i == rootfileStr.size()-6) {
        runNameRaw = "EGNAHC_eman_tluafed"; //"default_name_CHANGE";
      }
      break;
    }
    runNameRaw += rootfileStr.at(i);
  }

  string runName = "";
  for (int i=0; i<runNameRaw.size(); i++) {
    runName += runNameRaw.at(runNameRaw.size()-1-i);
  }
  // cout << "finished runName " << runName << endl;

  string runNumberString = "";
  for (int i=0; i<runName.size(); i++) {
    if (runName.at(i) == '_') break;
    runNumberString += runName.at(i);
  }

  // if not otherwise specified, files will be save to ../RootAnalysis/integralAnalysis/runName
  // HARDCODED PATH
  string saveFolder;
  if (argc < 3) { // || (string(argv[2])).empty() ) {
    saveFolder = "/mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/integralAnalysis/" + runName;
  } else {
    string saveFolder = argv[2]; // path specified when analysis is called (makeIntegralsPDF.sh)
  }

/*
  const Int_t nCh = 8;
  Int_t angles[nCh] = { 0,45,90,135,180,225,270,315 };
  Int_t channelOrder[nCh] = { 0,1,2,3,4,5,6,7 }; // XXX //0 deg --> channel 6, 45 deg --> channel 5, 90 deg --> channel 4, ... , 315 deg --> channel 7
*/
  int runNumber = stoi(runNumberString); //18; // DUMMY

  //Creating variables to read data from the TTree stored in the .file into our own TTree:
  int runPosition;
  float runEnergy;
  float photonCountPerEvent;
  int nEntries;
                      //eight channels for WOM in COSMICS
  /*
  //Style Settings:
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);
  gStyle->SetLineScalePS(1); */

  
  std::string location = saveFolder + "/integralMeans.txt";

  FILE* integralMeans = fopen(location.c_str(), "w"); // file to save unpropagated integral means

  // TLine* meanLineVec[nCh];                //array of lines for mean values
  


  // OPEN ROOTFILE
  printf("Analysing    file: %s\n", rootfileStr.c_str());
  TFile file(rootfileStr.c_str());
 

  if (file.IsZombie())
    {
        std::cout << "Problem with file " << rootfileStr << "; check if file path is correct!" << endl;
        exit(-1);
    }
  // std::cout << file.IsZombie() << endl;


  TTree* tree = new TTree;
  file.GetObject("T", tree);
  tree->GetEntry(1);

  //find out run parameters
  float rawAngleLimitLower, rawAngleLimitUpper;
  float rawTopPositionLimitLower, rawTopPositionLimitUpper;
  float rawBotPositionLimitLower, rawBotPositionLimitUpper;
  int centerChannel;

  tree->SetBranchAddress("angleIntTop", &rawAngleLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("angleIntBot", &rawAngleLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("posTopIntTop", &rawTopPositionLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("posTopIntBot", &rawTopPositionLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("posBotIntTop", &rawBotPositionLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("posBotIntBot", &rawBotPositionLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("CenterChannelPhiEw", &centerChannel);
  tree->GetEntry(1);


  float angleUpperLimit; // angle interval right (from normal) deg
  float angleLowerLimit; // angle interval left (from normal) deg
  float posTopLeftLimit; // left position limit upper plastic scint, in cm from 0
  float posTopRightLimit; // right position limit upper plastic scint, in cm from 0
  float posBotLeftLimit; // left position limit lower plastic scint, in cm from 0
  float posBotRightLimit; // right position limit lower plastic scint, in cm from 0
  
  if (rawAngleLimitUpper == 999) {
    angleUpperLimit = 90.0;
    angleLowerLimit = -90.0;
  } 
  else {
    angleUpperLimit = 90.0 - TMath::ATan(0.66 * 2 / (0.2998 * rawAngleLimitUpper)) / TMath::Pi() * 180.0;
    angleLowerLimit = -90.0 - TMath::ATan(0.66 * 2 / (0.2998 * rawAngleLimitLower)) / TMath::Pi() * 180.0;
  }

  if (rawTopPositionLimitLower == 999) {
    posTopLeftLimit = -100;
    posTopRightLimit = 100;
  }
  else {
    posTopLeftLimit = 0.2998 * 0.5 * 0.5 * rawTopPositionLimitLower * 100.0;
    posTopRightLimit = 0.2998 * 0.5 * 0.5 * rawTopPositionLimitUpper * 100.0;
  }

  if (rawBotPositionLimitLower == 999) {
    posBotLeftLimit = -100;
    posBotRightLimit = 100;
  }
  else {
    posBotLeftLimit = 0.2998 * 0.5 * 0.5 * rawBotPositionLimitLower * 100.0;
    posBotRightLimit = 0.2998 * 0.5 * 0.5 * rawBotPositionLimitUpper * 100.0;
  }

  float positionInfo[] = {(float) runNumber, posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit};

  // cout << " maybe raw angle limit upper " << TMath::ATan(0.66 * 2 / (0.2998 * rawAngleLimitUpper)) / TMath::Pi() * 180.0 << endl;

  cout << "Writing ..." << endl;

  //print Histograms of incidence times and amplitudes in Plastic Scintillators to PDF

  plasticScintTimes(tree, rootfileStr, saveFolder, positionInfo);
  cout << "- plasticScintTimes.pdf" << endl;
  plasticScintAmps(tree, rootfileStr, saveFolder, positionInfo);
  cout << "- plasticScintAmps.pdf" << endl;
  integralsSiPM(tree, rootfileStr, saveFolder, positionInfo);

  timeResDependencies(tree, rootfileStr, saveFolder, positionInfo);
  cout << "- timeResDependencies.pdf" << endl;

  return 0;
}
