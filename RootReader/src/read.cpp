//root
#include <TLine.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
#include <TText.h>
#include <TSpectrum.h>   // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h>      // root verbosity level
#include <TSystem.h>     // root verbosity level
#include <TLatex.h>      // root verbosity level

#include <sys/resource.h>
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
#include <string>
//specific
#include "geometry.h"
#include "analysis.h"
#include "misc.h"
#include "read.h"
#include <linux/limits.h>

// andrea
#include "meanAngleFuncs.h"

// andrea
int firstTrigger = 8; // first of 4 trigger channels. COSMICS
bool ANGLECUTS = false; //

bool POSITIONCUTS = true; // this automatically makes angle cut now!
bool INTEGRALCUT = false;
bool SCINTCUT = false;
bool FILTERWEIRD = true; // filter "weird" (may 2021) PMT signals by cutting PMT amp over 5mV

float lowAmpScint = -2000; // lowest amplitude in ch10 still counted
float highAmpScint = -60; // highest amplitude in ch10 still counted
float integralCut = 300; // 300 standard
float integralCutTop = 1500; // 1500.0 standard

float dTintervalTop = 1; // ANGLE cut upper limit (PMT) UNUSED
float dTintervalBot = -1; // angle cut lower limit (PMT) UNUSED

// pos,   top,            bottom
// pos0: (-1.75, 2.25), (-1.15, 2.85)
// pos1: (6.25; -),     (6.85, -)
// pos2: (-, -5.75),    (-, -5.15)
int position = 0; //0, 1, 2, 99 (none)


int angleMode = 0; // 

float diffTopIntervalTop; // POSITION cut upper limit (upper PMT)
float diffTopIntervalBot; // position cut lower limit (upper PMT)
float diffBotIntervalTop; // POSITION cut upper limit (lower PMT)
float diffBotIntervalBot; // position cut lower limit (lower PMT

int entriesChannelSum = 0;
float SiPMMaximumAverage;
int weirdEvents = 0;



Int_t angles[8] = { 0,45,90,135,180,225,270,315 };
Int_t channelOrder[8] = { 0,1,2,3,4,5,6,7 }; // XXX //e.g. 0 deg --> channel 6, 45 deg --> channel 5, 90 deg --> channel 4, ... , 315 deg --> channel 7

Float_t phi_ew[9];
float phiStd; // Standard deviation of phi_ew (no omission) for each event
Float_t threePhiEw[3];
int centerChannel = 0; // channel to center phi_ew around for no omissions
// end andrea


float SP = 0.3125; // ns per bin
float pe = 47.46;  //mV*ns

vector<float> calibrationCharges = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};      // dummy
vector<float> calibrationChargeErrors = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // dummy
string calibrationRunName = "dummy"; // "7_calib_vb58_tune8700_pcbd";
string dcIntegrationWindow = "";

vector<float> BL_const = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};                       // dummy
vector<float> integrationWindowsEntireSignal = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // dummy , left is always -20 measured to the max
vector<float> integrationWindowsPeakSignal = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};   // dummy , left is always -20 measured to the max
vector<float> correctionValues = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
;                                                                                                                                       // dummy , left is always -20 measured to the max
vector<float> correctionValueErrors = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // dummy , left is always -20 measured to the max

double coef = 2.5 / (4096 * 10);

// Run Parameter
Int_t runPosition = -999;
Float_t runEnergy = -999;
Int_t runAngle = -999;
Int_t runNumber = -999;
Int_t runChannelNumberWC = 12; //ANDREA
/*Declare & define the variables that are to be saved in the root-tree or that are used during the analysis.*/
Int_t EventNumber = -999;
Int_t LastEventNumber = -999;
unsigned long long int TDCsamIndex;
Float_t SamplingPeriod = -999;
Double_t EpochTime = -999;
Int_t Year = -999;
Int_t Month = -999;
Int_t Day = -999;
Int_t Hour = -999;
Int_t Minute = -999;
Int_t Second = -999;
Int_t Millisecond = -999;
Float_t trigT = -999; //t_trig = (t0+t1+t2+t3)/4
Float_t tSUMp = -999;
Float_t tSUMm = -999;
Int_t nCh = -1;
int womCount = 4;
int nActiveCh = -1;
int numberOfBinaryFiles = 0;
string workingDir;
char cwd[PATH_MAX];
#define btoa(x) ((x) ? "true" : "false")

//SETTINGS------------------------------------------------------------------------------------------
int defaultErrorLevel = kError;
bool print = true;
int headerSize = 328;
bool newerVersion = false;
bool switch_BL = false; // SWITCH dynamic <-> constant baseline // true = dyn, false = const
bool isDC = false;
bool automaticWindow = true;
bool useConstCalibValues = false; //IF the calibration values are correct, otherwise use dummies
float integralStart = 100;        //Testbeam: 100, 125 charge, Calib Nobember 2019: 165+25ns, Calib LED: 150  -> Very important: Integrationsstart and Time Window -> Determines the gain, slight changes result in big gain differences
float integralEnd = integralStart + 100;

int triggerChannel = 10; //starting from 0 -> Calib: 8?, Testbeam '18: 15, Important for timing tSipm,...
// cosmics: [10, 13] as coincidence trigger channels

int plotGrid = 3;

int maximalExtraPrintEvents = 50;
int printedExtraEvents = 0;
bool printExtraEvents = true;
//Event Skipping
bool skipThisEvent = false;
int skippedCount = 0;
// andrea
// event skipping cuts
bool skipThisEventInCut = false;
// andrea
bool weirdSkip;
//Skip events with bad baseline
bool allowBaselineEventSkipping = false;
int skipInChannel = 0;
//Skip events that exceed the WC maximum, does not include TriggerChannel
bool allowExceedingEventSkipping = true;
int exceedingThreshold = 600;
//Skip veto events -> if channel sees something-> Skip
bool allowVetoSkipping = false;
int vetoChannel = 9;
int vetoThreshold = 4; //abs Value -> compares with Amplitude

bool zoomedInWaves = false; //Zoom in the waves.pdf on the signal range

bool enableBaselineCorrection = true;
//Allow Force Printing individual events
bool allowForcePrintEvents = true;
bool forcePrintThisEvent = false;
int maximalForcePrintEvents = 50;
int forcePrintEvents = 0;

struct rusage r_usage;

/***
 *      ____  _____    _    ____    ____   ____ ____  ___ ____ _____ 
 *     |  _ \| ____|  / \  |  _ \  / ___| / ___|  _ \|_ _|  _ \_   _|
 *     | |_) |  _|   / _ \ | | | | \___ \| |   | |_) || || |_) || |  
 *     |  _ <| |___ / ___ \| |_| |  ___) | |___|  _ < | ||  __/ | |  
 *     |_| \_\_____/_/   \_\____/  |____/ \____|_| \_\___|_|    |_|  
 *      readscript                                                             
 */

void read(map<string, string> readParameters)
{

  /***
 *     __        __              ___ ___  ___  __   __  
 *    |__)  /\  |__)  /\   |\/| |__   |  |__  |__) /__` 
 *    |    /~~\ |  \ /~~\  |  | |___  |  |___ |  \ .__/ 
 *     paramaters                                                 
 */

  // andrea
  if (!ANGLECUTS) {
    dTintervalTop = 999;
    dTintervalBot = 999;
  }

   

  // assign PMT interval cuts according to chosen position and angle
  // float positions[4];
  // if (0 == angleMode) {
    float positionOptions[][4] = {{-1.75, 2.25, -1.15, 2.85},
                                    {6.25, 999, 6.85, 999},
                                    {-999, -5.75, -999, -5.15},
                                    {-999, 999, -999, 999}};
  // } 
  // else {
  //   float positions0[4] = {-999, 999, -999, 999};
  //   float positions1[4] = {-999, 999, -999, 999};
  //   float positions2[4] = {-999, 999, -999, 999};
  //   float positions99[4] = {-999, 999, -999, 999};
  // }

  float* positions;
  
  if (POSITIONCUTS) {
    switch(position) {
      case 0:
        positions = positionOptions[0];
        break;
      case 1:
        positions = positionOptions[1];
        break;
      case 2:
        positions = positionOptions[2];
        break;
      default:
        positions = positionOptions[3];
        break;
    }
  } else {
    positions = positionOptions[3];
  }
  diffTopIntervalBot = positions[0];
  diffTopIntervalTop = positions[1];
  diffBotIntervalBot = positions[2];
  diffBotIntervalTop = positions[3];
  
  // end andrea

  gErrorIgnoreLevel = defaultErrorLevel;
  std::string runName = readParameters["runName"];
  // cout << "runname in read" << runName << "\n" << endl;
  TString inFileList = readParameters["inFileList"];
  TString inDataFolder = readParameters["inDataFolder"];
  TString outFile = readParameters["outFile"];
  int headerSize = stoi(readParameters["headerSize"]);

  try
  {

    // std::cerr <<"#### runNumber:'" << readParameters["runNumber"] << "'\n";
    // std::cerr <<"#### pposition:'" << readParameters["runPosition"] << "'\n";
    // std::cerr <<"####ang:'" << readParameters["runAngle"] << "'\n";
    // std::cerr <<"#### runener:'" << readParameters["runEnergy"] << "'\n";
    // std::cerr <<"#### runNumberWC:'" << readParameters["runChannelNumberWC"] << "'\n";
    runNumber = stoi(readParameters["runNumber"]);
    runPosition = stoi(readParameters["runPosition"]);
    runAngle = stoi(readParameters["runAngle"]);
    runEnergy = stoi(readParameters["runEnergy"]);
    runChannelNumberWC = stoi(readParameters["runChannelNumberWC"]);
  }
  catch (const std::exception &e)
  {
    // std::cerr <<"Error at runNumber:" <<e.what() << '\n';
  }

  // HERE
  // cout << "RUNCHANNELWNUMBERWC " << runChannelNumberWC << endl;


  switch_BL = parseBoolean(readParameters["dynamicBL"]);
  isDC = parseBoolean(readParameters["isDC"]);
  useConstCalibValues = parseBoolean(readParameters["useConstCalibValues"]);
  automaticWindow = parseBoolean(readParameters["automaticWindow"]);
  string iWForceRun = readParameters["iWForceRun"];
  plotGrid = ceil(sqrt(runChannelNumberWC));

  if (getcwd(cwd, sizeof(cwd)) != NULL)
  {
    workingDir = cwd;
  }
  else
  {
    perror("getcwd() error");
    assert(0);
  }

  string charge_file = "/src/CalibrationValues.txt";
  string baseline_file = "/src/Baselines.txt";
  string integrationWindowPath = "/src/IntegrationWindows.txt";
  string correctionFactorPath = "/src/CorrectionValues.txt";

  string calib_path_charge = workingDir + charge_file;
  // cout << calib_path_charge << " CAlIBRATION PATH \n" << endl;
  string calib_path_bl = workingDir + baseline_file;
  string integrationWindowFile = workingDir + integrationWindowPath;
  string correctionValueFile = workingDir + correctionFactorPath;
  cout << "BL_path" << calib_path_bl << endl;

  if (useConstCalibValues)
  {
    pair<vector<float>, vector<float>> pairIW = readPair(calib_path_charge, calibrationRunName, 1, 0);
    calibrationCharges = pairIW.first;
    // cout << "calibration charges " << vectorToString(calibrationCharges) << "\n" << endl;
    calibrationChargeErrors = pairIW.second;
  }
  string iwSelection = runName;

  if (runName.find("cosmics") != std::string::npos)
  {
    allowVetoSkipping = true;
  }

  if (runName.find("dc") != std::string::npos)
  {
    isDC = true;
    if (automaticWindow)
    {
      iwSelection = dcIntegrationWindow;
    }
  }

  if (!automaticWindow)
  {
    if (iWForceRun.size() > 0)
    {
      iwSelection = iWForceRun;
    }
  }
  int iwRun = stoi(iwSelection.substr(0, iwSelection.find("_")));

  if (!switch_BL)
    BL_const = readVector(calib_path_bl, runName, 0);
  if (automaticWindow || iWForceRun.size() > 0)
  {
    pair<vector<float>, vector<float>> pairIW = readPair(integrationWindowFile, iwSelection, 0, 0);
    integrationWindowsPeakSignal = pairIW.first;
    integrationWindowsEntireSignal = pairIW.second;

    pair<vector<float>, vector<float>> pairCF = readPair(correctionValueFile, iwSelection, 1, 0);
    correctionValues = pairCF.first;
    correctionValueErrors = pairCF.second;
  }

  /***
 *    ___  __   ___  ___     __   ___ ___       __  
 *     |  |__) |__  |__     /__` |__   |  |  | |__) 
 *     |  |  \ |___ |___    .__/ |___  |  \__/ |    
 *       tree setup                                           
 */

  /*Create root-file and root-tree for data*/
  TFile *rootFile = new TFile(outFile, "RECREATE");
  if (rootFile->IsZombie())
  {
    if (numberOfBinaryFiles > 1)
    {
      cout << "PROBLEM with the initialization of the output ROOT ntuple "
           << outFile << ": check that the path is correct!!!"
           << endl;
    }
    exit(-1);
  }

  TTree *tree = new TTree("T", "USBWC Data Tree");
  TTree::SetBranchStyle(0);
  gStyle->SetLineScalePS(1); // export high resolution plots


/* // andrea
// set up additional  root-tree for data with applied cut on timing interval
  TTree *treeCut = new TTree("T", "USBWC Data Tree");
  TTree::SetBranchStyle(0);
  gStyle->SetLineScalePS(1); // export high resolution plots
 */


  Int_t ChannelNr[runChannelNumberWC];
  int WOMID[runChannelNumberWC];      //0=A, 1=B, 2=C, 3=D
  float chPE[runChannelNumberWC];     // single channel amplitude at sum signal
  float chPE_int[runChannelNumberWC]; // single channel integral
  std::vector<float> max(runChannelNumberWC, -999);
  std::vector<float> min(runChannelNumberWC, -999);
  Float_t t[runChannelNumberWC];

  Float_t tSiPM[(runChannelNumberWC - 1)]; //Minus Trigger
  float Integral[runChannelNumberWC];
  float IntegralErrorP[runChannelNumberWC]; //Integral with positive Error maxedout
  float IntegralErrorM[runChannelNumberWC];

  float IntegralDiff[runChannelNumberWC];
  float IntegralSum[runChannelNumberWC];
  float IntegralSumErrorP[runChannelNumberWC];
  float IntegralSumErrorM[runChannelNumberWC];
  
  //andrea
  Float_t IncidenceTime[4]; // for trig channels save time of signal incidence for coincidence cuts
  Float_t timeDifferenceTop; // calculate time differences between PMTs
  Float_t timeDifferenceBot;
  Float_t timeDifference; // time difference of top bottom time difference
  // these are old names, translates to the four trigger channels, in order
  Float_t Incidence10; // incidence time channel 10
  Float_t Incidence11; // incidence time channel 11
  Float_t Incidence12; // incidence time channel 12
  Float_t Incidence13; // incidence time channel 13

  Float_t timeMeanTop; // arithmetic mean
  Float_t timeMeanBot;
  Float_t timeResApprox; // (t1 + t2 - (t3 + t4))/4
  Float_t amplitudeMeanTop; // geometric mean
  Float_t amplitudeMeanBot;
  Float_t meanFlightTime;
  // these are old names, translates to the four trigger channels, in order
  Float_t Invert_Incidence10; // incidence time channel 10, calc by CFDInvert
  Float_t Invert_Incidence11; // incidence time channel 11, calc by CFDInvert
  Float_t Invert_Incidence12;
  Float_t Invert_Incidence13;
  Float_t inv_timeDifferenceTop;
  Float_t inv_timeDifferenceBot;
  Float_t inv_timeDifference;

  Float_t minCh10; // minimum of signal in channel 10
  Float_t minCh11;
  Float_t minCh12;
  Float_t minCh13;

  float integral_hist[8]; // save signal integrals for SiPM channels (COSMICs)

  // cartesian values for phiew calculations 
  float cartX;  // cartesian x value
  float cartY;  // cartesian y value
  float sumCartX; // sum of cart. x values
  float sumCartY; // sum of cart. y values
  float sigmaX; //std dev cart x values
  float sigmaY; //std dev cart y values
  float cartXarray[8]; // array of individual weighted x values
  float cartYarray[8]; // array of individual weighted y values
  // end andrea





  float Amplitude[runChannelNumberWC];
  float AmplitudeSum[runChannelNumberWC];
  float BL_output[4];                        //array used for output getBL-function
  Float_t BL_lower[runChannelNumberWC];      //store baseline for runChannelNumberWC channels for 0-75ns range
  Float_t BL_RMS_lower[runChannelNumberWC];  //store rms of baseline for runChannelNumberWC channels for 0-75ns range
  Float_t BL_Chi2_lower[runChannelNumberWC]; //store chi2/dof of baseline-fit for runChannelNumberWC channels for 0-75ns range
  Float_t BL_pValue_lower[runChannelNumberWC];
  Float_t BL_upper[runChannelNumberWC];      //store baseline for runChannelNumberWC channels for 220-320ns range
  Float_t BL_RMS_upper[runChannelNumberWC];  //store rms of baseline for runChannelNumberWC channels for 220-320ns range
  Float_t BL_Chi2_upper[runChannelNumberWC]; //store chi2/dof of baseline-fit for runChannelNumberWC channels for 220-320ns range
  Float_t BL_pValue_upper[runChannelNumberWC];
  Float_t BL_used[runChannelNumberWC];
  Float_t BL_Chi2_used[runChannelNumberWC];
  Float_t BL_pValue_used[runChannelNumberWC];
  float noiseLevel[runChannelNumberWC];
  int NumberOfBins;
  Int_t EventIDsamIndex[runChannelNumberWC];
  Int_t FirstCellToPlotsamIndex[runChannelNumberWC];
  Short_t amplValues[runChannelNumberWC][1024];

  TH1F hCh("hCh", "dummy;time [ns];voltage [mV]", 1024, -0.5 * SP, 1023.5 * SP);

  std::vector<TH1F *> hChSum;
  std::vector<TH1F> hChtemp;
  std::vector<TH1F *> hChShift;



  Float_t amplitudeChannelSumWOM[womCount];
  Float_t chargeChannelSumWOM[womCount];
  Float_t chargeChannelSumWOMErrorP[womCount];
  Float_t chargeChannelSumWOMErrorM[womCount];


  std::vector<TH1F *> histChannelSumWOM;
  int binNumber = 1024; //Default: 1024, change with caution

  if (print)
  {

    for (int i = 0; i < runChannelNumberWC; i++)
    {
      TString name("");
      name.Form("hChSum_%d", i);
      TH1F *h = new TH1F("h", ";time [ns];voltage [mV]", binNumber, -0.5 * SP, 1023.5 * SP);
      h->SetName(name);
      hChSum.push_back(h);
    }

    for (int i = 0; i < runChannelNumberWC; i++)
    {
      TString name("");
      name.Form("hChShift_%d", i);
      TH1F *h = new TH1F("h", ";time [ns];voltage [mV]", binNumber, -0.5 * SP, 1023.5 * SP);
      h->SetName(name);
      hChShift.push_back(h);
    }

    for (int i = 0; i < womCount; i++)
    {
      TString name("");
      name.Form("histChannelSumWOM%d", i);
      TH1F *h = new TH1F("h", ";time [ns];voltage [mV]", binNumber, -0.5 * SP, 1023.5 * SP);
      h->SetName(name);
      histChannelSumWOM.push_back(h);
    }

    //andrea
    // for (int i = 0; i < 8; i++)
    // {
    //   TString name("");
    //   name.Form("histChannelIntegral%d", i);
    //   TH1F *h = new TH1F("h", ";integral [ns*mV];counts", binNumber, -0.5 * SP, 1023.5 * SP);
    //   h->SetName(name);
    //   histChannelIntegral.push_back(h);
    // }
  }

  for (int i = 0; i < runChannelNumberWC; i++)
  {
    TString name("");
    name.Form("hChtemp_%d", i);
    TH1F h("h", ";time [ns];voltage [mV]", binNumber, -0.5 * SP, 1023.5 * SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }

  TString plotSaveFolder = outFile;
  plotSaveFolder.ReplaceAll((runName + ".root"), "");
  plotSaveFolder.ReplaceAll(("out.root"), "");

  TCanvas cWaves("cWaves", "cWaves", 1000, 1000);
  cWaves.Divide(plotGrid, plotGrid);
  TCanvas womCanvas("womCanvas", "womCanvas", 1000, 1000);
  womCanvas.Divide(2, 2);
  TCanvas cChSum("cChSum", "cChSum", 1500, 900);
  cChSum.Divide(plotGrid, plotGrid);

  // andrea
  TCanvas cIntegralsWOM("cIntegralsWOM", "cIntegralsWOM", 1000, 1000);
  cIntegralsWOM.Divide(2, 4);

  /*Create branches in the root-tree for the data.*/
  tree->Branch("EventNumber", &EventNumber, "EventNumber/I");
  tree->Branch("SamplingPeriod", &SamplingPeriod, "SamplingPeriod/F");
  tree->Branch("EpochTime", &EpochTime, "EpochTime/D");
  tree->Branch("Year", &Year, "Year/I");
  tree->Branch("Month", &Month, "Month/I");
  tree->Branch("Day", &Day, "Day/I");
  tree->Branch("Hour", &Hour, "Hour/I");
  tree->Branch("Minute", &Minute, "Minute/I");
  tree->Branch("Second", &Second, "Second_/I");
  tree->Branch("Millisecond", &Millisecond, "Millisecond/I");
  tree->Branch("trigT", &trigT, "trigT/F");
  tree->Branch("tSUMp", &tSUMp, "tSUMp/F");
  tree->Branch("tSUMm", &tSUMm, "tSUMm/F");

  //RUN PARAMETER
  tree->Branch("runPosition", &runPosition, "runPosition/I");
  tree->Branch("runEnergy", &runEnergy, "runEnergy/F");
  tree->Branch("runAngle", &runAngle, "runAngle/I");
  tree->Branch("runNumber", &runNumber, "runNumber/I");
  tree->Branch("runChannelNumberWC", &runChannelNumberWC, "runChannelNumberWC/I");
  tree->Branch("integrationWindowRun", &iwRun, "integrationWindowRun/I");

  // CHANNEL INFO (but everything that is nCH-dependend below)
  tree->Branch("nCh", &nCh, "nCh/I");
  tree->Branch("WOMID", WOMID, "WOMID[nCh]/I");
  tree->Branch("ch", ChannelNr, "ch[nCh]/I");
  // AMPLITUDE
  tree->Branch("Amplitude", Amplitude, "Amplitude[nCh]/F"); // calibrated
  //tree->Branch("amp_inRange", amp_inRange.data(), "amp_inRange[nCh]/F"); // calibrated
  tree->Branch("max", max.data(), "max[nCh]/F");
  tree->Branch("min", min.data(), "min[nCh]/F");
  // INTEGRAL
  //tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  //tree->Branch("Integral_inRange", Integral_inRange, "Integral_inRange[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("IntegralErrorP", IntegralErrorP, "IntegralErrorP[nCh]/F");
  tree->Branch("IntegralErrorM", IntegralErrorM, "IntegralErrorM[nCh]/F");

  // tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F"); // calibrated
  tree->Branch("IntegralDifference", IntegralDiff, "IntegralDifference[nCh]/F");

  // TIMING
  tree->Branch("t", t, "t[nCh]/F");
  tree->Branch("tSiPM", tSiPM, "tSiPM[nCh]/F");

  //andrea
  tree->Branch("Time_diff_top", &timeDifferenceTop, "Time_diff_top/F");
  tree->Branch("Time_diff_bot", &timeDifferenceBot, "Time_diff_bot/F");
  //tree->Branch("Trigger_time", IncidenceTime, "Trigger_time[4]/F");
  tree->Branch("Time_diff_all", &timeDifference, "Time_diff_all/F");
  tree->Branch("Incidence_time_ch10", &Incidence10, "Incidence_time_ch10/F");
  tree->Branch("Incidence_time_ch11", &Incidence11, "Incidence_time_ch11/F");
  tree->Branch("Incidence_time_ch12", &Incidence12, "Incidence_time_ch12/F");
  tree->Branch("Incidence_time_ch13", &Incidence13, "Incidence_time_ch13/F");

  tree->Branch("minimum_ch10", &minCh10, "minimum_ch10/F");
  tree->Branch("minimum_ch11", &minCh11, "minimum_ch11/F");
  tree->Branch("minimum_ch12", &minCh12, "minimum_ch12/F");
  tree->Branch("minimum_ch13", &minCh13, "minimum_ch13/F");
  tree->Branch("timeMeanTop", &timeMeanTop, "timeMeanTop/F");
  tree->Branch("timeMeanBot", &timeMeanBot, "timeMeanBot/F");
  tree->Branch("timeResApprox", &timeResApprox, "timeResApprox/F");
  tree->Branch("amplitudeMeanTop", &amplitudeMeanTop, "amplitudeMeanTop/F");
  tree->Branch("amplitudeMeanBot", &amplitudeMeanBot, "amplitudeMeanBot/F");
  tree->Branch("meanFlightTime", &meanFlightTime, "meanFlightTime/F");
 
  for (int i=0; i<8; i++) {
    tree->Branch(Form("integral_hist_%d", i), &(integral_hist[i]), Form("integral_hist_%d/F", i));
  }
  for (int i=0; i<8; i++) {
    tree->Branch(Form("Phi_ew_omit_ch%d", i), &(phi_ew[i]), Form("Phi_ew_omit_ch%d/F", i));
  }
  tree->Branch("Phi_ew_all_ch", &(phi_ew[8]), "Phi_ew_all_ch/F");
  tree->Branch("Std_Phi_ew_all", &phiStd, "Std_Phi_ew_all/F");
  tree->Branch("Phi_ew_shifted", threePhiEw, "Phi_ew_shifted[3]/F");
  tree->Branch("CenterChannelPhiEw", &centerChannel, "CenterChannelPhiEw/I");
  // tree->Branch("integral_hist_0", &integral_hist[0], "integral_hist_0/F");
  
  tree->Branch("angleIntTop", &dTintervalTop, "angleIntTop/F");
  tree->Branch("angleIntBot", &dTintervalBot, "angleIntBot/F");
  tree->Branch("posTopIntTop", &diffTopIntervalTop, "posIntTop/F");
  tree->Branch("posTopIntBot", &diffTopIntervalBot, "posIntBot/F");
  tree->Branch("posBotIntTop", &diffBotIntervalTop, "posIntTop/F");
  tree->Branch("posBotIntBot", &diffBotIntervalBot, "posIntBot/F");
  tree->Branch("integralCut", &integralCut, "integralCut/F");
  tree->Branch("integralCutTop", &integralCutTop, "integralCutTop/F");
  tree->Branch("lowAmpScint", &lowAmpScint, "lowAmpScint/F");
  tree->Branch("highAmpScint", &highAmpScint, "highAmpScint/F");

  tree->Branch("sumCartX", &sumCartX, "sumCartX/F");
  tree->Branch("sumCartY", &sumCartY, "sumCartY/F");
  tree->Branch("sigmaSumCartY", &sigmaY, "sigmaSumCartY/F");
  tree->Branch("sigmaSumCartX", &sigmaX, "sigmaSumCartX/F");
  for (int i=0; i<8; i++) {
    tree->Branch("cartXi", &cartXarray[i], "cartXi/F");
    tree->Branch("cartYi", &cartYarray[i], "cartYi/F");
  }

  int multiPeaks = 0;
  tree->Branch("multiPeaks", &multiPeaks, "multiPeaks/I");
  // end andrea



  
  // BASELINE
  tree->Branch("BL_lower", BL_lower, "BL_lower[nCh]/F");
  tree->Branch("BL_RMS_lower", BL_RMS_lower, "BL_RMS_lower[nCh]/F");
  tree->Branch("BL_Chi2_lower", BL_Chi2_lower, "BL_Chi2_lower[nCh]/F");
  tree->Branch("BL_pValue_lower", BL_pValue_lower, "BL_pValue_lower[nCh]/F");
  tree->Branch("BL_upper", BL_upper, "BL_upper[nCh]/F");
  tree->Branch("BL_RMS_upper", BL_RMS_upper, "BL_RMS_upper[nCh]/F");
  tree->Branch("BL_Chi2_upper", BL_Chi2_upper, "BL_Chi2_upper[nCh]/F");
  tree->Branch("BL_pValue_upper", BL_pValue_upper, "BL_pValue_upper[nCh]/F");
  tree->Branch("BL_used", BL_used, "BL_used[nCh]/F");
  tree->Branch("BL_Chi2_used", BL_Chi2_used, "BL_Chi2_used[nCh]/F");
  tree->Branch("BL_pValue_used", BL_pValue_used, "BL_pValue_used[nCh]/F");
  // CALIBRATED SUM
  tree->Branch("chargeChannelSumWOM", chargeChannelSumWOM, "chargeChannelSumWOM[4]/F");
  tree->Branch("chargeChannelSumWOMErrorP", chargeChannelSumWOMErrorP, "chargeChannelSumWOMErrorP[4]/F");
  tree->Branch("chargeChannelSumWOMErrorM", chargeChannelSumWOMErrorM, "chargeChannelSumWOMErrorM[4]/F");

  tree->Branch("amplitudeChannelSumWOM", amplitudeChannelSumWOM, "amplitudeChannelSumWOM[4]/F");
  tree->Branch("EventIDsamIndex", EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex", FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  tree->Branch("nSkipped", skippedCount, "nSkipped/I");




  /***
 *     __   ___       __          __      __  ___       __  ___ 
 *    |__) |__   /\  |  \ | |\ | / _`    /__`  |   /\  |__)  |  
 *    |  \ |___ /~~\ |__/ | | \| \__>    .__/  |  /~~\ |  \  |  
 *         reading start                                                     
 */

  /*Start reading the raw data from .bin files.*/
  int nitem = 1;
  ifstream inList;
  TString fileName;
  inList.open(inFileList);
  assert(inList.is_open());

  //Get Binary File Count
  ifstream countStream(inFileList);
  numberOfBinaryFiles = count(std::istreambuf_iterator<char>(countStream),
                              std::istreambuf_iterator<char>(), '\n');
  countStream.close();

  //Get First Line File Name-> for printing
  string tempFileName;
  ifstream tempName(inFileList);
  if (tempName.good())
  {
    getline(tempName, tempFileName);
  }
  tempName.close();

  //cout << "Number of Binary Files: " << numberOfBinaryFiles << endl;
  bool printParameterOverview = false;
  if (tempFileName.substr(tempFileName.find_last_of(".") + 1) == "bin")
  {
    printParameterOverview = true;
  }

  if (numberOfBinaryFiles > 1 || printParameterOverview)
  {
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
    cout << ":::::::::::::::::::RUN PARAMETER:::::::::::::::::::::::::::::::::::::::" << endl;

    // Iterate through all elements in std::map
    map<string, string>::iterator it = readParameters.begin();
    while (it != readParameters.end())
    {
      cout << it->first << " :: " << it->second << endl;
      it++;
    }
    cout << ":::::::::::::::::::CALIBRATION:::::::::::::::::::::::::::::::::::::::::" << endl;
    cout << "isDC: " << isDC << endl;
    cout << "allowVetoSkipping: " << allowVetoSkipping << endl;
    cout << "Baselines(Constant): " << vectorToString(BL_const) << endl;
    cout << "Charge Calibration: " << vectorToString(calibrationCharges) << endl;
    cout << "Charge CalibrationErr: " << vectorToString(calibrationChargeErrors) << endl;

    cout << "IntegrationWindowPeak: " << vectorToString(integrationWindowsPeakSignal) << endl;
    cout << "IntegrationWindowEntire: " << vectorToString(integrationWindowsEntireSignal) << endl;
    cout << "Is DarkCount: " << btoa(isDC) << " Dynamic Baseline: " << btoa(switch_BL) << " Is Calibrated: " << btoa(useConstCalibValues) << endl;

    cout << "CorrectionValues: " << vectorToString(correctionValues) << endl;
    cout << "CorrectionValuesErrors: " << vectorToString(correctionValueErrors) << endl;
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
  }

  int fileCounter = 0;
  int currentPrint = -1;
  int totalWeirdEvents = 0;
  int totalEvents = 0;

  /***
 *     ___         ___          __   __   __  
 *    |__  | |    |__     |    /  \ /  \ |__) 
 *    |    | |___ |___    |___ \__/ \__/ |    
 *               file loop                             
 */


  while (inList >> fileName)
  {
    int eventsInFile = 0;
    int weirdInFile = 0;
    fileName = inDataFolder + fileName;
    // cout << fileName << endl;
    FILE *pFILE = fopen(fileName.Data(), "rb");
    if (pFILE == NULL)
    {
      fputs("File error", stderr);
      assert(0);
    }
    fseek(pFILE, 0, SEEK_END);
    float totFileSizeByte = ftell(pFILE);
    rewind(pFILE);

    /***
 *          ___       __   ___  __                              __     __  
 *    |__| |__   /\  |  \ |__  |__)     /\  |\ |  /\  |    \ / /__` | /__` 
 *    |  | |___ /~~\ |__/ |___ |  \    /~~\ | \| /~~\ |___  |  .__/ | .__/ 
 * header analysis                                                                        
 */

    if (headerSize == -1)
    {
      char versionBuffer[700];
      fread(versionBuffer, 1, 700, pFILE);

      string softwareVersion;
      softwareVersion += versionBuffer[44]; //i = 44   2
      softwareVersion += versionBuffer[45]; //i = 45   .
      softwareVersion += versionBuffer[46]; //i = 46   9
      softwareVersion += versionBuffer[47]; //i = 47   .
      softwareVersion += versionBuffer[48]; //i = 48   1
      softwareVersion += versionBuffer[49]; //i = 49   0

      std::size_t found = string(versionBuffer).find("Correction"); //Determine the position of "INL Correction: 0" inside the header
      for (int o = 1; o <= 500; o++)
      { //test 500 more bytes after this position, if they are NOT PRINTABLE-> this is were the content starts
        if (isprint(versionBuffer[found + o]) == 0)
        {
          headerSize = found + o;
          break;
        }
      }
      if (headerSize > 328)
      {
        newerVersion = true;
      }
      if (numberOfBinaryFiles > 1 || printParameterOverview)
        cout << "VERSION: " << softwareVersion << " HEADERSIZE: " << headerSize << "  Newer Version: " << newerVersion << endl;
    }

    if (headerSize > 328)
    {
      newerVersion = true;
    }

    fseek(pFILE, 0, SEEK_SET);
    char header[headerSize];
    nitem = fread(header, 1, headerSize, pFILE);
    //cout << "Header:\n"
    //<< header << endl;

    char *word;
    word = strtok(header, " \n");
    while (word != NULL)
    {
      if (strcmp("ACQUIRED:", word) == 0)
      {
        word = strtok(NULL, " \n");
        nActiveCh = atoi(word);
        break;
      }
      word = strtok(NULL, " \n");
    }

    //Not sure why reading this byte, //if (nActiveCh > 9 || newerVersion), needs to be read in all versions 2.8.14+
    char dummy;
    nitem = fread(&dummy, 1, 1, pFILE);

    /***
 *    $$$$$$$$\ $$\    $$\ $$$$$$$$\ $$\   $$\ $$$$$$$$\       $$\       $$$$$$\   $$$$$$\  $$$$$$$\  
 *    $$  _____|$$ |   $$ |$$  _____|$$$\  $$ |\__$$  __|      $$ |     $$  __$$\ $$  __$$\ $$  __$$\ 
 *    $$ |      $$ |   $$ |$$ |      $$$$\ $$ |   $$ |         $$ |     $$ /  $$ |$$ /  $$ |$$ |  $$ |
 *    $$$$$\    \$$\  $$  |$$$$$\    $$ $$\$$ |   $$ |         $$ |     $$ |  $$ |$$ |  $$ |$$$$$$$  |
 *    $$  __|    \$$\$$  / $$  __|   $$ \$$$$ |   $$ |         $$ |     $$ |  $$ |$$ |  $$ |$$  ____/ 
 *    $$ |        \$$$  /  $$ |      $$ |\$$$ |   $$ |         $$ |     $$ |  $$ |$$ |  $$ |$$ |      
 *    $$$$$$$$\    \$  /   $$$$$$$$\ $$ | \$$ |   $$ |         $$$$$$$$\ $$$$$$  | $$$$$$  |$$ |      
 *    \________|    \_/    \________|\__|  \__|   \__|         \________|\______/  \______/ \__|  
                                                                                                                                                                                                                                                                                                
 */ // event loop

    int whileCounter = 0;
    while (nitem > 0)
    { //event loop

      eventsInFile += 1;

      skipThisEvent = false;
      skipThisEventInCut = false;
      weirdSkip = false;  // andrea
      forcePrintThisEvent = false;
      std::vector<TObject *> eventTrash;
      whileCounter++;
      nitem = fread(&EventNumber, sizeof(int), 1, pFILE);
      nitem = fread(&EpochTime, sizeof(double), 1, pFILE);
      nitem = fread(&Year, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Month, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Day, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Hour, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Minute, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Second, sizeof(unsigned int), 1, pFILE);
      nitem = fread(&Millisecond, sizeof(unsigned int), 1, pFILE);
      if (newerVersion)
        nitem = fread(&TDCsamIndex, 1, 8, pFILE); //New in 2.9.10

      nitem = fread(&nCh, sizeof(unsigned int), 1, pFILE); // since V2.8.14 the number of stored channels is written for each event

      

      if (EventNumber % 100 == 0)
      {
        printf("Percentage: %lf, EventNr: %d, nCh: %d+   \n", ftell(pFILE) / totFileSizeByte, EventNumber, nCh);
        // HERE
        // cout << "NCH " << nCh << endl;
      }

      float MeasuredBaseline[runChannelNumberWC];
      float AmplitudeValue[runChannelNumberWC];
      float ComputedCharge[runChannelNumberWC];
      float RiseTimeInstant[runChannelNumberWC];
      float FallTimeInstant[runChannelNumberWC];
      float RawTriggerRate[runChannelNumberWC];
      float floatR = -1;

      /*Loop over individual channels. For each event the data from every channel is 
      processed and analysed one by one in order*/

      for (int i = 0; i < nCh; i++)
      {
        // start old version
        // use for 17_cosmics_vb58_50k_2808!!!
        /*
        nitem = fread(&ChannelNr[i], sizeof(int), 1, pFILE);
        nitem = fread(&EventIDsamIndex[i], sizeof(int), 1, pFILE);
        nitem = fread(&FirstCellToPlotsamIndex[i], sizeof(int), 1, pFILE);
        nitem = fread(&floatR, 1, 4, pFILE);
        MeasuredBaseline[i] = floatR;
        nitem = fread(&floatR, 1, 4, pFILE);
        AmplitudeValue[i] = floatR;

        nitem = fread(&floatR, 1, 4, pFILE);
        ComputedCharge[i] = floatR;
        nitem = fread(&floatR, 1, 4, pFILE);
        RiseTimeInstant[i] = floatR;
        nitem = fread(&floatR, 1, 4, pFILE);
        FallTimeInstant[i] = floatR;
        nitem = fread(&floatR, 1, 4, pFILE);
        RawTriggerRate[i] = floatR;
        ChannelNr[i] = i;
        */ 
        // end old version

        //new version
        // use for measurements 04/2021 onward
        nitem = fread(&ChannelNr[i], sizeof(int), 1, pFILE);
        nitem = fread(&EventIDsamIndex[i], sizeof(int), 1, pFILE);
        nitem = fread(&FirstCellToPlotsamIndex[i], sizeof(int), 1, pFILE);
        // end new version

        /***
 *          __              __  
 *    |  | /  \  |\/|    | |  \ 
 *    |/\| \__/  |  |    | |__/ 
 *             wom id                 
 */

        if (i < 8)
        {
          WOMID[i] = 3;
        }
        else if (i < 16)
        {
          WOMID[i] = 2;
        }
        else if (i < 24)
        {
          WOMID[i] = 0;
        }
        else if (i < 31)
        {
          WOMID[i] = 1;
        }
        else
        {
          WOMID[i] = -1;
        }
        if (i == triggerChannel)
          WOMID[i] = -1;
        TString title("");
        title.Form("channel: %d, event: %d", i, EventNumber);
        hCh.Reset();
        hCh.SetTitle(title);

        /*
        __ Waveform Histogram _______________________________________________
        Writing the signal amplitude values into the root-histogram hCh.
        */

        for (int j = 0; j < 1024; j++)
        {
          nitem = fread(&amplValues[i][j], sizeof(short), 1, pFILE);
          hCh.SetBinContent(j + 1, (amplValues[i][j] * coef * 1000));
        }

        /*The error of each value in each bin is set to 0.5 mV.*/
        for (int j = 1; j <= hCh.GetXaxis()->GetNbins(); j++)
        {
          hCh.SetBinError(j, 0.5);
        }

        /*Analysis if the event/signal starts.*/
        max[i] = hCh.GetMaximum();
        min[i] = hCh.GetMinimum();

        /*Saving the histogram of that event into a temporary histogram hChtemp. These histograms are available outside of the channel-loop. If analysis using the signals/events of multiple channels needs to be done, this can be accomplished by using hChtemp after the channel-loop.*/
        hChtemp.at(i) = hCh;

        //Exceeding Events Skipping

        if (i != triggerChannel && max[i] > exceedingThreshold)
        {
          // forcePrintThisEvent = true;

          if (allowExceedingEventSkipping && !skipThisEvent)
          {
            skipThisEvent = true;
          }
          else
          {
            //  skipThisEvent = false;
          }
        }

        /***
 *     __        __   ___              ___ 
 *    |__)  /\  /__` |__  |    | |\ | |__  
 *    |__) /~~\ .__/ |___ |___ | | \| |___ 
 *         baseline                                
 */

        // BL_fit(&hChtemp.at(i), BL_output, 0.0, 75.0);
        BL_fit(&hChtemp.at(i), BL_output, 0.0, 30.0);
        BL_lower[i] = BL_output[0];
        BL_RMS_lower[i] = BL_output[1];
        BL_Chi2_lower[i] = BL_output[2];
        BL_pValue_lower[i] = BL_output[3];
        // BL_fit(&hChtemp.at(i), BL_output, 220.0, 320.0);

        BL_fit(&hChtemp.at(i), BL_output, 290.0, 320.0);
        BL_upper[i] = BL_output[0];
        BL_RMS_upper[i] = BL_output[1];
        BL_Chi2_upper[i] = BL_output[2];
        BL_pValue_upper[i] = BL_output[3];

        // determine "best" baseline
        if (BL_Chi2_upper[i] <= BL_Chi2_lower[i] && false) //always use lower -> upper might not be representing the signal correctly since APO
        {
          BL_used[i] = BL_upper[i];
          BL_Chi2_used[i] = BL_Chi2_upper[i];
          BL_pValue_used[i] = BL_pValue_upper[i];
        }
        else
        {
          BL_used[i] = BL_lower[i];
          BL_Chi2_used[i] = BL_Chi2_lower[i];
          BL_pValue_used[i] = BL_pValue_lower[i];
        }

        if (i == skipInChannel && allowBaselineEventSkipping)
        {
          if (BL_Chi2_used[i] > 1)
          {
            skipThisEvent = true;
          }
          else
          {
            // skipThisEvent = false;
          }
        }

        // SWITCH dynamic <-> constant baseline
        float BL_shift;
        if (switch_BL)
        {
          BL_shift = BL_used[i];
        }
        else
        {
          BL_shift = BL_const[i];
        }
        if (!enableBaselineCorrection)
          BL_shift = 0;

        //Baseline Correction

        TF1 *f_const = new TF1("f_const", "pol0", 0, 320);
        f_const->SetParameter(0, BL_shift);
        hCh.Add(f_const, -1); //this = this + c1*h1

        /***
 *    ___                 __      __     __         __  
 *     |  |  |\/| | |\ | / _`    /__` | |__)  |\/| /__` 
 *     |  |  |  | | | \| \__>    .__/ | |     |  | .__/ 
 *                               timing sipms                       
 */
        /* Deactivated for COSMICS from 04/2021 since four trigger channels now
        if (i == triggerChannel)
        {                                //trigger
          t[i] = CFDNegative(&hCh, 0.5); //Negative Trigger
        }
        */

        // andrea
		    if ((i==firstTrigger || i==firstTrigger+1 || i==firstTrigger+2 || i==firstTrigger+3)){ 
          float threshold = 0.5;
          int offset = 120;
          if (EventNumber % 10000 == 0) {
            cout << "Threshold for CFD in timing: " << threshold << "Offset: " << offset << endl;
          }

          // andrea
          Float_t inc = CFDNegativeCustom(&hChtemp.at(i), threshold, offset); // search half peak from peak - offset up

          // cout << "channel number inc time" << i << endl;

          Float_t invert_inc = CFDNegativeInvert(&hChtemp.at(i), threshold); // search half peak from peak down
			    IncidenceTime[i-8] = inc; // want to record signal time from PMTs
          
          
          
          Float_t signalMinimum = hChtemp.at(i).GetMinimum();
          Float_t signalMaximum = hChtemp.at(i).GetMaximum();
          
          if (signalMinimum > (-20.0)) {
            // skipThisEvent = true;
            // continue;
          } 
          if (SCINTCUT && (signalMinimum < lowAmpScint || signalMinimum > highAmpScint) ) { //&& (i==firstTrigger || i==firstTrigger+1)) {
              skipThisEvent = true;
          }

          if (FILTERWEIRD && (signalMaximum > 5)) {
            skipThisEvent = true;
            weirdSkip = true;
            
          } else {
            if (FILTERWEIRD && (EventNumber%50 == 0) ) forcePrintThisEvent = true; // ATTENTION
          }  

          // if (!skipThisEvent) {
          //   cout<< "signalMaximum" << signalMaximum << endl;
          // }

          
          // cout << "minimum signal " << signalMinimum << endl;
          // bla bla blaaa
          if (firstTrigger == i) {
            Incidence10 = inc;
            Invert_Incidence10 = invert_inc;
            minCh10 = signalMinimum;
            
            if (inc > 350 || inc < 0) { 
              //forcePrintThisEvent = true;
              skipThisEvent = true;
              cout << "Inc time weird " << whileCounter << endl;
              cout << inc << minCh10 << endl;
            }
          }   
          if (firstTrigger+1 == i) {
            Incidence11 = inc;
            Invert_Incidence11 = invert_inc;
            minCh11 = signalMinimum;
          }
          if (firstTrigger+2 == i) {
            Incidence12 = inc;
            Invert_Incidence12 = invert_inc;
            minCh12 = signalMinimum;
          }
          if (firstTrigger+3 == i) {
            Incidence13 = inc;
            Invert_Incidence13 = invert_inc;
            minCh13 = signalMinimum;           
          }
		    } //end andrea
        
        else
        { //SiPMs
          t[i] = CFDInRange(&hCh, 0.35, integralStart, integralEnd);

          if (t[i] < 95)
          {
            t[i] = CFDinvertInRange(&hCh, 0.35, integralStart, integralEnd);
          }
        }

		//andrea
        if (firstTrigger+1 == i) {
          timeDifferenceTop = Incidence11 - Incidence10;
          timeMeanTop = 0.5*(Incidence11 + Incidence10);
          amplitudeMeanTop = sqrtf(fabs(minCh10 * minCh11));

          inv_timeDifferenceTop = Invert_Incidence11 - Invert_Incidence10;
        }	

        if (firstTrigger+3 == i) {
          timeDifferenceBot = Incidence13 - Incidence12;
          inv_timeDifferenceBot = Invert_Incidence13 - Invert_Incidence12;
          inv_timeDifference = inv_timeDifferenceTop - inv_timeDifferenceBot;
          timeDifference = timeDifferenceTop - timeDifferenceBot;
          
          timeMeanBot = 0.5*(Incidence12 + Incidence13);
          amplitudeMeanBot = sqrtf(fabs(minCh12 * minCh13));
          timeResApprox = 0.25*(Incidence11 + Incidence10 - (Incidence12 + Incidence13));
          meanFlightTime = 0.5*(Incidence11 + Incidence10 - (Incidence12 + Incidence13));
          float absTimeDifference = fabs(timeDifference);

        // cout << " uh oh " << absTimeDifference << endl;
        if (15 < absTimeDifference) {
          // skipThisEvent = true;
          //cout << "skip this event? " << skipThisEvent << endl;
          //cout << " greater 15! " << endl;
          //skipThisEventInCut = true;
        }

        // cut on TimeDifference -> angle:
        // only events in Interval get counted
        if (ANGLECUTS && !(timeDifference <= dTintervalTop && timeDifference >= dTintervalBot)) {
          skipThisEvent = true;
        }

        // cut on timeDiffTop -> Position and angle
        // only events in interval get counted:
        if (POSITIONCUTS && !(timeDifferenceTop <= diffTopIntervalTop && timeDifferenceTop >= diffTopIntervalBot)) {
          skipThisEvent = true; // Dt top in interval positions[2,3]
        }
        if (POSITIONCUTS && !(timeDifferenceBot <= diffBotIntervalTop && timeDifferenceBot >= diffBotIntervalBot)) {
          skipThisEvent = true; // Dt bot in interval positions[0,1]
        }

        
        //end andrea

          
        }

        /***
 *           ___  ___  __   __               |                __         ___       __   ___ 
 *    | |\ |  |  |__  / _` |__)  /\  |       |     /\   |\/| |__) |    |  |  |  | |  \ |__  
 *    | | \|  |  |___ \__> |  \ /~~\ |___    |    /~~\  |  | |    |___ |  |  \__/ |__/ |___ 
 *    integral                                       |                                              
 */
        hCh.SetStats(0);

        float t_amp = t_max_inRange(&hCh, integralStart, integralEnd);
        int integrationLeftOffset = 20;
        if (runName.find("calib") != std::string::npos)
        {
          //correctionValues[i] = 1;      //exclude for calib runs
          //  correctionValueErrors[i] = 0; //exclude for calib runs
          integrationLeftOffset = 10;
        }
        bool dynamicDCWindow = true; //DO NOT CHANGE, Unless you know what you do
        if (isDC)
        {
          if (!dynamicDCWindow)
            t_amp = 120;
          correctionValues[i] = 1;      //exclude for calib runs
          correctionValueErrors[i] = 0; //exclude for calib runs
        }

        float integralStartShifted = t_amp - integrationLeftOffset;
        float integralEndShifted = t_amp + integrationWindowsPeakSignal[i];
        float integralEndShiftedAll = t_amp + integrationWindowsEntireSignal[i];

        int shiftedIndex = i + 0 * 8; //calib values are ordered D C A B; if one wants to measure data taken with SIPM A -> Shift index by 2*8
        if (shiftedIndex == 32)
          shiftedIndex = 31;

        Amplitude[i] = AmplitudeHist(&hCh, integralStartShifted, integralEndShifted, 0);

        float calibrationError = calibrationChargeErrors[shiftedIndex];
        float correctionValueError = correctionValueErrors[shiftedIndex];

        float effectivFactor = correctionValues[shiftedIndex] / calibrationCharges.at(shiftedIndex);
        float effectivFactorError = sqrt(pow((correctionValueError / calibrationCharges.at(shiftedIndex)), 2) + pow((correctionValues[shiftedIndex] * calibrationError / pow(calibrationCharges.at(shiftedIndex), 2)), 2));

        Integral[i] = IntegralHist(&hCh, integralStartShifted, integralEndShifted, 0) ; //* effectivFactor;
       
        // andrea
        // Integral stuff only for the 8 SiPM channels
        if (i < 8) {
          integral_hist[i] = Integral[i];
        // cout << Integral[i] << endl;

          if (INTEGRALCUT && i<8 && ((Integral[i] < integralCut) || (Integral[i] > integralCutTop))) { // (Integral[i] < integralCut) ||
            skipThisEvent = true;
          }

          float SiPMMaximum = hChtemp.at(i).GetMaximum();
          // cout << SiPMMaximum << endl;
          if (0==i) {
            SiPMMaximumAverage = 0;
            bool manyPeaks = false;
          }
         
          SiPMMaximumAverage += SiPMMaximum;

          // ATTENTION
          // if (0==i && !skipThisEvent){
          //   TSpectrum *s = new TSpectrum(2*5);
          //   Int_t nfound = s->Search(&hChtemp.at(i),2,"",0.80);
          //   // cout << "nfound " << nfound << endl;
          //   if (nfound != 1) multiPeaks += 1;
          //   }

          if (FILTERWEIRD && i == 7 && (SiPMMaximumAverage/8.0 < 10)) {

            // skipThisEvent = true;
          }
        }

        // end andrea


        IntegralErrorP[i] = IntegralHist(&hCh, integralStartShifted, integralEndShifted, 0) * (effectivFactor + effectivFactorError);
        IntegralErrorM[i] = IntegralHist(&hCh, integralStartShifted, integralEndShifted, 0) * (effectivFactor - effectivFactorError);
        IntegralDiff[i] = IntegralDifference(&hCh, integralStartShifted, integralEndShifted, integralEndShiftedAll, Amplitude[shiftedIndex], 0);

        float ampForVeto = 0.0;
        if (allowVetoSkipping && i == vetoChannel)
        {
          ampForVeto = AmplitudeHistAlternative(&hCh, 0, 320, 0); //Search Everywhere
          if (ampForVeto > vetoThreshold)
          {
            skipThisEvent = true;
            //forcePrintEvent=true;
          } 
        }

        // if (!skipThisEvent && (maximalForcePrintEvents >= forcePrintEvents)) {
        //   forcePrintThisEvent = true;
        //   forcePrintEvents++;
        // }

        // cout<<"DD: "<<IntegralDiff[i]<<"  "<<i<<endl;
        // skipThisEvent=true;

        if (WOMID[i] >= 0)
          histChannelSumWOM[WOMID[i]]->Add(&hCh);

        /***
 *     __   __         ___         __                     ___  __  
 *    |__) |__) | |\ |  |  | |\ | / _`    |  |  /\  \  / |__  /__` 
 *    |    |  \ | | \|  |  | | \| \__>    |/\| /~~\  \/  |___ .__/ 
 *                                                                 
 */

        if (print)
        {
          //  if(!skipThisEvent || (skipThisEvent && (forcePrintThisEvent || (printedExtraEvents < maximalExtraPrintEvents)))){
          //  if ((forcePrintThisEvent && (forcePrintEvents < maximalForcePrintEvents)) || ((currentPrint != fileCounter) || (printedExtraEvents < maximalExtraPrintEvents)))

          if ( (allowForcePrintEvents || ((currentPrint != fileCounter) || (printedExtraEvents < maximalExtraPrintEvents)))) // !skipThisEvent && 
          {

            cWaves.cd(i + 1);
            if (zoomedInWaves)
              hCh.GetXaxis()->SetRange((integralStartShifted - 50) / SP, (integralEndShiftedAll + 50) / SP);

            if (i == 8)
              hCh.GetXaxis()->SetRange(0, 150.0 / SP);

            hCh.DrawCopy("HIST"); //No error bars pls
            hCh.GetXaxis()->SetRange((t[i] - 20) / SP, (t[i] + 30) / SP);
            int max_bin = hCh.GetMaximumBin();
            int lower_bin = max_bin - 20.0 / SP;
            int upper_bin = max_bin + 30.0 / SP;
            // double x = h->GetXaxis()->GetBinCenter(binmax);
            float max_time = hCh.GetXaxis()->GetBinCenter(max_bin);
            float lower_time = hCh.GetXaxis()->GetBinCenter(lower_bin);
            float upper_time = hCh.GetXaxis()->GetBinCenter(upper_bin);
            hCh.GetXaxis()->SetRange(0, 1024);
            TLine *ln4 = new TLine(0, BL_lower[i], 30, BL_lower[i]);
            TLine *ln5 = new TLine(290, BL_upper[i], 320, BL_upper[i]);
            TText *text = new TText(.5, .5, Form("%f | %f", BL_lower[i], BL_upper[i]));
            ln4->SetLineColor(2);
            ln5->SetLineColor(2);
            ln4->SetLineWidth(3);
            ln5->SetLineWidth(3);

            TLine *baselineUsed = new TLine(0, 0, 320, 0); //The entire waveform is already moved by this BL_shift amount -> this will be at 0
            baselineUsed->SetLineColor(3);
            baselineUsed->SetLineWidth(2);

            float minY = hCh.GetMinimum();
            float maxY = hCh.GetMaximum();
            TLine *leftInterval;
            TLine *rightInterval;
            TLine *endInterval;

            leftInterval = new TLine(integralStartShifted, minY, integralStartShifted, maxY);
            rightInterval = new TLine(integralEndShifted, minY, integralEndShifted, maxY);
            endInterval = new TLine(integralEndShiftedAll, minY, integralEndShiftedAll, maxY);

            leftInterval->SetLineWidth(3);
            rightInterval->SetLineWidth(3);
            endInterval->SetLineWidth(3);

            leftInterval->SetLineColor(2);
            rightInterval->SetLineColor(2);
            endInterval->SetLineColor(2);

            baselineUsed->Draw("same");
            //ln4->Draw("same"); //always use the same baseline -> already shown in baselineused
            //ln5->Draw("same"); //Upper baseline isnt used, dont draw it

            leftInterval->Draw("same");
            rightInterval->Draw("same");
            endInterval->Draw("same");
            TLegend *h_leg = new TLegend(0.50, 0.65, 0.90, 0.90);
            bool lightMode = false;
            if (allowVetoSkipping && i == vetoChannel)
            {
              lightMode = true;
              h_leg->AddEntry((TObject *)0, Form("Amplitude VETO: %1.2f", ampForVeto), "");
              float sumD = Amplitude[0] + Amplitude[1] + Amplitude[2] + Amplitude[3] + Amplitude[4] + Amplitude[5] + Amplitude[6] + Amplitude[7];
              h_leg->AddEntry((TObject *)0, Form("SumIntD: %1.2f", sumD), "");
            }

            h_leg->SetTextSize(0.020);
            if (allowForcePrintEvents)
              h_leg->AddEntry((TObject *)0, Form("ForcePrinted: %d, Count: %d, Skip: %d", forcePrintThisEvent, forcePrintEvents, skipThisEvent), "");
            h_leg->AddEntry((TObject *)0, Form("Integral: %1.2f (+%1.2f/-%1.2f)", Integral[i], IntegralErrorP[i], IntegralErrorM[i]), "");
            h_leg->AddEntry((TObject *)0, Form("Amplitude: %1.2f", Amplitude[i]), "");
            if (!lightMode)
            {
              h_leg->AddEntry((TObject *)0, Form("Calibration Value: %1.2f +- %1.2f", calibrationCharges.at(shiftedIndex), calibrationChargeErrors.at(shiftedIndex)), "");
              h_leg->AddEntry((TObject *)0, Form("CF: %1.2f +- %1.2f", correctionValues[i], correctionValueErrors[i]), "");
              h_leg->AddEntry(baselineUsed, Form("BL: %1.2f (Histogram shifted)", BL_shift), "l");
              h_leg->AddEntry((TObject *)0, Form("Percentage: %1.2f", IntegralDiff[i]), "");
              h_leg->AddEntry((TObject *)0, Form("Eff. Factor: %1.2f+-%1.2f", effectivFactor, effectivFactorError), "");
              h_leg->AddEntry(leftInterval, Form("Window: Left: %d, Right %1.2f, %1.2f", integrationLeftOffset, integrationWindowsPeakSignal[i], integrationWindowsEntireSignal[i]), "l");
              h_leg->AddEntry(leftInterval, Form("Window Size: %1.0f", abs(integralStartShifted - integralEndShifted)), "l");
            }

            // h_leg->Clear();
            /*h_leg->AddEntry(ln4, Form("baseline start: %1.2f mV",  BL_lower[i]), "l");
            h_leg->AddEntry((TObject *)0, Form("reduced chi^{2}: %1.2f",  BL_Chi2_lower[i]), "");
           h_leg->AddEntry((TObject *)0, "", "");

           h_leg->AddEntry(ln5, Form("baseline end: %1.2f mV",  BL_upper[i]), "l");
          h_leg->AddEntry((TObject *)0, Form("reduced chi^{2}: %1.2f",  BL_Chi2_upper[i]), "");*/
            // h_leg->AddEntry((TObject *)0, Form("position: %d",  runPosition), "");
            //h_leg->AddEntry((TObject *)0, Form("energy: %1.1f GeV",  runEnergy*0.1), "");
            //h_leg->AddEntry((TObject *)0, Form("box rotation: %d deg",  runAngle), "");
            //h_leg->AddEntry((TObject *)0, "WOM: D", "");
            //  h_leg->AddEntry((TObject *)0, Form("f_W: %1.2f", IntegralDiff[i]), "");

            h_leg->Draw();

            //text->Draw("same");
          }
          hCh.GetXaxis()->SetRange(1, 30 / SP);
          noiseLevel[i] = hCh.GetMaximum() - hCh.GetMinimum();
          hCh.GetXaxis()->SetRange(1, 1024);
          // End of loop over inividual channels
        }

        if (IntegralDiff[i] > -88 && !skipThisEvent) {
          entriesChannelSum += 1;
          // cout << "added in channel sum" << endl;
          // cout << entriesChannelSum << endl;
          hChSum.at(i)->Add(&hCh, 1); //Dont sum empty waveforms into your sum histogram
        }
      } // end of channel loop

      // andrea
      // ATTENTION HERE
      
      Double_t pi = TMath::Pi();

      

      for (int i=0; i<9; i++) { // do loop 9 times, omit every channel once, omit no channel once
        sumCartX = 0;
        sumCartY = 0;
        
        float individualPhi[8];

        for (int k=0; k<8; k++) { // sum weighted cart values from channels
          if (k != i) { // omits channel i. when i=8, no channel omitted
            int kk = channelOrder[k]; // take corresponding angle value for each channel
            cartX = TMath::Cos(angles[kk] * (pi/180.0)) * Integral[k];
            cartY = TMath::Sin(angles[kk] * (pi/180.0)) * Integral[k];           
            sumCartY += cartY;
            sumCartX += cartX;
          }        
        }

        if (8==i) {
          for (int j=0; j<8; j++) {
            int kk = channelOrder[j]; // take corresponding angle value for each channel
            cartX = TMath::Cos(angles[kk] * (pi/180.0)) * Integral[j];
            cartY = TMath::Sin(angles[kk] * (pi/180.0)) * Integral[j];
            cartXarray[j] = cartX;
            cartYarray[j] = cartY;
            sumCartY += 1/8.0 * cartY;
            sumCartX += 1/8.0 * cartX;
          }
          float sigmaX_raw;
          float sigmaY_raw;
          for (int j=0; j<8; j++) {
            sigmaX_raw += (sumCartX - cartXarray[j])*(sumCartX - cartXarray[j]);
            sigmaX_raw += (sumCartY - cartYarray[j])*(sumCartY - cartYarray[j]);
          }
          sigmaX = TMath::Sqrt(sigmaX_raw * 1/8.0);
          sigmaY = TMath::Sqrt(sigmaY_raw * 1/8.0);
        }

        // now to convert sumCartY and sumCartX back to polar coordinates, minding ATan periodicity
        // if (sumCartX > 0 && sumCartY >= 0) {
        //   phi_ew[i] = (TMath::ATan(sumCartY / sumCartX) * 180.0 / pi);
        // } 
        // else if (sumCartX < 0) {
        //   phi_ew[i] = (TMath::ATan(sumCartY / sumCartX)  + pi) * 180.0 / pi;
        // } 
        // else if (sumCartX > 0 && sumCartY < 0) {
        //   phi_ew[i] = (TMath::ATan(sumCartY / sumCartX)  + 2*pi) * 180.0 / pi;
        // } 
        // else if (sumCartX = 0 && sumCartY > 0) {
        //   phi_ew[i] = 0.5 * 180.0;
        // }
        // else if (sumCartX = 0 && sumCartY < 0) {
        //   phi_ew[i] = -1.5* 180.0;
        // }

        // now to convert sumCartY and sumCartX back to polar coordinates, minding ATan periodicity
        phi_ew[i] = cartesianToPolar(sumCartX, sumCartY);

        
        if (8 == i) {
          // cout << " xes " << endl;
          // printArray(cartXarray, 8);
          // cout << " ys " << endl;
          // printArray(cartYarray, 8);
          // cout << "convert " << endl;

          // also convert indivudual weighted angles to cartesian again
          // this is equivalent to making another angle list that is ordered according to which angle is assigned to which detector
          cartesianToPolar(8, cartXarray, cartYarray, individualPhi);

          // printArray(individualPhi, 8);

          // calculate the standard deviation of phi_ew for each event, only when no channel was omitted:
          float sigma_raw = 0;
          float difference = 0;
          for (int m=0; m<8; m++) {
            difference = phi_ew[i] - individualPhi[m];
            if (difference > 180) difference -= 180;
            else if (difference < -180) difference += 180;
            sigma_raw += difference * difference; // XXX
            // cout << sigma_raw << endl;
          }
          phiStd = sqrt(1/8.0 * sigma_raw);
        }
        

        // translate angles from [0, 360) range to (-180, 180] range centered around omitted channel (channel k for no omissions):
        // COSMICS specific! like everything else too basically
        int k;

        if (i > 3 && 8 != i) k = (i + 4) % (4);
        else if (8 != i) k = i + 4;
        else k = centerChannel; // no channels omitted


        phi_ew[i] = translateAngle(phi_ew[i], angles[k]);
        threePhiEw[0] = phi_ew[i] - 360.0;
        threePhiEw[1] = phi_ew[i];
        threePhiEw[2] = phi_ew[i] + 360.0;
      }



      
     // ATTENTION HERE
      //end andrea

      

      /***
 *          __            __              __  
 *    |  | /  \  |\/|    /__` |  |  |\/| /__` 
 *    |/\| \__/  |  |    .__/ \__/  |  | .__/ 
 *                                            
 */

      for (int i = 0; i < 4; i++)
      {
        womCanvas.cd(i + 1);
        histChannelSumWOM[i]->DrawCopy();

        if (i == 3)
        {
          //WOM D
          IntegralSum[i] = Integral[0] + Integral[1] + Integral[2] + Integral[3] + Integral[4] + Integral[5] + Integral[6] + Integral[7];
          IntegralSumErrorP[i] = IntegralErrorP[0] + IntegralErrorP[1] + IntegralErrorP[2] + IntegralErrorP[3] + IntegralErrorP[4] + IntegralErrorP[5] + IntegralErrorP[6] + IntegralErrorP[7];
          IntegralSumErrorM[i] = IntegralErrorM[0] + IntegralErrorM[1] + IntegralErrorM[2] + IntegralErrorM[3] + IntegralErrorM[4] + IntegralErrorM[5] + IntegralErrorM[6] + IntegralErrorM[7];

          AmplitudeSum[i] = Amplitude[0] + Amplitude[1] + Amplitude[2] + Amplitude[3] + Amplitude[4] + Amplitude[5] + Amplitude[6] + Amplitude[7];

          if (IntegralSum[3] < 6)
          {
            //fore if(allowForcePrintEvents)
            //  forcePrintThisEvent = true;
          }
        }
        else if (i == 2)
        {
          //WOM C
          IntegralSum[i] = Integral[8] + Integral[9] + Integral[10] + Integral[11] + Integral[12] + Integral[13] + Integral[14] + Integral[15];
          IntegralSumErrorP[i] = IntegralErrorP[8] + IntegralErrorP[9] + IntegralErrorP[10] + IntegralErrorP[11] + IntegralErrorP[12] + IntegralErrorP[13] + IntegralErrorP[14] + IntegralErrorP[15];
          IntegralSumErrorM[i] = IntegralErrorM[8] + IntegralErrorM[9] + IntegralErrorM[10] + IntegralErrorM[11] + IntegralErrorM[12] + IntegralErrorM[13] + IntegralErrorM[14] + IntegralErrorM[15];

          AmplitudeSum[i] = Amplitude[8] + Amplitude[9] + Amplitude[10] + Amplitude[11] + Amplitude[12] + Amplitude[13] + Amplitude[14] + Amplitude[15];
        }
        else if (i == 1)
        {
          //WOM B
          IntegralSum[i] = (8.0 / 7.0) * (Integral[24] + Integral[25] + Integral[26] + Integral[27] + Integral[28] + Integral[29] + Integral[30]);                                                 //1 Missing -> Trigger
          IntegralSumErrorP[i] = (8.0 / 7.0) * (IntegralErrorP[24] + IntegralErrorP[25] + IntegralErrorP[26] + IntegralErrorP[27] + IntegralErrorP[28] + IntegralErrorP[29] + IntegralErrorP[30]); //1 Missing -> Trigger
          IntegralSumErrorM[i] = (8.0 / 7.0) * (IntegralErrorM[24] + IntegralErrorM[25] + IntegralErrorM[26] + IntegralErrorM[27] + IntegralErrorM[28] + IntegralErrorM[29] + IntegralErrorM[30]); //1 Missing -> Trigger

          AmplitudeSum[i] = (8.0 / 7.0) * (Amplitude[24] + Amplitude[25] + Amplitude[26] + Amplitude[27] + Amplitude[28] + Amplitude[29] + Amplitude[30]); //1 Missing -> Trigger
        }
        else if (i == 0)
        {
          //WOM A
          IntegralSum[i] = Integral[16] + Integral[17] + Integral[18] + Integral[19] + Integral[20] + Integral[21] + Integral[22] + Integral[23];
          IntegralSumErrorP[i] = IntegralErrorP[16] + IntegralErrorP[17] + IntegralErrorP[18] + IntegralErrorP[19] + IntegralErrorP[20] + IntegralErrorP[21] + IntegralErrorP[22] + IntegralErrorP[23];
          IntegralSumErrorM[i] = IntegralErrorM[16] + IntegralErrorM[17] + IntegralErrorM[18] + IntegralErrorM[19] + IntegralErrorM[20] + IntegralErrorM[21] + IntegralErrorM[22] + IntegralErrorM[23];

          AmplitudeSum[i] = Amplitude[16] + Amplitude[17] + Amplitude[18] + Amplitude[19] + Amplitude[20] + Amplitude[21] + Amplitude[22] + Amplitude[23];
        }

        chargeChannelSumWOM[i] = IntegralSum[i];
        chargeChannelSumWOMErrorP[i] = IntegralSumErrorP[i];
        chargeChannelSumWOMErrorM[i] = IntegralSumErrorM[i];

        amplitudeChannelSumWOM[i] = AmplitudeSum[i];
      }

      if (skipThisEvent)
      {
        forcePrintThisEvent = false;
      }

      if (forcePrintThisEvent)
      {
        forcePrintEvents++;
      }

      /***
 *    ___                 __     ___  __     __   __   ___  __  
 *     |  |  |\/| | |\ | / _`     |  |__) | / _` / _` |__  |__) 
 *     |  |  |  | | | \| \__>     |  |  \ | \__> \__> |___ |  \ 
 *                                                              
 */

      trigT = t[triggerChannel];
      for (int i = 0; i < runChannelNumberWC; i++)
      {
        if (i != triggerChannel)
          tSiPM[i] = t[i] - trigT;
      }

      /***
 *     __       ___  __       ___ 
 *    /  \ |  |  |  |__) |  |  |  
 *    \__/ \__/  |  |    \__/  |  
 *                                
 */
      


      if (print)
      {

        // if (!skipThisEvent || (skipThisEvent && (forcePrintThisEvent || (printedExtraEvents < maximalExtraPrintEvents))))
        if (!skipThisEvent && (printedExtraEvents < maximalExtraPrintEvents)) // ANDREA
        {
          if ((forcePrintThisEvent && (forcePrintEvents < maximalForcePrintEvents)) || ((currentPrint != fileCounter) || (printedExtraEvents < maximalExtraPrintEvents)))
          {
            if (printExtraEvents)
              printedExtraEvents++;
            currentPrint = fileCounter;

            if (fileCounter == 0 && ((fileCounter != (numberOfBinaryFiles - 1)) || forcePrintThisEvent))
            {
              cWaves.Print((TString)(plotSaveFolder + "/waveforms.pdf("), "pdf");
              //cWaves.Print((TString)(plotSaveFolder + "/waveforms2.pdf("), "pdf");

              womCanvas.Print((TString)(plotSaveFolder + "/waveforms_womSum.pdf("), "pdf");
              // cout << "print in if filecounter" << endl;
              
            } else
            {
              // HERE
              // cWaves.Print((TString)(plotSaveFolder + "/waveforms.pdf"), "pdf");
              womCanvas.Print((TString)(plotSaveFolder + "/waveforms_womSum.pdf"), "pdf");
              // cout << "print in else" << endl;
            }
          }
        }
      }
      cWaves.Clear("D");

      /*Writing the data for that event to the tree.*/
      if (!skipThisEvent)
      {
        tree->Fill();
      }
      else
      {
        skippedCount = skippedCount + 1;
      }

      if (weirdSkip) weirdInFile += 1;
    } // end of event loop


    auto nevent = tree->GetEntries();

    cout << "Events in Tree:  " << nevent << " Skipped:  " << skippedCount << endl;
    fclose(pFILE);
    fileCounter++;
    // andrea
    totalEvents += eventsInFile;
    totalWeirdEvents += weirdInFile;
    // end andrea

  } // end of file loop

  

  if (print)
  {
    /*Clearing objects and saving files.*/
    inList.close();

    if ((numberOfBinaryFiles != 1 || forcePrintEvents > 0))
    {
      cWaves.Print((TString)(plotSaveFolder + "/waveforms.pdf)"), "pdf");


      womCanvas.Print((TString)(plotSaveFolder + "/waveforms_womSum.pdf)"), "pdf");

    }
    cWaves.Clear();
    womCanvas.Clear();

    for (int i = 0; i < runChannelNumberWC; i++)
    {
      cChSum.cd(i + 1);
      // hChSum.at(i)->GetXaxis()->SetRange(100,250);
      hChSum.at(i)->GetXaxis()->SetLabelSize(0.04);
      hChSum.at(i)->GetYaxis()->SetLabelSize(0.04);

      hChSum.at(i)->SetStats(0);

      hChSum.at(i)->Draw("HIST");
      //   hChSum.at(i)->SetFillColorAlpha(4, 0.8);
    }
    cChSum.Print((TString)(plotSaveFolder + "/waveforms_chSum.pdf"), "pdf");


  }

  gErrorIgnoreLevel = kWarning;

  rootFile = tree->GetCurrentFile();
  

  //int skipcount= &skippedCount;
  //gDirectory->WriteObject(skipcount,"skippedCount");
  rootFile->Write();
  rootFile->Close();

  // andrea
  FILE* cut_log = fopen("/mnt/d/Work_SHK_Bachelor/RootReader/runlogs/cut_log.txt", "a");

  FILE* weirdEventCounter = fopen("/mnt/d/Work_SHK_Bachelor/RootReader/runlogs/weirdEvents.txt", "a");

  fprintf(cut_log, "\n#Run log of timing cuts \n");
  fprintf(cut_log, "#Runname: \n");
  fprintf(cut_log, "%s", runName.c_str());
  fprintf(cut_log, "\nPosition cuts: %s\nAngle cuts (DEPR): %s\nIntegral cuts: %s\nPMT amp cuts: %s", POSITIONCUTS ? "true" : "false", ANGLECUTS ? "true" : "false", INTEGRALCUT ? "true" : "false", SCINTCUT ? "true" : "false");
  fprintf(cut_log, "\n#POS top min\tPOS top max\t POS bot min\t bot max\tANGLE interval left\tANGLE interval right\tmin. N_pe\tmax. N_pre\tmin PMT amp\tmax PMT amp\n");
  fprintf(cut_log, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", diffTopIntervalBot, diffTopIntervalTop, diffBotIntervalBot, diffBotIntervalTop, dTintervalBot, dTintervalTop, integralCut, integralCutTop, lowAmpScint, highAmpScint);
  fclose(cut_log);

  fprintf(weirdEventCounter, "%s\t%d\t%d\t%f\n", runName.c_str(), totalWeirdEvents, totalEvents, totalWeirdEvents*1.0/totalEvents);
  fclose(weirdEventCounter);
/*   
  TFile *rootFileCut = new TFile("cuts.root", "RECREATE");
  if (rootFileCut->IsZombie())
  {
    if (numberOfBinaryFiles > 1)
    {
      cout << "PROBLEM with the initialization of the output ROOT ntuple "
           << outFile << ": check that the path is correct!!!"
           << endl;
    }
    exit(-1);
  }

  rootFileCut = treeCut->GetCurrentFile();
  rootFileCut->Write();
  rootFileCut->Close();>*/
} 