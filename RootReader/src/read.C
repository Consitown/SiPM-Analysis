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
//specific
#include "geometry.h"
#include "analysis.h"
#include "misc.h"
#include "read.h"
#include <linux/limits.h>

float SP = 0.3125; // ns per bin
float pe = 47.46;  //mV*ns

vector<float> calibrationCharges = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};      // dummy
vector<float> calibrationChargeErrors = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // dummy
string calibrationRunName = ""; //7_calib_vb58_tune8700_pcbd
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
Int_t runChannelNumberWC = 32; //maximum
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

int triggerChannel = 8; //starting from 0 -> Calib: 8?, Testbeam '18: 15, Important for timing tSipm,...
int plotGrid = 3;

int maximalExtraPrintEvents = 0;
int printedExtraEvents = 0;
bool printExtraEvents = false;
//Event Skipping
bool skipThisEvent = false;
int skippedCount = 0;
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
bool allowForcePrintEvents = false;
bool forcePrintThisEvent = false;
int maximalForcePrintEvents = 2;
int forcePrintEvents = 0;

struct rusage r_usage;

/***
 *      ____  _____    _    ____    ____   ____ ____  ___ ____ _____ 
 *     |  _ \| ____|  / \  |  _ \  / ___| / ___|  _ \|_ _|  _ \_   _|
 *     | |_) |  _|   / _ \ | | | | \___ \| |   | |_) || || |_) || |  
 *     |  _ <| |___ / ___ \| |_| |  ___) | |___|  _ < | ||  __/ | |  
 *     |_| \_\_____/_/   \_\____/  |____/ \____|_| \_\___|_|    |_|  
 *                                                                   
 */

void read(map<string, string> readParameters)
{

  /***
 *     __        __              ___ ___  ___  __   __  
 *    |__)  /\  |__)  /\   |\/| |__   |  |__  |__) /__` 
 *    |    /~~\ |  \ /~~\  |  | |___  |  |___ |  \ .__/ 
 *                                                      
 */

  gErrorIgnoreLevel = defaultErrorLevel;
  string runName = readParameters["runName"];
  TString inFileList = readParameters["inFileList"];
  TString inDataFolder = readParameters["inDataFolder"];
  TString outFile = readParameters["outFile"];
  int headerSize = stoi(readParameters["headerSize"]);

  try
  {
    runNumber = stoi(readParameters["runNumber"]);
    runPosition = stoi(readParameters["runPosition"]);
    runAngle = stoi(readParameters["runAngle"]);
    runEnergy = stoi(readParameters["runEnergy"]);
    runChannelNumberWC = stoi(readParameters["runChannelNumberWC"]);
  }
  catch (const std::exception &e)
  {
    //  std::cerr <<"Error at runNumber:" <<e.what() << '\n';
  }

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
  string calib_path_bl = workingDir + baseline_file;
  string integrationWindowFile = workingDir + integrationWindowPath;
  string correctionValueFile = workingDir + correctionFactorPath;

  if (useConstCalibValues)
  {
    pair<vector<float>, vector<float>> pairIW = readPair(calib_path_charge, calibrationRunName, 1, 0);
    calibrationCharges = pairIW.first;
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
 *                                                  
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
 *                                                              
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

  /***
 *     ___         ___          __   __   __  
 *    |__  | |    |__     |    /  \ /  \ |__) 
 *    |    | |___ |___    |___ \__/ \__/ |    
 *                                            
 */

  while (inList >> fileName)
  {

    fileName = inDataFolder + fileName;
    cout << fileName << endl;
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
 *                                                                         
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
 */

    int whileCounter = 0;
    while (nitem > 0)
    { //event loop

      skipThisEvent = false;
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

        /***
 *          __              __  
 *    |  | /  \  |\/|    | |  \ 
 *    |/\| \__/  |  |    | |__/ 
 *                              
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
 *                                         
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
 *                                                      
 */

        if (i == triggerChannel)
        {                                //trigger
          t[i] = CFDNegative(&hCh, 0.5); //Negative Trigger
        }
        else
        { //SiPMs
          t[i] = CFDInRange(&hCh, 0.35, integralStart, integralEnd);

          if (t[i] < 95)
          {
            t[i] = CFDinvertInRange(&hCh, 0.35, integralStart, integralEnd);
          }
        }

        /***
 *           ___  ___  __   __               |                __         ___       __   ___ 
 *    | |\ |  |  |__  / _` |__)  /\  |       |     /\   |\/| |__) |    |  |  |  | |  \ |__  
 *    | | \|  |  |___ \__> |  \ /~~\ |___    |    /~~\  |  | |    |___ |  |  \__/ |__/ |___ 
 *                                           |                                              
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

        Integral[i] = IntegralHist(&hCh, integralStartShifted, integralEndShifted, 0) * effectivFactor;

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

          if (allowForcePrintEvents || ((currentPrint != fileCounter) || (printedExtraEvents < maximalExtraPrintEvents)))
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

        if (IntegralDiff[i] > -88 && !skipThisEvent)
          hChSum.at(i)->Add(&hCh, 1); //Dont sum empty waveforms into your sum histogram
      }

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

        if (!skipThisEvent || (skipThisEvent && (forcePrintThisEvent || (printedExtraEvents < maximalExtraPrintEvents))))
        {
          if ((forcePrintThisEvent && (forcePrintEvents < maximalForcePrintEvents)) || ((currentPrint != fileCounter) || (printedExtraEvents < maximalExtraPrintEvents)))
          {
            if (printExtraEvents)
              printedExtraEvents++;
            currentPrint = fileCounter;

            if (fileCounter == 0 && ((fileCounter != (numberOfBinaryFiles - 1)) || forcePrintThisEvent))
            {
              cWaves.Print((TString)(plotSaveFolder + "/waveforms.pdf("), "pdf");
              womCanvas.Print((TString)(plotSaveFolder + "/waveforms_womSum.pdf("), "pdf");
            }
            else
            {

              cWaves.Print((TString)(plotSaveFolder + "/waveforms.pdf"), "pdf");
              womCanvas.Print((TString)(plotSaveFolder + "/waveforms_womSum.pdf"), "pdf");
            }
          }
        }
      }

      /*Writing the data for that event to the tree.*/
      if (!skipThisEvent)
      {
        tree->Fill();
      }
      else
      {
        skippedCount = skippedCount + 1;
      }
    }
    auto nevent = tree->GetEntries();

    cout << "Events in Tree:  " << nevent << " Skipped:  " << skippedCount << endl;
    fclose(pFILE);
    fileCounter++;
  }

  if (print)
  {
    /*Clearing objects and saving files.*/
    inList.close();

    if (numberOfBinaryFiles != 1 || forcePrintEvents > 0)
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
}