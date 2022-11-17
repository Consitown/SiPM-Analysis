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
//#include <TStyle.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
//#include <stdlib.h>
//#include <string>
//#include <iomanip>

//specific
#include "geometry.h"
#include "analysis.h"
#include "readFull.h"
#define btoa(x) ((x) ? "true" : "false")

float SP = 0.3125;
float pe = 47.46; //mV*ns
//vector<float> pe_SiPM = {32.14, 39.33, 34.20, 30.79, 34.09, 29.99, 30.69, 29.95}; //a,b,c,d,e,f,g,h  -  Gain-Baseline from fit
vector<float> pe_SiPM = {42.01, 34.67, 34.28, 33.84, 37.55, 34.68, 33.81, 38.84}; //sorted by Wavecatcher-Channel
vector<float> SiPM_shift = {2.679, 2.532, 3.594, 3.855, 3.354, 3.886, 3.865, 4.754};
// vector<float> calib_amp_AB = {10.0024,9.24254,9.08902,10.0149,9.95047,9.55901,10.1483,10.4179,10.0141,9.92513,10.4975,10.422,10.1208,10.1884,10.1682,1};
vector<float> calib_amp_AB = {6.748, 6.16313, 6.07082, 6.68036, 6.65783, 6.37541, 6.7711, 6.85418, 6.68469, 6.58283, 6.98329, 6.97906, 6.76493, 6.75924, 6.78279, 1};
vector<float> calib_amp_AB_max = {9.91652, 8.86927, 8.88147, 9.57771, 9.58071, 9.14965, 9.53239, 10.1344, 9.62728, 9.62879, 10.0288, 10.3354, 9.75948, 9.53048, 9.68774, 1};

/*
__ CALIBRATION FACTORS of indiv. arrays _______________________________
*/
vector<float> calib_amp_8SiPM_17 = {1, 3.90, 3.35, 3.25, 3.25, 3.60, 3.25, 3.25, 3.55};
vector<float> calib_amp_8SiPM_1 = {1, 6.748, 6.16313, 6.14954, 6.07082, 6.68036, 6.65783, 6.37541, 6.7711};
vector<float> calib_amp_8SiPM_1_v2 = {1, 6.522099, 5.992113, 5.895799, 5.845471, 6.190067, 6.401603, 6.188478, 6.546069};

vector<float> calib_amp_8SiPM_2 = {1, 6.85418, 6.68469, 6.58283, 6.98329, 6.97906, 6.76493, 6.75924, 6.78279};
vector<float> calib_amp_8SiPM_2_v2 = {1, 6.813743, 6.468777, 6.242779, 6.607681, 6.813801, 6.603786, 6.476896, 6.634614};

vector<float> calib_amp_40SiPM_1_sw1 = {1, 5.051830, 4.945279, 4.925690, 4.749162, 4.949097, 4.723241, 4.707982, 4.868597};
vector<float> calib_amp_40SiPM_1_sw2 = {1, 4.819562, 4.740458, 4.655669, 4.563672, 4.697865, 4.539159, 4.613854, 4.748918};
vector<float> calib_amp_40SiPM_1_sw4 = {1, 4.682508, 4.566132, 4.469709, 4.511808, 4.556756, 4.569269, 4.513887, 4.730196};

vector<float> calib_amp_40SiPM_2_sw1 = {1, 4.504408, 4.321436, 4.349675, 4.474370, 4.773048, 4.893937, 4.735639, 4.620717};
vector<float> calib_amp_40SiPM_2_sw2 = {1, 4.427750, 4.289465, 4.360091, 4.501094, 4.655731, 4.804050, 4.808499, 4.657945};
vector<float> calib_amp_40SiPM_2_sw4 = {1, 4.396243, 4.217127, 4.344094, 4.416440, 4.678121, 4.678319, 4.633572, 4.705655};

// new procedure:
// 20 ns amplitude acceptance window
// constant baseline
// continuous gain fit

// __ AMP CALIB ___
vector<float> g_8SiPM_1_constBL = {1, 6.643101, 5.937994, 5.910624, 5.840343, 6.142987, 6.347942, 6.203394, 6.606222};

vector<float> g_8SiPM_2_constBL = {1, 6.838216, 6.477332, 6.248899, 6.621680, 6.813820, 6.551395, 6.484221, 6.685122};

vector<float> g_8SiPM_17_constBL = {1, 3.496053, 3.189863, 3.542560, 3.110981, 3.210311, 3.576601, 3.231899, 3.257059};

vector<float> g_8SiPM_17_constBL_HV58 = {1, 2.472151, 2.216387, 2.595906, 2.191874, 2.224226, 2.598138, 2.246901, 2.250155};
vector<float> g_8SiPM_17_constBL_HV59 = {1, 2.902316, 2.627994, 3.016038, 2.613367, 2.629431, 3.043153, 2.671950, 2.697508};
vector<float> g_8SiPM_17_constBL_HV59_alt = {1, 2.888785, 2.615228, 2.966872, 2.555194, 2.590216, 3.037215, 2.637975, 2.638252};
vector<float> g_8SiPM_17_dynBL_HV59 = {1, 2.937356, 2.624025, 3.017423, 2.613488, 2.639697, 3.069505, 2.650986, 2.712761};
vector<float> g_8SiPM_17_constBL_HV60 = {1, 3.477412, 3.172458, 3.543589, 3.105156, 3.209005, 3.558541, 3.210445, 3.255698};
vector<float> g_8SiPM_17_constBL_HV61 = {1, 3.828382, 3.523702, 3.865647, 3.438725, 3.538159, 3.862610, 3.525762, 3.600609};

vector<float> g_40SiPM_1_sw4_constBL = {1, 4.761221, 4.719708, 4.592255, 4.580244, 4.561942, 4.508876, 4.519476, 4.744483};
vector<float> g_40SiPM_1_sw2_Int590_constBL = {1, 4.800180, 4.712084, 4.578888, 4.521083, 4.783934, 4.477283, 4.640081, 4.785927};
vector<float> g_40SiPM_1_sw1_constBL = {1, 5.057210, 4.952975, 4.914601, 4.745170, 5.006344, 4.695190, 4.671079, 4.845346};

vector<float> g_40SiPM_2_sw4_constBL = {1, 4.367106, 4.278874, 4.324155, 4.502402, 4.770954, 4.676074, 4.684658, 4.699855};
vector<float> g_40SiPM_2_sw2_constBL = {1, 4.384197, 4.294646, 4.372352, 4.479390, 4.659250, 4.775097, 4.811665, 4.681038};
vector<float> g_40SiPM_2_sw1_constBL = {1, 4.484410, 4.304996, 4.311017, 4.451629, 4.779645, 4.873369, 4.759235, 4.637783};

// new, dedicated gain parameter fit
vector<float> g_8SiPM_0_constBL_amp_run65 = {1, 5.964026, 5.752018, 5.692648, 6.003136, 5.952826, 5.724358, 5.983295, 5.922015};
vector<float> g_8SiPM_1_constBL_amp_run68 = {1, 6.225833, 5.681876, 5.617020, 5.520674, 5.982826, 6.179563, 6.041097, 6.416068};
vector<float> g_8SiPM_2_constBL_amp_run69 = {1, 6.072533, 5.697019, 5.452204, 5.798762, 6.023438, 5.794798, 5.796922, 5.869892};

vector<float> g_8SiPM_17_constBL_amp_run82 = {1, 3.446521, 3.110892, 3.505680, 3.059208, 3.133277, 3.482466, 3.157271, 3.207298};
vector<float> g_8SiPM_17_constBL_amp_runah = {1, 3.256826, 3.336919, 3.907744, 3.564930, 3.305927, 3.617756, 3.259159, 3.260260};

vector<float> g_40SiPM_0_sw4_constBL_amp_run73 = {1, 4.083619, 4.010255, 3.864083, 4.043048, 4.044682, 4.110684, 4.044110, 4.130935};
vector<float> g_40SiPM_0_sw2_constBL_amp_run77 = {1, 4.478691, 4.436281, 4.234903, 4.442311, 4.609338, 4.369185, 4.386523, 4.407896};
vector<float> g_40SiPM_0_sw1_constBL_amp_run81 = {1, 4.429597, 4.342752, 4.202926, 4.356503, 4.468076, 4.301858, 4.375093, 4.394943};

vector<float> g_40SiPM_1_sw1_constBL_amp_run42 = {1, 4.983981, 4.816440, 4.796575, 4.629867, 4.833348, 4.602363, 4.607787, 4.662207};
vector<float> g_40SiPM_1_sw2_constBL_amp_run41 = {1, 4.723942, 4.647454, 4.632952, 4.420412, 4.712581, 4.396582, 4.549371, 4.741479};
vector<float> g_40SiPM_1_sw4_constBL_amp_run38 = {1, 4.697122, 4.561701, 4.429458, 4.466335, 4.486256, 4.545566, 4.497920, 4.660397};
vector<float> g_40SiPM_1_sw4_constBL_amp_run84 = {1, 4.661835, 4.543407, 4.393203, 4.356386, 4.440221, 4.389425, 4.484233, 4.662783};

vector<float> g_40SiPM_2_sw1_constBL_amp_run50 = {1, 4.454238, 4.233100, 4.244721, 4.388974, 4.683845, 4.767173, 4.663702, 4.586292};
vector<float> g_40SiPM_2_sw2_constBL_amp_run47 = {1, 4.407262, 4.213307, 4.280317, 4.388974, 4.584915, 4.629858, 4.728094, 4.627459};
vector<float> g_40SiPM_2_sw4_constBL_amp_run44 = {1, 4.452656, 4.115571, 4.178642, 4.432074, 4.627349, 4.619523, 4.515398, 4.479550};
vector<float> g_40SiPM_2_sw4_constBL_amp_run88 = {1, 3.939226, 3.737358, 3.876855, 3.971836, 4.207425, 4.142826, 4.179445, 4.226949};

// __ CHARGE CALIB ___
// integration window: 25 ns, Generalized Poisson fit
vector<float> g_8SiPM_0_constBL_charge_run65 = {1, 57.878024, 56.174099, 55.690767, 57.385766, 55.208148, 55.010960, 56.358564, 54.551059};
vector<float> g_8SiPM_1_constBL_charge_run68 = {1, 54.339372, 51.120311, 50.323608, 48.323768, 51.724367, 53.002368, 51.895161, 53.368556};
vector<float> g_8SiPM_2_constBL_charge_run69 = {1, 54.160940, 50.392792, 48.624219, 52.848405, 52.114772, 51.153844, 50.862783, 50.617176};

vector<float> g_8SiPM_17_constBL_charge_run82 = {1, 38.484672, 34.691483, 39.675637, 33.549637, 33.375852, 37.475023, 34.309823, 33.993859};
vector<float> g_8SiPM_17_constBL_charge_runah = {1, 32.295190, 33.063051, 39.319845, 35.306201, 33.005156, 36.094271, 32.217953, 32.376389};

vector<float> g_40SiPM_0_sw4_constBL_charge_run73 = {1, 43.573183, 42.598886, 42.130642, 42.583962, 41.787534, 41.779299, 42.792170, 42.464986};
vector<float> g_40SiPM_0_sw2_constBL_charge_run77 = {1, 44.797443, 45.136919, 43.947454, 45.311786, 44.719003, 43.283286, 43.864817, 44.257726};
vector<float> g_40SiPM_0_sw1_constBL_charge_run81 = {1, 42.185309, 43.821183, 42.567506, 43.421509, 42.907376, 41.566102, 43.004158, 43.267250};

vector<float> g_40SiPM_1_sw4_constBL_charge_run38 = {1, 45.144409, 44.690920, 43.004946, 42.404838, 42.360279, 41.824444, 42.600214, 42.955418};
vector<float> g_40SiPM_1_sw2_constBL_charge_run41 = {1, 45.511759, 45.039899, 43.961472, 42.708219, 43.968406, 41.951216, 43.145589, 43.788582};
vector<float> g_40SiPM_1_sw1_constBL_charge_run42 = {1, 47.756255, 46.178977, 46.196831, 43.313765, 45.733900, 43.988610, 43.394693, 44.405452};
vector<float> g_40SiPM_1_sw4_constBL_charge_run84 = {1, 44.267965, 43.887981, 42.866710, 42.386506, 42.066467, 41.266592, 42.462270, 42.703238};

vector<float> g_40SiPM_2_sw4_constBL_charge_run44 = {1, 41.686783, 39.788297, 41.507382, 41.046722, 43.547554, 42.832502, 43.277914, 42.493972};
vector<float> g_40SiPM_2_sw2_constBL_charge_run47 = {1, 41.426082, 39.606853, 40.685430, 40.492471, 43.440455, 43.619321, 44.750959, 42.655855};
vector<float> g_40SiPM_2_sw1_constBL_charge_run50 = {1, 41.749177, 39.963841, 40.124565, 40.850420, 44.759179, 44.674908, 44.686569, 43.093194};
vector<float> g_40SiPM_2_sw4_constBL_charge_run88 = {1, 37.792835, 36.457600, 37.816670, 37.611428, 39.822824, 39.078728, 39.895177, 39.592268};

//vector<float> charge_calib_vb56 = {};
vector<float> charge_calib_vb57 = {23.048371, 21.741578, 23.774282, 23.083660, 25.054184, 24.319185, 24.334173, 24.853267, 1};
vector<float> charge_calib_vb58 = {28.798605, 27.162817, 29.150066, 28.549446, 30.820398, 29.892134, 30.319136, 30.365344, 1};
vector<float> charge_calib_vb591 = {34.185421, 32.560238, 34.443731, 33.495544, 36.136513, 35.077469, 35.733271, 35.739918, 1};
vector<float> charge_calib_vb60 = {38.916981, 37.187023, 39.095598, 38.365510, 41.048950, 40.473847, 40.615454, 40.649090, 1};
vector<float> charge_calib_vb61 = {43.269773, 41.377581, 43.197679, 42.415529, 45.438782, 44.406172, 45.059923, 45.089259, 1};
vector<float> charge_calib_vb62 = {46.764727, 44.902501, 46.419689, 45.536503, 49.301141, 48.813421, 48.574390, 48.369936, 1};

vector<float> amp_calib_vb57 = {2.427135, 2.258541, 2.552066, 2.536811, 2.690861, 2.650630, 2.676413, 1.491813, 1};
vector<float> amp_calib_vb58 = {3.029364, 2.893179, 3.043664, 3.044143, 3.273233, 3.215380, 3.236108, 3.271589, 1};
vector<float> amp_calib_vb591 = {3.552258, 3.405994, 3.556098, 3.538558, 3.811056, 3.747463, 3.765869, 3.812031, 1};
vector<float> amp_calib_vb60 = {4.099530, 3.917634, 4.046946, 4.044494, 4.345492, 4.275376, 4.298305, 4.364655, 1};
vector<float> amp_calib_vb61 = {4.513963, 4.305106, 4.364560, 4.426066, 4.815006, 4.724756, 4.760836, 4.816021, 1};
vector<float> amp_calib_vb62 = {4.855707, 4.646365, 4.665608, 4.671475, 5.076868, 4.996825, 5.038053, 5.142009, 1};

vector<float> calib_dummy(16, 1);

vector<float> calib_amp;

vector<float> calib_charge;

/*
__ BASELINE SHIFT VALUES of indiv. arrays _______________________________
*/

//JAN, Channel 0-7 + Trigger

vector<float> BL_calib_vb56 = {0.131471, -0.062083, -0.125199, -0.053618, -0.062996, -0.107044, -0.092104, 0.024049, 0};
vector<float> BL_calib_vb57 = {0.133201, -0.092415, -0.163606, -0.078015, -0.104844, -0.151757, -0.126055, -0.007578, 0};
vector<float> BL_calib_vb58 = {0.068561, -0.137188, -0.204772, -0.118647, -0.177451, -0.226627, -0.179798, -0.077161, 0};
vector<float> BL_calib_vb591 = {-0.009987, -0.198874, -0.276064, -0.198530, -0.408925, -0.409971, -0.355443, -0.265315, 0};
vector<float> BL_calib_vb60 = {-0.248426, -0.381159, -0.511187, -0.409334, -0.617221, -0.603323, -0.513275, -0.419863, 0};
vector<float> BL_calib_vb61 = {-0.417435, -0.516524, -0.720246, -0.567681, -0.833585, -0.814164, -0.728049, -0.627355, 0};
vector<float> BL_calib_vb62 = {-0.615492, -0.677449, -0.965486, -0.784299, -1.169815, -1.071159, -0.998565, -0.912138, 0};

vector<float> BL_calib_vb60_tune90 = {-0.381951, -0.301294, -0.999940, -1.005571, -1.112136, -1.225693, -1.708072, -1.112632, 0};
vector<float> BL_calib_vb60_tune90_schablone = {-0.436188, -1.267119, -1.494210, -1.493263, -1.475991, -1.889365, -2.526676, -2.012294, 0};

vector<float> BL_dummy(9, 0);

vector<float> BL_const;

// SWITCH dynamic <-> constant baseline
bool switch_BL = false; // true = dyn, false = const
//Integration Window
bool isDC = true;
//IF the calibration values are correct, otherwise use dummies
bool isCalib = true;

int amp_array_printRate = 1000;
int wavesPrintRate = 1000;
// clibrated sum

int ch0PrintRate = 1000;
int trigPrintRate = 1000;   //100
int signalPrintRate = 1000; //100
double coef = 2.5 / (4096 * 10);
string WCHU("AB"), WCAlexander("CD");

//External Variables - mostly definded in main.C
extern string WCVersion;
extern int runNr;
extern float horizontal;
extern float vertical;
extern float angle;
extern int pdgID;
extern float energy;
extern int isSP;
extern int mp;
extern int safPMT2;
extern int safPMT1;
extern int safSiPM;
extern int trackL;

void readFull(TString _inFileList, TString _inDataFolder, TString _outFile, string _runName, string dynamicBL_, string isDC_, string useCalibValues_)
{

  if (dynamicBL_ == "0")
  {
    switch_BL = false;
  }
  else
  {
    switch_BL = true;
  }
  if (isDC_ == "0")
  {
    isDC = false;
  }
  else
  {
    isDC = true;
  }
  if (useCalibValues_ == "0")
  {
    isCalib = false;
  }
  else
  {
    isCalib = true;
  }
  cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
  cout << "DYNAMIC BASELINE: " << switch_BL << endl;
  cout << "isDC: " << isDC << endl;
  cout << "use Calib Values: " << isCalib << endl;
  cout << "runName: " << _runName << endl;

  cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;

  /*
  SWITCH ANPASSEN
  Number of Parameters has to be defined in read.h, main.C and runList.sh AND in the number of parameters in main.c (if case)


  /*
  Function that is used by main(). read() does most of the analysis of the data. The raw .bin files 
  are read, converted into root-histograms and then the analysis calls are done. Most of the longer
  functions are defined in alanysis.C.
  The read() function the saves all the events (event by event) of that particular run to a root tree
  which is the saved in /runs/runName/out.root.
  */
  TF1 *f_const;
  vector<float> calib_amp = calib_dummy;

  vector<float> calib_charge = calib_dummy;

  if ((strcmp(_runName.c_str(), "calib_vb60_tune90_schablone") == 0))
  {
    BL_const = BL_calib_vb60_tune90_schablone;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb60_tune90") == 0) )
  {
    BL_const = BL_calib_vb60_tune90;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb56") == 0) || (strcmp(_runName.c_str(), "dc_vb56") == 0))
  {
    BL_const = BL_calib_vb56;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb57") == 0) || (strcmp(_runName.c_str(), "dc_vb57") == 0))
  {
    BL_const = BL_calib_vb57;
    calib_amp = amp_calib_vb57;
    calib_charge = charge_calib_vb57;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb58") == 0) || (strcmp(_runName.c_str(), "dc_vb58") == 0))
  {
    BL_const = BL_calib_vb58;
    calib_amp = amp_calib_vb58;
    calib_charge = charge_calib_vb58;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb591") == 0) || (strcmp(_runName.c_str(), "dc_vb591") == 0))
  {
    BL_const = BL_calib_vb591;
    calib_amp = amp_calib_vb591;
    calib_charge = charge_calib_vb591;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb60") == 0) || (strcmp(_runName.c_str(), "dc_vb60") == 0))
  {
    BL_const = BL_calib_vb60;
    calib_amp = amp_calib_vb60;
    calib_charge = charge_calib_vb60;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb61") == 0) || (strcmp(_runName.c_str(), "dc_vb61") == 0))
  {
    BL_const = BL_calib_vb61;
    calib_amp = amp_calib_vb61;
    calib_charge = charge_calib_vb61;
  }
  else if ((strcmp(_runName.c_str(), "calib_vb62") == 0) || (strcmp(_runName.c_str(), "dc_vb62") == 0))
  {
    BL_const = BL_calib_vb62;
    calib_amp = amp_calib_vb62;
    calib_charge = charge_calib_vb62;
  }
  else
  {
    //Falscher Name
    BL_const = BL_dummy;
    calib_amp = calib_dummy;
    calib_charge = calib_dummy;
  }
  if (switch_BL)
  {
    BL_const = BL_dummy;
  }
  if (!isCalib)
  {
    calib_amp = calib_dummy;
    calib_charge = calib_dummy;
  }
  //Aktivieren, falls keine Baselines gemessen, z.B. für das erstmalige Auslesen der ROOT Daten
  // BL_const=BL_dummy;

  if (isDC)
  {
    amp_array_printRate = 1000;
    wavesPrintRate = 1000;
  }
  else
  {
    amp_array_printRate = 2000;
    wavesPrintRate = 2000;
  }

  printf(("USED BL ARRAY: %f,%f,%f,%f,%f,%f,%f,%f FOR: " + _runName + "\n").c_str(), BL_const[0], BL_const[1], BL_const[2], BL_const[3], BL_const[4], BL_const[5], BL_const[6], BL_const[7], BL_const[8]);
  printf(("USED AMP CALIB ARRAY: %f,%f,%f,%f,%f,%f,%f,%f FOR: " + _runName + "\n").c_str(), calib_amp[0], calib_amp[1], calib_amp[2], calib_amp[3], calib_amp[4], calib_amp[5], calib_amp[6], calib_amp[7], calib_amp[8]);

  printf("IS DARKCOUNT: %s | Dynamic Baseline: %s | Is Calibrated: %s ", btoa(isDC), btoa(switch_BL), btoa(isCalib));

  /*Create root-file and root-tree for data*/
  TFile *rootFile = new TFile(_outFile, "RECREATE");
  if (rootFile->IsZombie())
  {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
         << _outFile << ": check that the path is correct!!!"
         << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");

  tree->SetAutoSave(0LL);
  tree->SetAutoFlush(0LL);

  TTree::SetBranchStyle(0);

  gStyle->SetLineScalePS(1); // high resolution plots

  /*Declare & define the variables that are to be saved in the root-tree or that are used during the analysis.*/
  Int_t EventNumber = -999;
  Int_t LastEventNumber = -999;
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
  Float_t tPMT1 = -999;
  Float_t tPMT2 = -999;
  Float_t tPMT2i = -999;
  Float_t tSUMp = -999;
  Float_t tSUMm = -999;
  Float_t trigTp = -999; //t_trig' = [(t0+t1)-(t2+t3)]/4
  Float_t t0t1 = -999;   //t0t1 = [(t0-t1)]
  Float_t t2t3 = -999;   //t2t3 = [(t2-t3)]
  Int_t isVeto = -999;   //variable to define veto, 1 if veto, 0 if not, -999 if undefined
  Int_t isTrig = -999;

  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  Int_t WOMID[16]; //1=A, 2=B, 3=C, 4=D

  float amp_array;       // calibrated amplitude, baseline-shifted sum signal
  float charge_array;    // calibrated charge, baseline-shifted sum signal
  float t_amp_array;     // point in time of sum signal
  float chPE_amp[16];    // single channel amplitude at sum signal
  float chPE_charge[16]; // single channel charge at sum signal

  Float_t t_amp_array_invCFD = -999;
  Float_t t_amp_array_invCFD_wrtTrig = -999;

  std::vector<float> amp(16, -999);
  std::vector<float> amp_inRange(16, -999);
  std::vector<float> max(16, -999);
  std::vector<float> min(16, -999);
  Float_t t[16];
  Float_t tSiPM[16];

  float Integral_0_300[16]; //array used to store Integral of signal from 0 to 300ns
  float Integral_inRange[16];
  float Integral[16];
  float Integral_mVns[16];

  float BL_output[4];        //array used for output getBL-function
  Float_t BL_lower[16];      //store baseline for 16 channels for 0-75ns range
  Float_t BL_RMS_lower[16];  //store rms of baseline for 16 channels for 0-75ns range
  Float_t BL_Chi2_lower[16]; //store chi2/dof of baseline-fit for 16 channels for 0-75ns range
  Float_t BL_pValue_lower[16];
  Float_t BL_upper[16];      //store baseline for 16 channels for 220-320ns range
  Float_t BL_RMS_upper[16];  //store rms of baseline for 16 channels for 220-320ns range
  Float_t BL_Chi2_upper[16]; //store chi2/dof of baseline-fit for 16 channels for 220-320ns range
  Float_t BL_pValue_upper[16];

  Float_t BL_used[16];
  Float_t BL_Chi2_used[16];
  Float_t BL_pValue_used[16];

  int nPeaks = 4; // maximum number of peaks to be stored by peakfinder; has to be set also when creating branch
  Double_t peakX[16][nPeaks];
  Double_t peakY[16][nPeaks];

  int NumberOfBins;
  Int_t EventIDsamIndex[16];
  Int_t FirstCellToPlotsamIndex[16];
  std::vector<TH1F *> hChSum;
  for (int i = 0; i < 16; i++)
  {
    TString name("");
    name.Form("hChSum_%d", i);
    TH1F *h = new TH1F("h", ";ns;Amplitude, mV", 1024, -0.5 * SP, 1023.5 * SP);
    h->SetName(name);
    hChSum.push_back(h);
  }
  std::vector<TH1F *> hChShift;
  for (int i = 0; i < 16; i++)
  {
    TString name("");
    name.Form("hChShift_%d", i);
    TH1F *h = new TH1F("h", ";ns;Amplitude, mV", 1024, -0.5 * SP, 1023.5 * SP);
    h->SetName(name);
    hChShift.push_back(h);
  }
  std::vector<TH1F> hChtemp;
  for (int i = 0; i < 16; i++)
  {
    TString name("");
    name.Form("hChtemp_%d", i);
    TH1F h("h", ";ns;Amplitude, mV", 1024, -0.5 * SP, 1023.5 * SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }
  std::vector<TH1F> hChShift_temp;
  for (int i = 0; i < 16; i++)
  {
    TString name("");
    name.Form("hChShift_temp_%d", i);
    TH1F h("h", ";ns;Amplitude, mV", 1024, -0.5 * SP, 1023.5 * SP);
    h.SetName(name);
    hChShift_temp.push_back(h);
  }
  Short_t amplValues[16][1024];
  TH1F hCh("hCh", "dummy;ns;Amplitude, mV", 1024, -0.5 * SP, 1023.5 * SP);

  TString plotSaveFolder = _outFile;
  plotSaveFolder.ReplaceAll((_runName + ".root").data(), "");

  TCanvas cWaves("cWaves", "cWaves", 1000, 700);
  cWaves.Divide(4, 4);
  TCanvas csumWOMA("csumWOMA", "csumWOMA", 1000, 700);
  csumWOMA.Divide(4, 2);
  TCanvas csumWOMB("csumWOMB", "csumWOMB", 1000, 700);
  csumWOMB.Divide(3, 3);
  TCanvas cCh0("cCh0", "cCh0", 1500, 900);
  cCh0.Divide(2, 2);
  TCanvas cTrig("cTrig", "cTrig", 1500, 900);
  cTrig.Divide(2, 2);
  TCanvas cSignal("cSignal", "cSignal", 1500, 900);
  cSignal.Divide(2, 2);

  // clibrated sum
  TCanvas C_amp_array("C_amp_array", "C_amp_array", 1000, 700);
  C_amp_array.Divide(3, 3);

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
  tree->Branch("tPMT1", &tPMT1, "tPMT1/F");
  tree->Branch("tPMT2", &tPMT2, "tPMT2/F");
  tree->Branch("tPMT2i", &tPMT2i, "tPMT2i/F");
  tree->Branch("tSUMp", &tSUMp, "tSUMp/F");
  tree->Branch("tSUMm", &tSUMm, "tSUMm/F");
  tree->Branch("runNr", &runNr, "runNr/I");      //run number in google table
  tree->Branch("horiz", &horizontal, "horiz/F"); // horizontal position of the box units: [cm]
  tree->Branch("vert", &vertical, "vert/F");     //vertical position of the box, units: [cm]
  tree->Branch("angle", &angle, "angle/F");
  tree->Branch("pdgID", &pdgID, "pdgID/I");
  tree->Branch("energy", &energy, "energy/F");
  tree->Branch("isSP", &isSP, "isSP/I");
  tree->Branch("mp", &mp, "mp/I");
  tree->Branch("safPMT2", &safPMT2, "safPMT2/I"); //solid angle factor
  tree->Branch("safPMT1", &safPMT1, "safPMT1/I"); //solid angle factor
  tree->Branch("safSiPM", &safSiPM, "safSiPM/I"); //solid angle factor
  tree->Branch("trackL", &trackL, "trackL/I");    //track length
  tree->Branch("isLastEvt", &isLastEvt, "isLastEvt/I");
  tree->Branch("trigGate", &trigGate, "trigGate/F");
  tree->Branch("trigTp", &trigTp, "trigTp/F");
  tree->Branch("t0t1", &t0t1, "t0t1/F"); //t0t1 = [(t0-t1)]
  tree->Branch("t2t3", &t2t3, "t2t3/F");
  tree->Branch("isVeto", &isVeto, "isVeto/I");
  tree->Branch("isTrig", &isTrig, "isTrig/I");
  tree->Branch("isGoodSignal_5", &isGoodSignal_5, "isGoodSignal_5/I");

  // CHANNEL INFO (but everything that is nCH-dependend below)
  tree->Branch("nCh", &nCh, "nCh/I");
  tree->Branch("WOMID", WOMID, "WOMID[nCh]/I");
  tree->Branch("ch", ChannelNr, "ch[nCh]/I");
  // AMPLITUDE
  tree->Branch("amp", amp.data(), "amp[nCh]/F");
  tree->Branch("amp_inRange", amp_inRange.data(), "amp_inRange[nCh]/F");
  tree->Branch("max", max.data(), "max[nCh]/F");
  tree->Branch("min", min.data(), "min[nCh]/F");
  // INTEGRAL
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral_inRange", Integral_inRange, "Integral_inRange[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F");
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
  // PEAKFINDER
  tree->Branch("peakX", peakX, "peakX[nCh][4]/D");
  tree->Branch("peakY", peakY, "peakY[nCh][4]/D");
  // CALIBRATED SUM
  tree->Branch("t_amp_array_invCFD", &t_amp_array_invCFD, "t_amp_array_invCFD/F");
  tree->Branch("t_amp_array_invCFD_wrtTrig", &t_amp_array_invCFD_wrtTrig, "t_amp_array_invCFD_wrtTrig/F");
  tree->Branch("amp_array", &amp_array, "amp_array/F");
  tree->Branch("t_amp_array", &t_amp_array, "t_amp_array/F");
  tree->Branch("charge_array", &charge_array, "charge_array/F");
  tree->Branch("chPE_amp", chPE_amp, "chPE_amp[nCh]/F");
  tree->Branch("chPE_charge", chPE_charge, "chPE_charge[nCh]/F");

  tree->Branch("EventIDsamIndex", EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex", FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  /*Start reading the raw data from .bin files.*/
  int nitem = 1;
  ifstream inList;
  TString fileName;
  inList.open(_inFileList);
  assert(inList.is_open());

  // clibrated sum
  int amp_array_PrintStatus = -1;

  int wavePrintStatus = -1;
  int sumWOMAPrintStatus = -1;
  int sumWOMBPrintStatus = -1;
  int ch0PrintStatus = -1;
  int trigPrintStatus = -1;
  int signalPrintStatus = -1;
  while (inList >> fileName)
  {
    fileName = _inDataFolder + fileName;
    cout << endl;
    cout << fileName << endl;
    FILE *pFILE = fopen(fileName.Data(), "rb");
    if (pFILE == NULL)
    {
      fputs("File error", stderr);
      assert(0);
    }
    fseek(pFILE, 0, SEEK_END);
    int totFileSizeByte = ftell(pFILE);
    rewind(pFILE);
    cout << "totFileSizeByte = " << totFileSizeByte << endl;
    int size_of_header;
    /*During 2018 testbeam measurements two WaveCatchers were used. One from the Berlin group
    and one from Alexander from Geneva. As these two Wavecatchers had two different versions
    there are two types of raw data files that have different header lengths.*/
    if (WCVersion == WCHU)
    {
      size_of_header = 328;
    }
    else if (WCVersion == WCAlexander)
    {
      size_of_header = 327;
    }
    char header[size_of_header];
    nitem = fread(header, 1, size_of_header, pFILE);

    cout << "Header:\n"
         << header << endl;

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

    if (nActiveCh > 9)
    {
      cout << endl;
      char dummy;
      nitem = fread(&dummy, 1, 1, pFILE);
    }

    int whileCounter = 0;
    /*Loop over events. Events are processed and analysed one by one in order.*/
    while (nitem > 0)
    { //event loop
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
      if (WCVersion == WCHU)
      {
        nitem = fread(&nCh, sizeof(unsigned int), 1, pFILE); // since V2.8.14 the number of stored channels is written for each event
      }
      else if (WCVersion == WCAlexander)
      {
        nCh = 16;
      }

      if (EventNumber % 100 == 0)
      {
        printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber, Year, Month, Day, Hour, Minute, Second, Millisecond, nActiveCh, nCh);
      }

      float MeasuredBaseline[16];
      float AmplitudeValue[16];
      float ComputedCharge[16];
      float RiseTimeInstant[16];
      float FallTimeInstant[16];
      float RawTriggerRate[16];
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

        /*
        __ Set WOMID _________________________________________________________
        The labeling of the WOMs in the box was done using the letters A,B,C,D. For convinience these letters are here replaced by the numbers 1-4 which is stored in the root-tree for every channel and every event.
        */
        if (WCVersion == WCAlexander)
        {
          if (i <= 6)
          {
            WOMID[i] = 3;
          }
          else if (i >= 7 && i <= 14)
          {
            WOMID[i] = 4;
          }
        }
        else
        {
          if (i <= 6)
          {
            WOMID[i] = 1;
          }
          else if (i >= 7 && i <= 14)
          {
            WOMID[i] = 2;
          }
        }

        TString title("");
        title.Form("ch %d, ev %d", i, EventNumber);
        hCh.Reset();
        hCh.SetTitle(title);

        /*
        __ Waveform Histogram _______________________________________________
        Writing the signal amplitude values into the root-histogram hCh.
        */
        if (i == 15)
        {
          for (int j = 0; j < 1024; j++)
          {
            nitem = fread(&amplValues[i][j], sizeof(short), 1, pFILE);
            hCh.SetBinContent(j + 1, -(amplValues[i][j] * coef * 1000));
          }
        }
        else
        {
          for (int j = 0; j < 1024; j++)
          {
            nitem = fread(&amplValues[i][j], sizeof(short), 1, pFILE);
            hCh.SetBinContent(j + 1, (amplValues[i][j] * coef * 1000));
          }
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

        /*
        __ Baseline Fit _______________________________________________________
        Calculate baseline values infront and after the triggered signal
        Triggered signal is expected in the range fromm 100 to 150 ns
        */
        BL_fit(&hChtemp.at(i), BL_output, 0.0, 75.0);
        BL_lower[i] = BL_output[0];
        BL_RMS_lower[i] = BL_output[1];
        BL_Chi2_lower[i] = BL_output[2];
        BL_pValue_lower[i] = BL_output[3];
        BL_fit(&hChtemp.at(i), BL_output, 220.0, 320.0);
        BL_upper[i] = BL_output[0];
        BL_RMS_upper[i] = BL_output[1];
        BL_Chi2_upper[i] = BL_output[2];
        BL_pValue_upper[i] = BL_output[3];

        // determine "best" baseline
        if (BL_Chi2_upper[i] <= BL_Chi2_lower[i])
        {
          BL_used[i] = BL_upper[i];
          BL_Chi2_used[i] = BL_Chi2_upper[i];
          BL_pValue_used[i] = BL_pValue_upper[i];
        }
        else
        {
          BL_used[i] = BL_lower[i];
          BL_pValue_used[i] = BL_pValue_lower[i];
          BL_Chi2_used[i] = BL_Chi2_lower[i];
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

        /*
        __ Peakfinder _________________________________________________________
        Implemented to search double-muon-event candiates
        Set maximum number of peaks stored in beginning of script -> nPeaks
        peakX/Yarray[nCh][nPeaks] stores peak coordinates as branches in tree
        Switch on/off with pfON
        -> when off:  set peakX/Yarray[nCh][nPeaks] to zero
        */
        gErrorIgnoreLevel = kError; // suppress root terminal output

        bool pfON = false;
        if (i < 15)
        {
          pfON = false;
        }                     // switch on/off peakfinder
        int sigma = 10;       // sigma of searched peaks
        Double_t thrPF = 0.1; // peakfinder threshold
        TPolyMarker pm;       // store polymarker showing peak position, print later
        peakfinder(&hCh, 0, 130, nPeaks, sigma, thrPF, peakX[i], peakY[i], &pm, pfON);

        gErrorIgnoreLevel = kUnset; // return to normal terminal output

        // baseline-correct Y-values and convert to units of p.e.
        if (pfON)
        {
          for (int j = 0; j < nPeaks; ++j)
          {
            peakY[i][j] = amp2pe(peakY[i][j], calib_amp[i], BL_shift);
          }
        }

        // printf("X: %d %f %f %f %f \n",i,peakX[i][0],peakX[i][1],peakX[i][2],peakX[i][3]);
        // printf("Y: %d %f %f %f %f \n",i,peakY[i][0],peakY[i][1],peakY[i][2],peakY[i][3]);

        /*
        __ CFD _____________________________________________________________
        Setting the signal time by using a constant fraction disriminator method.
        The SiPM and the trigger sinals are handled differently using different thresholds.
        */
        if (i == 9)
        { //trigger
          t[i] = CDF(&hCh, 0.5);
        }
        else
        { //SiPMs
          t[i] = CFD2(&hCh, 0.35);
          if (t[i] < 95)
          {
            t[i] = CFDinvert2(&hCh, 0.35);
          }
        }

        /*
        __Print Raw Data to .txt ______________________________________________
        Select channel. Prints histogram raw data in a two column text file
        */

        // if (i==4 && BL_chi2[4]<1.7 && BL_chi2[4]>0.7)
        // if (i==4)
        // {
        //   TString histDataName;
        //   histDataName.Form("Ch%d_hist_data.txt",i);
        //   TString path2hist_data;
        //   path2hist_data.Form("%s/%s",(const char*)plotSaveFolder,(const char*)histDataName);
        //   FILE * histOut;
        //   histOut = fopen(path2hist_data,"a"); // produces overhead, maybe put this infront of loop

        //   Int_t nbins_x = hCh.GetNbinsX(); // bins at k==0 and k==nbins_x seem to have BinContent==0
        //   for (Int_t k=1; k<=nbins_x; k++)
        //   {
        //     fprintf(histOut,"%.4f %.8f\n",
        //     hCh.GetBinLowEdge(k)+hCh.GetBinWidth(k)/2,
        //     hCh.GetBinContent(k));
        //   }
        //   fclose(histOut);
        // }

        // clibrated sum
        // if(EventNumber%amp_array_printRate==0 && (i!=9)){
        if (EventNumber % amp_array_printRate == 0 && (i < 8))
        {
          C_amp_array.cd(i + 1);
          hCh.DrawCopy();
        }

        /*
        __ Integral & Amplitude ________________________________________
        There are several definitions of the integral of a signal used here. Those are:
        - Integral_0_300: Integration over the entire time window (~320ns)
        - Integral: Integration over a smaller time window (~50ns) relative to the trigger
        Additionally the number of p.e. is now calculated using the amplitude
        and the calibration factors in the calib_amp-vactor. The function 'PE' calculates the amplitude of the signal, subtracts the better BL value and divides by the calibration factor.
        */
        // Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-0.0*1024*SP);

       // float integralStart = 145;
        //float integralEnd = 185;

        Integral[i] = Integrate_50ns(&hCh, BL_shift) / calib_charge.at(i); // difined 50 ns window

        // calibrated, BL-shifted charge
        if (isDC)
          Integral_inRange[i] = integral(&hCh, 50, 75, BL_shift) / calib_charge.at(i); // for DC runs
        else
          //Integral_inRange[i] = integral(&hCh, integralStart, integralEnd, BL_shift) / calib_charge.at(i);
         Integral_inRange[i] = integral(&hCh, 145, 170, BL_shift) / calib_charge.at(i);

        // calibrated, BL-shifted amplitude at maximum in window
        if (isDC)
          amp[i] = PE(&hCh, calib_amp.at(i), BL_shift, 50.0, 100.0); // for DC runs
        else
        //  amp[i] = PE(&hCh, calib_amp.at(i), BL_shift, integralStart, integralEnd);
         amp[i] = PE(&hCh, calib_amp.at(i), BL_shift, 135.0, 185.0);
        // reduced window
        if (isDC)
          amp_inRange[i] = PE(&hCh, calib_amp.at(i), BL_shift, 50.0, 75.0); // for DC runs
        else
         // amp_inRange[i] = PE(&hCh, calib_amp.at(i), BL_shift, integralStart, integralEnd);
        amp_inRange[i] = PE(&hCh, calib_amp.at(i), BL_shift, 145.0, 165.0);

        /*
        __ Printing Wafevorms ____________________________________________
        The signals for events can be printed to a .pdf file called waves.pdf. The rate at which the events are drawn to waves.pdf is set via the variable wavesPrintRate. Additional requirements can be set in the if-statement to look at specific events only.
        The entire if-statement so far also plots lines at the found signal maximum, the corresponding integration limit, as well as the BL values to each of the histograms.
        */
        if (EventNumber % wavesPrintRate == 0)
        {
          cWaves.cd(1 + 4 * (i % 4) + (i) / 4);
          hCh.DrawCopy();
          hCh.GetXaxis()->SetRange((t[i] - 20) / SP, (t[i] + 30) / SP);
          int max_bin = hCh.GetMaximumBin();
          int lower_bin = max_bin - 20.0 / SP;
          int upper_bin = max_bin + 30.0 / SP;
          // double x = h->GetXaxis()->GetBinCenter(binmax);
          float max_time = hCh.GetXaxis()->GetBinCenter(max_bin);
          float lower_time = hCh.GetXaxis()->GetBinCenter(lower_bin);
          float upper_time = hCh.GetXaxis()->GetBinCenter(upper_bin);
          hCh.GetXaxis()->SetRange(0, 1024);
          TLine *ln4 = new TLine(0, BL_lower[i], 75, BL_lower[i]);
          TLine *ln5 = new TLine(220, BL_upper[i], 320, BL_upper[i]);
          TText *text = new TText(.5, .5, Form("%f %f", BL_lower[i], BL_upper[i]));
          ln4->SetLineColor(2);
          ln5->SetLineColor(2);
          ln4->Draw("same");
          ln5->Draw("same");
          text->Draw("same");
          if (pfON)
          {
            pm.Draw();
          } // print peakfinders polymarker
        }
        // End of loop over inividual channels
      }

      // continue;

      /*
      __ TIMING _____
      */
      trigT = t[9];
      for (int i = 0; i <= 15; i++)
      {
        tSiPM[i] = t[i] - trigT;
        /*
        if (tSiPM[i+7] < -66){
          t[i+7] = CDFinvert(&hChtemp.at(i+7),0.33);
          tSiPM[i+7] = t[i+7] - trigT;
        }
        */
      }

      /*
      __ FILLING SUM HISTOGRAMS ________________________
      __ Calibrated SUM SIGNAL of all channels ______________
      */
      TH1F hSum("hSum", "Calibrated Sum ;ns;Amplitude in npe", 1024, -0.5 * SP, 1023.5 * SP);
      for (int hSumIndex = 0; hSumIndex < 8; hSumIndex++)
      {
        // SWITCH dynamic <-> constant baseline
        float BL_shift;
        if (switch_BL)
        {
          BL_shift = BL_used[hSumIndex];
        }
        else
        {
          BL_shift = BL_const[hSumIndex];
        }
        free(f_const);
        f_const = new TF1("f_const", "pol0", 0, 320);
        f_const->SetParameter(0, BL_shift);

        hChtemp.at(hSumIndex).Add(f_const, -1);
        hChtemp.at(hSumIndex).Scale(1.0 / calib_amp.at(hSumIndex));

        hSum.Add(&hChtemp.at(hSumIndex), 1);
      }

      // get point of amplitude maximum in 50 ns window
      amp_array = max_inRange(&hSum, 50., 100.0);
      t_amp_array = t_max_inRange(&hSum, 50.0, 100.0);
      // amp_array = 8.0/7.0 * max_inRange(&hSum, 100.0, 150.0);

      C_amp_array.cd(9);
      hSum.DrawCopy();

      // get single channel amplitude/charge at time of sum maximum
      for (int i = 0; i < 8; i++)
      {

        chPE_amp[i] = amp_atTime(&hChtemp.at(i), t_amp_array);
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
        // reverse amplitude calibration before integration
        hChtemp.at(i).Scale(calib_amp.at(i));
        chPE_charge[i] = integral(&hChtemp.at(i), t_amp_array - 10, t_amp_array + 15, 0) / calib_charge.at(i);
        // chPE_charge[i] = integral(&hChtemp.at(i), t_amp_array-10, t_amp_array+15, BL_shift) / calib_charge.at(i);
      }
      // summed charge
      charge_array = chPE_charge[0] + chPE_charge[1] + chPE_charge[2] + chPE_charge[3] + chPE_charge[4] + chPE_charge[5] + chPE_charge[6] + chPE_charge[7];
      // timing
      t_amp_array_invCFD = CFDinvert2(&hSum, 0.4);
      t_amp_array_invCFD_wrtTrig = trigT - t_amp_array_invCFD;

      /* end */

      /*Saving the plotted signals/events to a new page in the .pdf file.*/
      if (EventNumber % wavesPrintRate == 0)
      {
        if (wavePrintStatus < 0)
        {
          cWaves.Print((TString)(plotSaveFolder + "/waves.pdf("), "pdf");
          wavePrintStatus = 0;
        }
        else
          cWaves.Print((TString)(plotSaveFolder + "/waves.pdf"), "pdf");
      }
      if (EventNumber % trigPrintRate == 0)
      {
        if (trigPrintStatus < 0)
        {
          cTrig.Print((TString)(plotSaveFolder + "/trig.pdf("), "pdf");
          trigPrintStatus = 0;
        }
        else
          cTrig.Print((TString)(plotSaveFolder + "/trig.pdf"), "pdf");
      }
      if (EventNumber % signalPrintRate == 0)
      {
        if (signalPrintStatus < 0)
        {
          cSignal.Print((TString)(plotSaveFolder + "/signal.pdf("), "pdf");
          signalPrintStatus = 0;
        }
        else
          cSignal.Print((TString)(plotSaveFolder + "/signal.pdf"), "pdf");
      }

      // clibrated sum
      if (EventNumber % amp_array_printRate == 0)
      {
        if (amp_array_PrintStatus < 0)
        {
          C_amp_array.Print((TString)(plotSaveFolder + "/amp_array.pdf("), "pdf");
          amp_array_PrintStatus = 0;
        }
        else
          C_amp_array.Print((TString)(plotSaveFolder + "/amp_array.pdf"), "pdf");
      }

      /*Writing the data for that event to the tree.*/
      tree->Fill();

      // cout<<"SIZE OF: "<<sizeof(amp_array)<< endl;
      //tree->Print();
    }
    auto nevent = tree->GetEntries();

    cout << "EVENTS:  " << nevent << endl;

    fclose(pFILE);
  }
//tree->Print();

  /*Clearing objects and saving files.*/
  inList.close();
  cWaves.Clear();
  cWaves.Print((TString)(plotSaveFolder + "/waves.pdf)"), "pdf");
  cCh0.Print((TString)(plotSaveFolder + "/ch0.pdf)"), "pdf");
  cTrig.Print((TString)(plotSaveFolder + "/trig.pdf)"), "pdf");
  cSignal.Print((TString)(plotSaveFolder + "/signal.pdf)"), "pdf");

  // clibrated sum
  C_amp_array.Clear();
  C_amp_array.Print((TString)(plotSaveFolder + "/amp_array.pdf)"), "pdf");

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}