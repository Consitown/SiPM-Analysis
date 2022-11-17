#ifndef ANALYSIS
#define ANALYSIS

#include <TList.h>
#include <TString.h>
#include <TPolyMarker.h>

using namespace std;

string checkFilename(TString filename);
// TIMING
float CFD(TH1F* hWave,float thr);
float CFDNegative(TH1F* hWave,float thr);
float CFDNegativeCustom(TH1F* hWave,float thr, int negativeOffset);
float CFDInRange(TH1F* hWave,float thr,float start,float end);
float CFDinvert(TH1F* hWave,float thr);
float CFDinvertInRange(TH1F* hWave,float thr,float start,float end);
float CFDNegativeInvert(TH1F *hWave, float thr);
// INTEGRAL
float IntegralHist(TH1F* hWave,float t1,float t2,float BL);
float IntegralHistCFD(TH1F* hWave,float trheshold,float window,float BL);
float IntegralDifference(TH1F* hWave,float leftStart ,float rightEnd,float rightEndAll,float Amplitude,float BL);
float IntegralTimeOverThreshold(TH1F* hWave,float startSearch,float endSearch,float threshold,float BL);

// BASELINE
float* getBL(TH1F* hWave, float* BL, float t1, float t2);
float* BL_fit(TH1F* hWave, float* BL_chi2, float t1, float t2);
// AMPLITUDE
float AmplitudeHist(TH1F* hWave,  float t1, float t2,float BL);
float AmplitudeHistAlternative(TH1F* hWave, float t1, float t2,float BL);
float max_inRange(TH1F* hWave,float t1, float t2);
float t_max_inRange(TH1F* hWave,float t1, float t2);
float amp_atTime(TH1F* hWave,float t_max);
double amp2pe_u_l(double y, float calib_factor, float BL_upper, float BL_lower, float BL_Chi2_upper, float BL_Chi2_lower);
double amp2pe(double y, float calib_factor, float BL_used);
// OTHER
double correction_function(double x);
void peakfinder(TH1F *hWave, float t1, float t2, int nPeaks, int sigma, double thr, double *Xarray, double *Yarray, TPolyMarker *pfMarker, bool pfON);
float estimateNL(TH1F *hWave, float t);




#endif