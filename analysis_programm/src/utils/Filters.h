/// Helper functions
#ifndef _Filters
#define _Filters

#include <TVirtualFFT.h>
#include <TComplex.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH1.h>

#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "Helpers.h"

using namespace std;

class Filters {
public:
	// convolution for filtering waveforms
	static void Convolute(double*&, double*, double*, int);
	static void Convolute(double*&, double*, int);
	static void Deconvolute(double*&, double*, double*, int, double, double = 0., double = 0.);

	// smoothing
	static void SmoothArray(double*&, int, double = .625, string = "Gaus", double = .3125, double = 1.5);
	static void SmoothArray(double*&, int, double, int, double = .3125, double = 1.5);

	static void BoxFilter(double*&, int, int);
	static void GausFilter(double*&, int, double, double);
	static void GausFFTFilter(double*&, int, double, double);

	static void BilateralFilter(double*&, int, double, double, double);
	static void Bilateral2Filter(double*&, int, int, double, double, double);
	static void MedianFilter(double*& ar, int, int);

	// custom filter emulating primitive response function
	static void ResponseFilter(double*&, int, double = .4, double= 1.2, double= .25, double = .3125);
	static void SecondOrderUnderdampedFilter(double*&, int, double, double, double, double = .3125, bool = false);
	static void SecondOrderUnderdampedFilter(TH1F*&, int, double, double, double, double = .3125, bool = false);
};
#endif