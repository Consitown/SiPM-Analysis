#ifndef _ReadRun
#define _ReadRun

// contact: christian.scharf@cern.ch
// some includes are probably not needed anymore

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TClonesArray.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TComplex.h>
#include <TObjString.h>
#include <TVirtualFFT.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
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
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TSpline.h> 
#include <TSpectrum.h>   // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h>      // root verbosity level
#include <TSystem.h>     // root verbosity level
#include <TLatex.h>      // root verbosity level

//#include <sys/resource.h>
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

using namespace std;

class ReadRun/* : public TObject*/ {
private:
	/// @brief Stores data
	TClonesArray* rundata;

	/// @brief Collects sums of all waveforms for each channel
	double** amplValuessum;

	/// @brief Events will be stored here in the order they have been read
	vector<unsigned int> eventnr_storage;

	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrum_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrumPMT_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrumPMTthreshold_cnt;


#pragma pack(1) // padding suppression
	// struct copied from
	// WaveCatcher binary -> root converter
	// by manu chauveau@cenbg.in2p3.fr
	// see https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2

	struct event_data
	{
		int EventNumber;
		double EpochTime;
		unsigned int Year;
		unsigned int Month;
		unsigned int Day;
		unsigned int Hour;
		unsigned int Minute;
		unsigned int Second;
		unsigned int Millisecond;
		unsigned long long int TDCsamIndex;
		int nchannelstored;
	};

	struct channel_data_without_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		short waveform[1024];
	};

	struct channel_data_with_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		float	MeasuredBaseline;
		float	AmplitudeValue;
		float	ComputedCharge;
		float	RiseTimeInstant;
		float	FallTimeInstant;
		float	RawTriggerRate;
		short waveform[1024];
	};
#pragma pack() // padding suppression

public:

	// plots amplValuessum
	void PlotChannelSums(bool = true, bool = false, double = 0., double = 4., bool = false);

	// baseline correction (shifts all waveforms individually)
	void CorrectBaseline(float, float = -999);
	void CorrectBaseline_function(TH1F*, float, float, int);

	void CorrectBaselineMinSlopeRMS(int = 100, bool = false, double = 10, int = 0, int = 0, bool = false, bool = false, int = 8);

	void CorrectBaselineMin(int = 100, bool = false, double = 10, int = 0, int = 0, bool = false, int = 8);

	// get timing of peaks
	void GetTimingCFD(float = .3, float = 100, float = 140, double = 0, bool = false, bool = false, bool = false);
	void SkipEventsTimeDiffCut(int, int, double, double, bool = false);

	void FractionEventsAboveThreshold(float = 4, bool = true, bool = true, double = 0., double = 0., bool = false);

	// average all waveforms to simplify peak ID
	void SmoothAll(double = 5, bool = false);
	void DerivativeAll();

	// functions for charge spectrum
	int* GetIntWindow(TH1F*, float, float, float, float, int);
	float GetPeakIntegral(TH1F*, float, float, float, float, int = 0);
	void PrintChargeSpectrumWF(float, float, float = 0, float = 300, int = 1, float = 0., float = 0.);
	TH1F* ChargeSpectrum(int, float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	void PrintChargeSpectrum(float, float, float = 0, float = 300, float = -50, float = 600, int = 750, float = 0., float = 0., int = 8, int = 0);
	/// @brief Starting values of the fit parameters for PrintChargeSpectrum()
	vector<float> PrintChargeSpectrum_pars;
	void PrintChargeSpectrumPMT(float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	/// @brief Starting values of the fit parameters for PrintChargeSpectrumPMT()
	vector<float> PrintChargeSpectrumPMT_pars;
	void PrintChargeSpectrumPMTthreshold(float = 0, float = 0, float = 0, float = 300, int = 750, double = 4, bool = false);

	// SiPM specific
	void PrintDCR(float = 15, float = 85, float = 0, float = 300, double = 3); // based on PrintChargeSpectrumPMTthreshold()

	// functions for time distribution
	TH1F* TimeDist(int, float = 0, float = 300, float = 0, float = 300, int = 100, int = 0, float = .3);
	void PrintTimeDist(float = 0, float = 300, float = 0, float = 300, int = 100, int = 0, float = .3);
	TGraph2D* MaxDist(int, float = 0, float = 300);
	void PrintMaxDist(float = 0, float = 300);

	TH1F* His_GetTimingCFD(int, float, float, int = -999);
	void Print_GetTimingCFD(float = 100, float = 140, int = 0, int = -999, string = "S");
	TH1F* His_GetTimingCFD_diff(vector<int>, vector<int>, float, float, int = -999);
	void Print_GetTimingCFD_diff(vector<int>, vector<int>, float = 100, float = 140, int = 0, int = -999, float = -999, float = -999, string = "RS");

	// print FFT
	void PrintFFTWF(int = 1, float = 0., float = 0., int = 1);

	// helper functions
	stringstream list_files(const char*, const char*);	// find data files
	TH1F* Getwf(int, int, int = 1);						// channel, eventnr, color
	double* getx(double = 0.);							// x values
	double* gety(int, int);								// y values for waveform(ch, event)
	double* gety(TH1F*);								// y values for histogram
	double* gety(TH1F*, int, int);						// y values for dedicated y range of a histogram 
	int rcolor(unsigned int);							// useful root colors
	
	float LinearInterpolation(float, float, float, float, float); // linear interpolation

	int GetEventIndex(int);										// get index of a triggered event (finds the correct event if files are not read sequentially)
	void SplitCanvas(TCanvas*&);								// split canvas into pads to display all active channels on one canvas
	void Convolute(double*&, double*, double*, int, int);		// convolution for filtering waveforms
	void SmoothArray(double*&, int, double = 1., bool = false);	// filtering

	/// @brief Constructor of the class with arguments to filter noise events in the cosmics setup. Default values do nothing 
	ReadRun(double = 0, int = 1);

	void ReadFile(string, bool = false, int = 9, string = "out.root", bool = false); // file name, bool whether or not to change sign of PMT channels (channel number>8)

	virtual ~ReadRun();

	/// @brief Path to data
	/// 
	/// Can be used to save analysis results in the data folder
	string data_path;

	/// @brief Do analysis only for limited range of channels to reduce memory usage
	/// 
	/// For large datasets with many channels and many events \n
	/// Only read and analyze channels from ReadRun::start_read_at_channel to ReadRun::end_read_at_channel. \n
	/// The recorded channel with the lowest wavecatcher channel number is 0 (e.g. recorded channels 3 and 4, so start would be 0 and end 1). \n 
	/// If set to -1 (default) all channels will be read in one go. \n
	/// Else channels from "start_read_at_channel" to "end_read_at_channel" will be read. \n 
	/// If "end_read_at_channel" is not defined will only read channel specified in "start_read_at_channel".
	int start_read_at_channel = -1;	
	/// @brief See ReadRun::start_read_at_channel
	int end_read_at_channel = -1;

	/// @brief Number of triggered events in data
	int nevents;
	/// @brief Number of active channels in data
	int nchannels;
	/// @brief Total number of waveforms in data (nchannels*nacquisitions)
	int nwf;

	/// @brief ns per bin in data (has to be .3125 ns)
	float SP;
	/// @brief Conversion coefficient for wavecatcher
	/// 
	/// From https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2 
	double coef;
	/// @brief Number of bins (always 1024 samples per waveform)
	int binNumber;

	/// @brief Stores bin numbers where the sum of waveforms have their maximum
	/// 
	/// Can be used for fixed integration window relative to maximum of the sum of all waveforms per channel (ReadRun::amplValuessum)
	int* maxSumBin;

	/// @brief Stores the numbers of the active channels
	vector<int> active_channels;
	/// @brief Stores the numbers of the active channels which should be plotted
	///
	/// You can select the channels you want with plot_active_channels.push_back(channel_to_plot); to add them to the list. \n 
	/// If undefined all channels will be plotted.
	vector<int> plot_active_channels;

	/// @brief Stores the fit results of PrintChargeSpectrum() for all channels and all function calls in ascending order 
	vector<TFitResultPtr> fit_results; 

	/// @brief Stores the event numbers which should be skipped in the analysis
	/// 
	/// To identify events to be filtered use functions IntegralFilter(), SkipEventsPerChannel(), and SkipEventsTimeDiffCut().
	vector<bool> skip_event;

	/// @brief Special parameter for HU cosmics setup
	/// 
	/// Threshold (usually 4 mV) for PMT signal (hardcoded channel >8) to skip events where PMTs pick up radio frequency noise (NO BASELINE CORRECTION!).
	double skip_event_threshold;
	/// @brief Special parameter for HU cosmics setup
	///
	/// define how many PMT channels need to be above threshold to discard event (RF pick up should be seen by alls PMTs).
	int skip_event_threshold_nch; 

	void SkipEventsPerChannel(vector<double>, double = 0, double = 0, bool = false);  // in case you want to have indiviual thresholds in individual channels
	void IntegralFilter(vector<double>, vector<bool>, float, float, float = 50, float = 250, bool = false, bool = false); // Same as SkipEventsPerChannel() but filtering all events with integrals <(>) threshold
	void PrintSkippedEvents();

	/// @brief Stores baseline correction results for CorrectBaseline() and related functions
	vector<vector<float>> baseline_correction_result;

	/// @brief Store timing of peaks from GetTimingCFD()
	vector<vector<float>> timing_results;
	/// @brief Stores the fit results of Print_GetTimingCFD() for all channels
	vector<TFitResultPtr> timing_fit_results;

	/// @brief Stores results of analysis
	TFile* root_out;

	//other controls 

	/// @brief Set true for baseline correction during data reading
	/// Needs to be called before ReadFile()
	bool Using_BaselineCorrection_in_file_loop = false;
	/// @brief Start of time window for baseline correction when using ReadRun::Using_BaselineCorrection_in_file_loop
	float tCutg;
	/// @brief End of time window for baseline correction when using ReadRun::Using_BaselineCorrection_in_file_loop
	float tCutEndg;

	/// @brief Shift waveforms with CFD so that all events start at the same time
	/// Call after initializing class and before calling ReadFile().
	bool Shift_WFs_in_file_loop = false;
	/// @brief Constant fraction of maximum (between ~0.1 and 1) for ReadRun::Shift_WFs_in_file_loop
	float tWF_CF = 0.3;
	/// @brief Time bin all events will be shifted to for ReadRun::Shift_WFs_in_file_loop
	/// Needs to be 300<"tWF_CF_bin"<500 ("tWF_CF_bin"=375 means all peaks will be shifted to 375*.3125 ns=117.1875 ns)
	int tWF_CF_bin = 375;
	/// @brief Start of range of bins where the signal is expected for ReadRun::Shift_WFs_in_file_loop
	int tWF_CF_lo = 320;
	/// @brief End of range of bins where the signal is expected for ReadRun::Shift_WFs_in_file_loop
	int tWF_CF_hi = 500;

	ClassDef(ReadRun, 1)
};
#endif

class Fitf {
public:
	/// @brief Default fit function for SiPMs missing after-pulses and dark counts
	/// 
	/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
	/// See https://root.cern/manual/fitting/ for ROOT fitting.
	/// 
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// @return Function value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};

class Fitf_full {
public:
	/// @brief Default fit function for SiPMs with after-pulses missing dark counts
	/// 
	/// Still missing dark counts in integration window (3.3 in paper). \n 
	/// Seems to be working now, but please check for possible bugs. \n \n 
	/// 
	/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
	/// See https://root.cern/manual/fitting/ for ROOT fitting.
	/// 
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// 7 - alpha: after-pulsing probability \n 
	/// 8 - beta: the inverse of the exponential slope of the after-pulse PH distribution
	/// @return Function value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double alpha = p[7];
			double beta = p[8];

			//numbe of fired cells
			double k = static_cast<double>(kint);
			//pulse width
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			double gausnormsigmak = 1. / (sqrt(2. * TMath::Pi()) * sigmaK);
			//gauss peak
			double gauss = gausnormsigmak * TMath::Exp(-1. * TMath::Power(x[0] - (k * G + B), 2.) / (sigmaK * sigmaK * 2.));

			//generalized poisson envelope
			double gp = p[0] * mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += (gp * gauss);
			}
			else {
				//after-pulse + delayed cross talk
				double bk0 = TMath::Power(1. - alpha, k);
				double bk1 = TMath::Factorial(kint) / TMath::Factorial(kint - 1) * alpha * TMath::Power(1. - alpha, k - 1);
				double pk1 = TMath::Exp(-1. * (x[0] - (k * G + B)) / beta) * gausnormsigmak / beta * sigmaK * sqrt(TMath::Pi() / 2) * (TMath::Erf((x[0] - (k * G + B)) / (sqrt(2) * sigmaK)) + 1.);


				if (kint == 1) sum += gp * (bk0 * gauss + bk1 * pk1);
				else {
					double api2k = 0.;
					for (int ii = 2; ii < kint; ii++) {
						double iid = static_cast<double>(ii);
						double bkialpha = TMath::Factorial(kint) / (TMath::Factorial(ii) * TMath::Factorial(kint - ii)) * TMath::Power(alpha, iid) * TMath::Power((1. - alpha), (k - iid));
						double dpkidph = 0;
						if (x[0] > (k * G + B)) dpkidph = TMath::Power((x[0] - (k * G + B)), (iid - 1.)) / (TMath::Factorial(ii - 1) * TMath::Power(beta, iid)) * TMath::Exp(-1. * (x[0] - (k * G + B)) / beta);
						api2k += bkialpha * dpkidph;
					}
					sum += gp * (bk0 * gauss + bk1 * pk1 + api2k);
				}
			}
		}
		return sum;
	};
};

class Fitf_biased {
public:
	/// @brief Fit function for SiPMs missing after-pulses and dark counts but with biased pedestal
	/// 
	/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
	/// See https://root.cern/manual/fitting/ for ROOT fitting.
	/// 
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Virtual pedestal shift of pe peaks \n 
	/// 7 - Pedestal scaling for biased pedestal \n 
	/// 8 - Position of biased pedestal
	/// @return Function value 
	/// @return Func value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double a_ped = p[7];
			double x_ped = p[8];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += a_ped * p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - x_ped) / sqrt(2) / sigmaK), 2));
			}
			else {
				sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
			}
		}
		return sum;
	};
};

class Fitf_PMT {
public:
	/// @brief Gauss-Poisson distribution for fit of PMT charge spectra
	/// 
	/// See https://doi.org/10.1016/0168-9002(94)90183-X 
	/// 
	/// @param x 
	/// @param p 
	/// 0 - A:		normalization to number of events in fit region \n 
	/// 1 - w:		probability for type II BG \n 
	/// 2 - alpha:	coefficient of exponential decrease of typ II BG \n 
	/// 3 - sigma0:	sigma of pedestal \n 
	/// 4 - Q0:		position of pedestal \n 
	/// 5 - mu:		mean number of PE \n 
	/// 6 - sigma1:	width of 1 PE peak \n 
	/// 7 - Q1:		position of 1 PE peak
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;
		int kmax = static_cast<int>(ceil(p[5])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else if (x[0] != p[4]) gn_term *= 0.; // delta function

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

class Fitf_PMT_pedestal {
public:
	/// @brief Gauss-Poisson distribution for fit of PMT charge spectra
	/// 
	/// With biased pedestal peak. \n
	/// See https://doi.org/10.1016/0168-9002(94)90183-X 
	/// 
	/// @param x 
	/// @param p
	/// 0 - A:		normalization to number of events in fit region \n 
	/// 1 - w:		probability for type II BG \n 
	/// 2 - alpha:	coefficient of exponential decrease of typ II BG \n 
	/// 3 - sigma0:	sigma of pedestal \n 
	/// 4 - Q0:		position of pedestal \n 
	/// 5 - mu:		mean number of PE \n 
	/// 6 - sigma1:	width of 1 PE peak \n 
	/// 7 - Q1:		position of 1 PE peak \n 
	/// 8 - norm0:	norm of 0 PE peak 
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;
		int kmax = static_cast<int>(ceil(p[5])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			else normgn_term *= p[8] / (p[3] * TMath::Sqrt(2. * TMath::Pi()));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else  gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - p[4], 2.) / (2. * p[3] * p[3]));

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

class Fitf_PMT_ideal {
public:
	// Gauss-Poisson
	// https://doi.org/10.1016/0168-9002(94)90183-X 

	/// @brief Ideal PMT charge spectrum
	/// 
	/// Gives very good fit but will not describe pedestal well. \n 
	/// See https://doi.org/10.1016/0168-9002(94)90183-X 
	/// 
	/// @param x 
	/// @param p 
	/// 0 - A_s:		norm. of PE spectrum \n 
	/// 1 - mu:		mean number of PE \n 
	/// 2 - sigma1:	width of 1st PE peak \n 
	/// 3 - Q1:		position (gain*e) of 1st peak
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;

		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 1; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[1], k) * TMath::Exp(-1. * p[1]) / TMath::Factorial(kint);

			double norm = 1. / (p[2] * TMath::Sqrt(2. * TMath::Pi() * k));

			double gauss = TMath::Exp(-1. * TMath::Power(x[0] - k * p[3], 2) / (2. * k * p[2] * p[2]));

			pmt_charge_spectrum += p[0] * poiss * norm * gauss;
		}
		return pmt_charge_spectrum;
	};
};

class Fitf_langaus {
public:
	/// @brief Landau-Gauss-convolution
	/// 
	/// Used for large photon yields. \n 
	/// From https://root.cern.ch/doc/master/langaus_8C.html => \n
	/// In the Landau distribution (represented by the CERNLIB approximation), \n 
	/// the maximum is located at x=-0.22278298 with the location parameter=0. \n 
	/// This shift is corrected within this function, so that the actual \n 
	/// maximum is identical to the MP parameter. 
	/// 
	/// @param x 
	/// @param par
	/// par[0]=Width (scale) parameter of Landau density \n 
	/// par[1]=Most Probable (MP, location) parameter of Landau density \n 
	/// par[2]=Total area (integral -inf to inf, normalization constant) \n 
	/// par[3]=Width (sigma) of convoluted Gaussian function
	/// @return Func value
	double operator() (double* x, double* par) {
		// Numeric constants
		Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
		Double_t mpshift = -0.22278298;       // Landau maximum location

		// Control constants
		Double_t np = 100.0;      // number of convolution steps
		Double_t sc = 5.0;      // convolution extends to +-sc Gaussian sigmas

		// Variables
		Double_t xx;
		Double_t mpc;
		Double_t fland;
		Double_t sum = 0.0;
		Double_t xlow, xupp;
		Double_t step;
		Double_t i;


		// MP shift correction
		mpc = par[1] - mpshift * par[0];

		// Range of convolution integral
		xlow = x[0] - sc * par[3];
		xupp = x[0] + sc * par[3];

		step = (xupp - xlow) / np;

		// Convolution integral of Landau and Gaussian by sum
		for (i = 1.0; i <= np / 2; i++) {
			xx = xlow + (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);

			xx = xupp - (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);
		}

		return (par[2] * step * sum * invsq2pi / par[3]);
	};
};

class Fitf_plus_DC {
public:
	/// @brief Fit function for SiPMs missing after-pulses and dark counts but including additional dark count spectrum
	/// 
	/// Sum of two spectra for event spectrum + dark count (background trigger) spectrum. 
	/// To be used if the trigger and filter selections are not perfect and there are still empty events.
	/// 
	/// @param x 
	/// @param p
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// 7 - mu_dk: mu for dark count rate in dark events (dark events = noise/background triggers without photons in SiPMs) \n 
	/// 8 - N0_dk: fraction (dark events)/(total number of events) 
	/// @return 
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double mu_dk = p[7];
			double N0_dk = p[8];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = ((1. - N0_dk) * (mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda))) + N0_dk * (mu_dk * TMath::Power((mu_dk + k * lambda), k - 1) * TMath::Exp(-(mu_dk + k * lambda)))) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};