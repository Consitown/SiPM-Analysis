#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>

// you'll need to compile the analysis (ReadRun.cc) first! --> use command 'make -f makefile' in the directory where ReadRun.cc and ReadRun.h are

using namespace std;

void read_everyother_group(int which) // main
{
	string path;

	// edit for your work-directory
	path = "/users/eel/lab/SHiP/ConstantinWCData/";
	int run = 0;

	// WARNING: the name of the folder where the .bin-files are stored must be the same as the name of the bin files
	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(6): {
		path += "37_nor_PS_coinc_w_ort_PS+SiPMs_C0_120k/";
		run = 37;
		break;
	}
	case(7): {
		path += "38_nor_PS_coinc_w_ort_PS+SiPMs_R0_120k/";
		run = 38;
		break;
	}
	case(8): {
		path += "39_nor_PS_coinc_w_ort_PS+SiPMs_L0_120k/";
		run = 39;
		break;
	}
	case(9): {
		path += "Andrea_C0/"; // Data from Andrea at the C0 position for comparison
		break;
	}
	case(10): {
		path += "Andrea_R0/"; // Data from Andrea at the R0 position for comparison
		break;
	}
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);
	mymeas.Shift_WFs_in_file_loop = true;
	mymeas.tWF_CF = 0.5;
	mymeas.tWF_CF_lo = 230;
	mymeas.tWF_CF_hi = 400;

	// Syntax:...(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debu)
	// read data; mymeas.ReadFile(path, true, 0, path + "cal_results.root") for an explizit output file
	mymeas.ReadFile(path, true, 8, path+"everyother_results.root");

	// only plot channels specified below. Leaving it empty will plot all channels
	mymeas.plot_active_channels={0,1,2,3,4,5,6,7};
	//mymeas.plot_active_channels={8,9,10,11,12,13};

	//use SmoothAll() for general smoothing of waveforms --> find good parameters; then use the internal smoothig of GetTimingCFD
	//mymeas.SmoothAll(3, false); //Syntax: ...(double sigma, bool doconv) ; doconv - If false use running average (default). If true use gaussian smoothing (slower).

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.CorrectBaseline(0., 50.); // use mean from 0 ns to 50 ns
	// Syntax: ...(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, int smooth_method, int skip_channel)
	mymeas.CorrectBaselineMin(50, false, 0, 352, 192); // use minimal mean method

	// Syntax: ...(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	//mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 10, 0, false, false, 9);
/*
	// filter out events where all PMT's have fired. small PS are channel 12 and 13
	// Syntax: ...(vector<double> thresholds, double rangestart, double rangeend, bool verbose)
	vector<double> thresholds2 = {0, 0, 0, 0, 0, 0, 0, 0, 4, 0}; //skip all events where ch8 fires
	mymeas.SkipEventsPerChannel(thresholds2, 110, 150, false);
	vector<bool> first_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = false; //reset the vector
	vector<double> thresholds3 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 4}; //skip all events where ch9 fires
	mymeas.SkipEventsPerChannel(thresholds3, 110, 150, false);
	vector<bool> second_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = second_ch_skip[i] && first_ch_skip[i]; //skip all events where ch8 AND ch9 fires --> skip all orthogonal events
	int counter = 0;
	for (int i = 0; i < mymeas.skip_event.size(); i++) {if (mymeas.skip_event[i]) ++counter;} //a bit debugging; necessary for event number computation
	cout << "Number of orthogonal events: " << counter << endl;
	mymeas.skip_event.flip(); //all non ortho events get skipped --> only ortho events get plotted
*/
	//skip weird events: -5 mV for identifying out burst/weird events (high frequency oscillations) in large PS (-7 only if signals inverted)
	vector<double> thresholds_weird = {0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5, -5, -5};
	vector<double> thresholds_weird_andrea = {0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5};
	mymeas.SkipEventsPerChannel(thresholds_weird, 100, 200, false);

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Syntax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, bool doconv, bool use_spline) --> fifth argument is selecting inverse/normal cfd: true results in normal cfd
	mymeas.GetTimingCFD(cfd_x, 100, 150, 3, true, false, false); // this creates the timing_results matrix/vector

	// apply cut for time difference between two channels
	mymeas.SkipEventsTimeDiffCut(10, 11, 4, 100, false); //skip all events which time diff is not of interest to us
	mymeas.SkipEventsTimeDiffCut(12, 13, 5, 100, false); //skip all events which time diff is not of interest to us
	//mymeas.skip_event.flip();

	//**********//
	// Plotting + Fitting
	//**********//

	// plot sums of all events per channel --> get parameters
	// Syntax: ...(bool doaverage, bool normalize, double shift, double sigma, bool doconv)
	//mymeas.PlotChannelSums(false);

	// PMT's are at channels 10-13; 10 - right upper PMT, 11 - left upper PMT, 12 - right lower PMT, 13 - left lower PMT (perspective from door)
	// Syntax: ...(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int do_fit, int nbins, float fitrangestart, float fitrangeend, string fitoption, bool set_errors)
	//from entering lab: upper right: 10, upper left: 11, lower right: 12, lower left: 13; channel 8 is upper orthogonal PMT, 9 is lower
	// do_fit: 0 - no fit, 1 - gauss, 2 - gauss*exp decay; fitoptions: use "LRS" for log likelihood and "RS" for chi-squared
	//mymeas.Print_GetTimingCFD_diff({10}, {11}, -15, 15, 1, 200, -8, 8, "RS", true);
	//mymeas.Print_GetTimingCFD_diff({12}, {13}, -15, 15, 2, 200, -8, 8, "RS", true);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 20.;	// lower integration window in ns rel. to max
	float intwindowplus = 10.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	float plotrangestart = 0; 	// decide on the
	float plotrangeend = 6000;	// plotrange of x-axis
	float fitstart = 15;		// start of the fit (x-axis)
	float fitend = 400;		// end of the fit
	int channels_to_fit = 8; 	// numbers of channels to apply the fit to (counts like this: i=0;i<channels_to_fit;i++); thats why one should use the first channels for data and the later channels for triggering
	int which_fit = 0;			// decide on which fit-function you want; options: default (no value besides 1-6) - default SiPM fit function
								// 1 - landau gauss convolution for large number of photons
								// 2 - biased: if pedestal is biased because of peak finder algorithm
								// 3 - SiPM fit function with exponential delayed afterpulsing
								// 4 - ideal PMT fit function (5 is similar)
								// 6 - PMT fit function with biased pedestal
	
	// Plotting the phi_ew spectrum
	// Syntax: ...(vector<int> phi_chx, vector<float> ly_C0, vector<int> SiPMchannels, float windowmin, float windowmax, float maxfrom, float maxto, int nbins, bool corr, bool triple_gauss)
	vector<int> phi_chx = {225-4, 270-4, 315-4, 0-4, 45-4, 90-4, 135-4, 180-4}; //ordered from channel 0 to channel 7; my channel alignment was a bit different from Alex's
	//vector<int> phi_chx_even = {225, 315, 45, 135}; vector<int> phi_chx_odd = {270, 0, 90, 180};
	//vector<float> ly_C0 = {106, 59.01, 105, 59.27, 107.6, 67.13, 111.2, 58.48}; //mean lightyields from PrintChargeSpectrum for run 37 (C0), selection with small PS
	//vector<float> ly_C0 = {124.4, 63.9418, 120.342, 59.6804, 121.875, 69.4333, 128.262, 62.1027}; //selection with 2ns window (time_cut_average_2ns_window.png)
	//vector<float> ly_C0 = {150.545, 77.3498, 147.05, 73.9891, 150.324, 85.745, 156.004, 74.9042}; //selection with 1ns window (time_cut_average_1ns_window.png)
	//vector<float> ly_C0 = {1107.14, 532.202, 1086.62, 521.439, 1120.01, 608.955, 1159.26, 520.146}; //selection with 1ns window (time_cut_average_1ns_window.png) + larger int-window (larger_intw_catch_prim_gamma.png)
	//vector<float> ly_C0 = {1071.33, 581.204, 1064.49, 584.104, 1097.8, 677.876, 1136.27, 576.991}; //selection with small PS + larger int-window (larger_intw_catch_prim_gamma.png)
	vector<float> ly_C0 = {1361.92, 768.23, 1354.06, 765.312, 1385.92, 873.016, 1433.86, 774.65}; //selection with small PS + larger int-window2 (larger_intw_catch_prim_gamma2.png)
	//vector<float> ly_C0 = {1399.68, 688.207, 1359.28, 667.615, 1411.16, 776.851, 1461.85, 671.025}; //selection with 1ns window + larger int-window (time_cut_run37_smallPS.png)
	//vector<float> ly_C0_even = {1071.33, 1064.49, 1097.8, 1136.27}; vector<float> ly_C0_odd = {581.204, 584.104, 677.876, 576.991};
	vector<int> SiPMchannels = {0, 1, 2, 3, 4, 5, 6, 7};
	//vector<int> SiPMchannels = {1, 3, 5, 7};
	//vector<int> phi_chx_andrea = {0, 315, 270, 225, 45, 90, 135, 180};
	//vector<float> ly_C0_andrea = {120.2, 126.9, 111.4, 108.6, 110.9, 113.8, 119.6, 108};
	mymeas.Print_Phi_ew(phi_chx, ly_C0, SiPMchannels, intwindowminus, intwindowplus, findmaxfrom, findmaxto, 200, true, true);
	mymeas.Print_Phi_ew(phi_chx, ly_C0, SiPMchannels, intwindowminus, intwindowplus, findmaxfrom, findmaxto, 200, false, true);

	// plot results between t=110 ns and t=140 ns and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1, 200);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// for getting average lightyield, do which_fit=0 and look at the histogramms in the .root-file and manually extract the means
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotrangestart, plotrangeend, 200, fitstart, fitend, channels_to_fit, which_fit);
	//for(int i = 0; i < mymeas.mean_integral.size(); i++) cout << mymeas.mean_integral[i] << endl;

	// suppress graphic output
	gROOT->SetBatch(kTRUE); // TRUE enables batch-mode --> disables graphic output (all prints before this will still be shown)

	// plot waveforms of individual events
	//int event1 = 68;
	//int event2 = 79;

	//plot range for individual waveforms
	double ymin = -5;
	double ymax = 25;

	// plot waveforms for certain events with integration window
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event1, ymin, ymax);
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event2, ymin, ymax);

	// plot individual waveforms for some events for debugging
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 10)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	}
	gROOT->SetBatch(kFALSE);
}