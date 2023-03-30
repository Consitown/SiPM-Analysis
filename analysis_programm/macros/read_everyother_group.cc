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
	path = "/mnt/d/Work_SHK_Bachelor/analysis_programm/measurements/";
	int run = 0;

	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "19cosmics_pcbj_everyothergroup_vb42_PSleft/"; // only every other group of the SiPM's
		run = 19;
        break;
	}
	case(1): {
		path += "22cosmics_pcbj_everyothergroup_vb42_PSmid/"; // has optical coupling; apparently, channel 11 is empty
		run = 22;
		break;
	}
	case(2): {
		path += "23cosmics_pcbj_everyothergroup_vb42_PSright/"; // PS positions are left (closer to door), mid, right (farther to door); apparently, channel 11 is empty
		run = 23;
		break;
	}
	case(3): {
		path += "24cosmics_pcbj_everyothergroup_vb42_PSleft_new/"; // PS positions are left (closer to door), mid, right (farther to door); apparently, channel 11 is empty
		run = 24;
		break;
	}
	case(4): {
		path += "30_nor_PS_coinc_w_ort_PS+SiPMs_C0/"; // PS positions are left (closer to door), mid, right (farther to door); apparently, channel 11 is empty
		run = 30;
		break;
	}
	case(5): {
		path += "31_all_PS_coinc+SiPMs_C0/"; // PS positions are left (closer to door), mid, right (farther to door); apparently, channel 11 is empty
		run = 31;
		break;
	}
	case(6): {
		path += "37_nor_PS_coinc_w_ort_PS+SiPMs_C0_120k/";
		run = 37;
		break;
	}
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);

	// Syntax:...(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debu)
	// read data; mymeas.ReadFile(path, true, 0, path + "cal_results.root") for an explizit output file
	mymeas.ReadFile(path, true, 8, "everyother_results.root");

	// only plot channels specified below. Leaving it empty will plot all channels
	//mymeas.plot_active_channels={0,1,2,3,4,5,6,7};
	//mymeas.plot_active_channels={8,9,10,11,12,13};

	//use SmoothAll() for general smoothing of waveforms --> find good parameters; then use the internal smoothig of GetTimingCFD
	//mymeas.SmoothAll(3, false); //Syntax: ...(double sigma, bool doconv) ; doconv - If false use running average (default). If true use gaussian smoothing (slower).

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaseline(0., 50.); // use mean from 0 ns to 50 ns

	// Syntax: ...(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	//mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 10, 0, false, false, 9);

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

	//skip weird events: -5 mV for identifying out burst/weird events (high frequency oscillations) in large PS (-7 only if signals inverted)
	vector<double> thresholds_weird = {0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5, -5, -5};
	mymeas.SkipEventsPerChannel(thresholds_weird, 100, 200, false);
/*
	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Syntax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, bool doconv, bool use_spline) --> fifth argument is selecting inverse/normal cfd: true results in normal cfd
	mymeas.GetTimingCFD(cfd_x, 110, 150, 3, true, false, false); // this creates the timing_results matrix

	// apply cut for time difference between two channels
	mymeas.SkipEventsTimeDiffCut(10, 11, -1, 1, false); //skip all events which time diff is not of interest to us
	mymeas.SkipEventsTimeDiffCut(12, 13, -0.5, 1.5, false); //skip all events which time diff is not of interest to us
	//mymeas.skip_event.flip();
*/
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
	//mymeas.Print_GetTimingCFD_diff({10}, {11}, -15, 15, 2, 200, -8, 8, "RS", true);
	//mymeas.Print_GetTimingCFD_diff({12}, {13}, -15, 15, 2, 200, -8, 8, "RS", true);

	// print events above a threshold to identify interesting events
	// mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 1.;	// lower integration window in ns rel. to max
	float intwindowplus = 1.;	// upper integration window in ns rel. to max
	float findmaxfrom = 110.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	float plotrangestart = 0; 	// decide on the
	float plotrangeend = 400;	// plotrange of x-axis
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
	vector<int> phi_chx = {225, 270, 315, 0, 45, 90, 135, 180}; //ordered from channel 0 to channel 7; my channel alignment was a bit different from Alex's
	//vector<int> phi_chx_even = {225, 315, 45, 135}; vector<int> phi_chx_odd = {270, 0, 90, 180};
	vector<float> ly_C0 = {106, 59.01, 105, 59.27, 107.6, 67.13, 111.2, 58.48}; //mean lightyields from PrintChargeSpectrum for run 37 (C0)
	//vector<float> ly_C0_even = {106, 105, 107.6, 111.2}; vector<float> ly_C0_odd = {59.01, 59.27, 67.13, 58.48};
	vector<int> SiPMchannels = {0, 1, 2, 3, 4, 5, 6, 7};
	//vector<int> SiPMchannels = {0, 2, 4, 6};
	mymeas.Print_Phi_ew(phi_chx, ly_C0, SiPMchannels, intwindowminus, intwindowplus, findmaxfrom, findmaxto, 200, true, true);

	// old timing histogramm (should not use for cfd, so cfd_x should be 1)
	//mymeas.PrintTimeDist(110, 140, 105, 145, 200, 1, cfd_x);

	// plot results between t=110 ns and t=140 ns and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// for getting average lightyield, do which_fit=0 and look at the histogramms in the .root-file and manually extract the means
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotrangestart, plotrangeend, 200, fitstart, fitend, channels_to_fit, which_fit);
	//for(int i = 0; i < mymeas.mean_integral.size(), i++) cout << mymeas.mean_integral[i] << endl;
	
	// suppress graphic output
	//gROOT->SetBatch(kTRUE); // TRUE enables batch-mode --> disables graphic output (all prints before this will still be shown)

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
	//gROOT->SetBatch(kFALSE);
}