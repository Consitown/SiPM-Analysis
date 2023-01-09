#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>

// you'll need to compile the analysis (ReadRun.cc) first! --> use command 'make -f makefile' in the directory where ReadRun.cc and ReadRun.h are

using namespace std;

void read_cfd_PS(int which) // main
{
	string path;

	// edit for your work-directory
	path = "/mnt/d/Work_SHK_Bachelor/analysis_programm/measurements/";
	int run = 0;

	// WARNING: the name of the folder where the .bin-files are stored must be the same as the name of the bin files
	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "5time_testrun_PS_better/"; // PS run with the right trigger settings for CFD
		run = 5;
		break;
	}
	case(1): {
		path += "25_ortho_PS_test/"; // test run with orthogonal PS
		run = 25;
		break;
	}
	case(2): {
		path += "26_normal_PS_coinc_with_ortho_PS/"; // test run with orthogonal PS and normal PS; combined: C0 position; also blanket around lower PS, cable swicth of 11/13
		run = 26;
		break;
	}
	case(3): {
		path += "27_normal_PS_coinc_with_ortho_PS_C1/"; //as run 26, but orthogonal PS at C1 position
		run = 27;
		break;
	}
	case(4): {
		path += "28_normal_PS_coinc_with_ortho_PS_C2/"; //as run 26, but orthogonal PS at C2 position
		run = 28;
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
	// read data; mymeas.ReadFile(path, true, 0, path + "/cfd_results.root") for an explizit output file
	mymeas.ReadFile(path, true, 0, path + "/cfd_results.root");

	// only plot channels specified below. Leaving it empty or commenting it will plot all channels
	//mymeas.plot_active_channels={};

	//use SmoothAll() for general smoothing of waveforms --> find good parameters; then use the internal smoothig of GetTimingCFD
	//mymeas.SmoothAll(3, false); //Syntax: ...(double sigma, bool doconv) ; doconv - If false use running average (default). If true use gaussian smoothing (slower).

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaseline(0., 50.); // use mean from 0 ns to 50 ns

	// Syntax: ...(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	//mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 10, 0, false, false, 9);

	// Syntax: ...(vector<double> thresholds, double rangestart, double rangeend, bool verbose)
	//std::vector<double> thresholds = {0, 0, -7, -7, -7, -7}; //for identifying out burst events (high frequency oscillations), if run5: erase the 0's
	vector<double> thresholds2 = {4, 0}; //skip all events where ch8 fires
	mymeas.SkipEventsPerChannel(thresholds2, 110, 150, false);
	vector<bool> first_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = false; //reset the vector
	vector<double> thresholds3 = {0, 4}; //skip all events where ch9 fires
	mymeas.SkipEventsPerChannel(thresholds3, 110, 150, false);
	vector<bool> second_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = second_ch_skip[i] && first_ch_skip[i]; //skip all events where ch8 AND ch9 fires --> skip all orthogonal events
	int counter = 0;
	for (int i = 0; i < mymeas.skip_event.size(); i++) {if (mymeas.skip_event[i]) ++counter;} //a bit debugging
	cout << "Number of orthogonal events: " << counter << endl;
	mymeas.skip_event.flip();

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Syntax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, bool doconv, bool use_spline) --> fifth argument is selecting inverse/normal cfd: true results in normal cfd
	mymeas.GetTimingCFD(cfd_x, 110, 150, 3, true, false, false); // this creates the timing_results matrix

	// Positioncuts (trial)
	//mymeas.SkipEventsTimeDiffCut(10, 11, -2.5, 2.5, false);
	//mymeas.SkipEventsTimeDiffCut(12, 13, -2.5, 2.5, false);

	// Syntax: ...(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int do_fit, int nbins, float fitrangestart, float fitrangeend, string fitoption)
	//from entering lab: upper right: 10, upper left: 11, lower right: 12, lower left: 13; channel 8 is upper orthogonal PMT, 9 is lower
	mymeas.Print_GetTimingCFD_diff({10,11}, {12,13}, -15, 15, 0, 200, -8, 8, "RS");

	// Syntax: ... (float rangestart, float rangeend, int do_fit, int nbins, string fitoption)
	//mymeas.Print_GetTimingCFD(110,140,1,200,"S"); //for channelwise cfd

	// prints some stats for events above a threshold into the terminal to identify interesting events
	// Syntax: ...(float threshold, bool max, bool greater, double from, double to, bool verbose)
	//mymeas.FractionEventsAboveThreshold(5, true, true, 200, 250, false);

	//**********//
	// Plotting
	//**********//

	// plot sums of all events per channel --> get parameters
	// Syntax: ...(bool doaverage, bool normalize, double shift, double sigma, bool doconv)
	//mymeas.PlotChannelSums(false);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 3.;	// lower integration window in ns rel. to max
	float intwindowplus = 5.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	
	// old timing histogramm (should not use for cfd)
	//mymeas.PrintTimeDist(110, 140, 105, 145, 200, 1, cfd_x);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -10, 350, 300, 0, 0, 0, 0);

	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins)
	// PrintChargeSpectrumPMT will apply a fit automatically from rangestart to rangeend (which are also the boundaries for the plot)
	//mymeas.PrintChargeSpectrumPMT(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 5, 300, 202);

	// plot waveforms of individual events
	//int event1 = 2257;
	//int event2 = 79;

	//plot range for individual waveforms
	double ymin = -5;
	double ymax = 25;

	// plot waveforms for certain events with integration window
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event1, ymin, ymax);
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event2, ymin, ymax);

	// suppress graphic output
	gROOT->SetBatch(kTRUE); // TRUE enables batch-mode --> disables graphic output (all prints before this will still be shown)

	// plot individual waveforms for some events for debugging
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 10)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	}
	gROOT->SetBatch(kFALSE);
}