#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <iomanip> //for precision couting

// you'll need to compile the analysis (ReadRun.cc) first! --> use command 'make -f makefile' in the directory where ReadRun.cc and ReadRun.h are

using namespace std;

void read_everyother_group_server(int which) // main
{
	string path;

	// edit for your work-directory
	path = "/users/eel/lab/SHiP/ConstantinWCData/";
	int run = 0;

	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "26_normal_PS_coinc_with_ortho_PS/";
		run = 26;
		break;
	}
	case(1): {
		path += "29_normal_PS_coinc_with_ortho_PS_C1_20k/";
		run = 29;
		break;
	}
	case(2): {
		path += "28_normal_PS_coinc_with_ortho_PS_C2/";
		run = 28;
		break;
	}
	case(3): {
		path += "32_nor_PS_coinc_w_ort_PS_at_7toC1/";
		run = 32;
		break;
	}
	case(4): {
		path += "33_nor_PS_coinc_w_ort_PS_at_7toC2/";
		run = 33;
		break;
	}
	case(5): {
		path += "34_nor_PS_coinc_w_ort_PS_at_205toC2/";
		run = 34;
		break;
	}
	case(13): {
		path += "35_nor_PS_coinc_w_ort_PS_at_1925toleft/";
		run = 35;
		break;
	}
	case(6): {
		path += "37_nor_PS_coinc_w_ort_PS+SiPMs_C0_120k/"; //used for C1/C2 and the C0-couping-correction 
		run = 37;
		break;
	}
	case(7): {
		path += "38_nor_PS_coinc_w_ort_PS+SiPMs_R0_120k/"; //used for R0
		run = 38;
		break;
	}
	case(8): {
		path += "39_nor_PS_coinc_w_ort_PS+SiPMs_L0_120k/"; //used for L0/L1/L2
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
	case(11): {
		path += "Andrea_L0/"; // Data from Andrea at the R0 position for comparison
		break;
	}
	case(12): {
		path += "40_nor_PS_coinc_w_ort_PS+SiPMs_R2_120k/"; //used for R1/R2 and the event selection comparison: short PS vs. time cuts
		break;
	}
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);
/*
	mymeas.Shift_WFs_in_file_loop = true;
	mymeas.tWF_CF = 0.5;
	mymeas.tWF_CF_lo = 230;
	mymeas.tWF_CF_hi = 480;
*/
	// Syntax:...(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debu)
	// read data; mymeas.ReadFile(path, true, 0, path + "cal_results.root") for an explizit output file
	mymeas.ReadFile(path, true, 8, path+"everyother_results.root");

	//mymeas.discard_original_eventnr = true; //false is default

	// only plot channels specified below. Leaving it empty will plot all channels
	//mymeas.plot_active_channels={0,1,2,3,4,5,6,7};
	//mymeas.plot_active_channels={10,11,12,13};
	//mymeas.plot_active_channels={0};

	//use SmoothAll() for general smoothing of waveforms --> find good parameters; then use the internal smoothig of GetTimingCFD
	//mymeas.SmoothAll(3, false); //Syntax: ...(double sigma, bool doconv) ; doconv - If false use running average (default). If true use gaussian smoothing (slower).

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.CorrectBaseline(0., 50.); // use mean from 0 ns to 50 ns
	//using minimal method:
	//Syntax: ...(vector<float> window, double sigma, int smooth_method, int increment)
	vector<float> window = {15.625, 60, 110};
	mymeas.CorrectBaselineMin(window, 0, 2, 3); // use minimal mean method

	// Syntax: ...(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	//mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 10, 0, false, false, 9);
/*
	// filter out events where all PMT's have fired. small PS are channel 12 and 13.
	// This is the didgital trigger for the small PS
	// Syntax: ...(vector<double> thresholds, double rangestart, double rangeend, bool verbose)
	vector<double> thresholds2 = {0, 0, 0, 0, 0, 0, 0, 0, 4, 0}; //skip all events where ch8 fires
	vector<double> thresholds2_PS_only = {4, 0};
	mymeas.SkipEventsPerChannel(thresholds2, 110, 150, false);
	vector<bool> first_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = false; //mymeas.UnskipAll(); //reset the vector
	vector<double> thresholds3 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 4}; //skip all events where ch9 fires
	vector<double> thresholds3_PS_only = {0, 4};
	mymeas.SkipEventsPerChannel(thresholds3, 110, 150, false);
	vector<bool> second_ch_skip = mymeas.skip_event; //save the info
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = second_ch_skip[i] && first_ch_skip[i]; //skip all events where ch8 AND ch9 fires --> skip all orthogonal events
	int counter = 0;
	for (int i = 0; i < mymeas.skip_event.size(); i++) {if (mymeas.skip_event[i]) ++counter;} //a bit debugging; necessary for event number computation
	cout << "Number of orthogonal events: " << counter << endl;
	mymeas.skip_event.flip(); //all non ortho events get skipped --> only ortho events get plotted
*/
/*
	// estimating events with low amplitude --> these events are hard for CFD algorithm
	float* maxima_first = mymeas.ChargeList(mymeas.GetChannelIndex(10), 0, 0, 90, 150, true);
	float* maxima_second = mymeas.ChargeList(mymeas.GetChannelIndex(11), 0, 0, 90, 150, true);
	float* maxima_third = mymeas.ChargeList(mymeas.GetChannelIndex(12), 0, 0, 90, 150, true);
	float* maxima_fourth = mymeas.ChargeList(mymeas.GetChannelIndex(13), 0, 0, 90, 150, true);
	int counter3 = 0;

	int max_threshold = 5;
	for (int j = 0; j < mymeas.nevents; j++) {
		if ((maxima_first[j] < max_threshold) || (maxima_second[j] < max_threshold) || (maxima_third[j] < max_threshold) || (maxima_fourth[j] < max_threshold)) {
			mymeas.skip_event[j] = true;
			counter3++;
		}
	}

	//comment this out, if you just want to skip low amplitude events
	//vector<bool> maxima_selection = mymeas.skip_event; // save the info
	//mymeas.UnskipAll(); // reset the vector

	cout << "\n\n" << counter3 << " events with at least 1 channel with less than " << max_threshold << " mV amplitude" << endl;
*/
	//skip weird events: -5 mV for identifying out burst/weird events (high frequency oscillations) in large PS (-7 only if signals inverted)
	vector<double> thresholds_weird = {0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5, -5, -5};
	vector<double> thresholds_weird_andrea = {0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5};
	vector<double> thresholds_PS_only = {0, 0, -5, -5, -5, -5};
	mymeas.SkipEventsPerChannel(thresholds_weird_andrea, 100, 200, false);
	//mymeas.skip_event.flip();

/*	//this is for plotting low amplitude events
	//comment this out, if you just want to skip low amplitude events
	for (int i = 0; i < mymeas.skip_event.size(); i++) mymeas.skip_event[i] = maxima_selection[i] && !mymeas.skip_event[i];
	mymeas.skip_event.flip();
*/
	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .5;
	// Syntax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, int smooth_method, bool use_spline) --> fifth argument is selecting inverse/normal cfd: true results in normal cfd
	// for high slope of the flank of the waveforms peak, suggested settings: smooth_method = 2, sigma = 0.3 (this is the sigma of the gauss, that will be used for smoothing --> small number)
	mymeas.GetTimingCFD(cfd_x, 90, 150, 1500, 0, true, 0, false); // this creates the timing_results matrix/vector
/*
	//skipping events with failed interpolation
	int first_channel = mymeas.GetChannelIndex(10);
	int second_channel = mymeas.GetChannelIndex(11);
	int third_channel = mymeas.GetChannelIndex(12);
	int fourth_channel = mymeas.GetChannelIndex(13);
	int counter2 = 0;

	for (int j = 0; j < mymeas.nwf; j += mymeas.nchannels) {
		if (!mymeas.skip_event[mymeas.GetCurrentEvent(j)]) {
			if (mymeas.timing_results[j + first_channel][7] == 1 || mymeas.timing_results[j + second_channel][7] == 1 || mymeas.timing_results[j + third_channel][7] == 1 || mymeas.timing_results[j + fourth_channel][7] == 1) {
				mymeas.skip_event[mymeas.GetCurrentEvent(j)] = true;
				counter2++;
			}
		}
	}
	cout << "\n" << counter2 << " events skipped due to failed interpolation in trigger channels" << endl;
*/
/*
	//writing out cfd-times to .txt files
	//open the .txt-files; edit for your directory
	ofstream myfile11; ofstream myfile13;
	string directory1 = path + "run" + to_string(run) + "_11-10_improv_cfd_times_ortho_only.txt";
	string directory2 = path + "run" + to_string(run) + "_13-12_improv_cfd_times_ortho_only.txt";
	myfile11.open(directory1.c_str());
	myfile13.open(directory2.c_str());
	myfile11 << "Delta t of cfd-times for channel 11-10\n"; myfile13 << "Delta t of cfd-times for channel 13-12\n";
	for (int i = 0; i < mymeas.nevents; i++) {
		if (!mymeas.skip_event[i]) {
			myfile11 << mymeas.timing_results[i*mymeas.nchannels+second_channel][1]-mymeas.timing_results[i*mymeas.nchannels+first_channel][1] << "\n";
			myfile13 << mymeas.timing_results[i*mymeas.nchannels+fourth_channel][1]-mymeas.timing_results[i*mymeas.nchannels+third_channel][1] << "\n";
		}
	}
*/
	// apply cut for time difference between two channels
	// for C1,L1,R1: 11-10 - 2.5 , 13-12 - 3 ; C2,L2,R2: 11-10 - 1.9 , 13-12 - 1.9 ;  andrea C0,L0,R0: 11-10 - [-0.2836,0.7164] , 13-12 - [-0.0892,0.9108]
	mymeas.SkipEventsTimeDiffCut(10, 11, 2.5, 100, false); //skip all events which time diff is not of interest to us
	mymeas.SkipEventsTimeDiffCut(12, 13, 3, 100, false); //skip all events which time diff is not of interest to us
	//mymeas.SkipEventsTimeDiffCut(10, 11, -100, -1.9, false); //skip all events which time diff is not of interest to us
	//mymeas.SkipEventsTimeDiffCut(12, 13, -100, -1.9, false); //skip all events which time diff is not of interest to us
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
	//mymeas.Print_GetTimingCFD_diff({10,11}, {12,13}, -7.5, 12.5, 1, 200, 0, 5.5, "RS", true);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 20.;	// lower integration window in ns rel. to max
	float intwindowplus = 30.;	// upper integration window in ns rel. to max
	float findmaxfrom = 90.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	float plotrangestart = 0; 	// decide on the
	float plotrangeend = 8000;	// plotrange of x-axis
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
	//vector<int> phi_chx = {225-4, 270-4, 315-4, 0-4, 45-4, 90-4, 135-4, 180-4}; //ordered from channel 0 to channel 7; my channel alignment was a bit different from Alex's/Andrea's
	vector<int> phi_chx_andrea = {0, 315, 270, 225, 45, 90, 135, 180}; //ordered from channel 0 to channel 7
	vector<int> phi_chx_andrea_even = {0, 270, 45, 135}; vector<int> phi_chx_andrea_odd = {315, 225, 90, 180};
	//vector<int> phi_chx_even = {225-4, 315-4, 45-4, 135-4}; vector<int> phi_chx_odd = {270-4, 0-4, 90-4, 180-4};
	//vector<float> ly_C0 = {1361.92, 768.23, 1354.06, 765.312, 1385.92, 873.016, 1433.86, 774.65}; //selection with small PS + larger int-window2 (larger_intw_catch_prim_gamma2.png) --> corrPSL2
	//vector<float> ly_C0 = {1558.22, 908.343, 1535.05, 909.704, 1578.06, 1045.79, 1656.23, 907.486}; //selection with small PS + larger int-window3 (larger_intw_catch_prim_gamma3.png) --> corrPSL3
	//vector<float> ly_C0 = {1526.43, 893.474, 1505.19, 894.504, 1548.95, 1019.63, 1614.68, 891.464}; //selection with small PS + larger int-window4 (larger_intw_catch_prim_gamma4.png) --> corrPSL4
	//vector<float> ly_C0 = {2093.183, 1065.871, 2054.536, 1004.251, 2104.203, 1195.281, 2176.630, 1039.790}; //trial
	//vector<float> ly_C0 = {1450.57, 871.348, 1436.27, 872.943, 1468.98, 976.537, 1532.12, 865.508}; //selection with small PS + larger int-window2 (larger_intw_catch_prim_gamma2.png) + improved baseline; best and final for my data --> corrPSL5; best for my data
	vector<float> ly_C0_andrea = {2193.156, 2291.847, 2022.790, 1966.711, 2053.464, 2045.395, 2177.515, 2021.000}; //andrea_C0: selection with 1ns window (time_cut_andrea_C0_new.png) + larger intw (intw_andrea_big.png); best for andreas data
	//vector<float> ly_C0_andrea = {2123.517, 2225.649, 1952.350, 1898.769, 1962.450, 1959.998, 2097.645, 1926.371}; //andrea_C0: selection with 1.5ns window (time_cut_andrea_C0_1.5_new.png) + larger intw (intw_andrea_big.png)
	vector<float> ly_C0_andrea_even = {2193.156, 2022.790, 2053.464, 2177.515}; vector<float> ly_C0_andrea_odd = {2291.847, 1966.711, 2053.464, 2021.000};
	//vector<float> ly_C0_even = {1450.57, 1436.27, 1468.98, 1532.12}; vector<float> ly_C0_odd = {871.348, 872.943, 976.537, 865.508};
	//vector<int> SiPMchannels = {0, 1, 2, 3, 4, 5, 6, 7};
	vector<int> SiPMchannels = {0, 2, 4, 6}; //for the odd/even analysis
	//mymeas.Print_Phi_ew(phi_chx, ly_C0, SiPMchannels, intwindowminus, intwindowplus, findmaxfrom, findmaxto, 180, false, false);
	mymeas.Print_Phi_ew(phi_chx_andrea_even, ly_C0_andrea_even, SiPMchannels, intwindowminus, intwindowplus, findmaxfrom, findmaxto, 180, true, true);

	// plot cfd-results (between t=110 ns and t=140 ns) for all channels and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1, 200);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// for getting average lightyield, do which_fit=0 and look at the histogramms in the .root-file and manually extract the means
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotrangestart, plotrangeend, 200, fitstart, fitend, channels_to_fit, which_fit);
	//for(int i = 0; i < mymeas.mean_integral.size(); i++) cout << fixed << setprecision(3) << mymeas.mean_integral[i] << endl;

	// suppress graphic output
	gROOT->SetBatch(kTRUE); // TRUE enables batch-mode --> disables graphic output (all prints before this will still be shown)

	//plot range for individual waveforms
	double ymin = -5;
	double ymax = 50;

	// plot waveforms for certain events with integration window
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event_number, ymin, ymax);

	// plot individual waveforms for some events for debugging
	vector<int> plotable_wfs = {};
	for (int i = 0; i < mymeas.nevents; i += 1) if (!mymeas.skip_event[i]) plotable_wfs.push_back(i+1); //apparently, there is some +1 bias in event numbering
	for (int i = 0; i < static_cast<int>(plotable_wfs.size()); i += static_cast<int>(plotable_wfs.size()) / 20) {
   		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotable_wfs[i], ymin, ymax);
	}
	gROOT->SetBatch(kFALSE);

}