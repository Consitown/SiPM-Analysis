#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>

// you'll need to compile the analysis (ReadRun.cc) first! --> use command 'make -f makefile' in the directory where ReadRun.cc and ReadRun.h are

using namespace std;

void read_sipm_cal_pcb(int which) // main
{
	string path;

	// edit for your work-directory
	path = "/mnt/d/Work_SHK_Bachelor/analysis_programm/measurements/";

	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "6sipm_cal_pcbj_vb43_tune8270_part/"; // vb = Voltage bias (applied voltage); note: the power supply actually outputs 891 mV less
		break;
	}
	case(1): {
		path += "7sipm_cal_pcbj_vb43_tune8230_part/"; // tune = laser tune
		break;
	}
	case(2): {
		path += "8sipm_cal_pcbj_vb41_tune8200_part/"; // part stands for partial events option in WaveCatcher enabled
		break;
	}
	case(3): {
		path += "9sipm_cal_pcbb_vb55_tune8200_part/"; // placeholder
		break;
	}
	case(4): {
		path += "13sipm_cal_pcba_vb551_tune8200_part/"; // placeholder
		break;
	}
	case(5): {
		path += "15sipm_cal_pcbc_vb56_tune8120_part/"; // placeholder
		break;
	}
	case(6): {
		path += "16sipm_cal_pcbc_vb56_tune8080_part/"; // placeholder
		break;
	}
	case(7): {
		path += "17sipm_cal_pcbd_vb57_tune8060_part/"; // placeholder
		break;
	}
	case(8): {
		path += "18sipm_cal_pcbd_vb571_tune8020_part/"; // placeholder
		break;
	}
	case(9): {
		path += "20sipm_cal_pcbb_vb565_tune8380_part/"; // placeholder
		break;
	}
	case(10): {
		path += "21sipm_cal_pcbb_vb566_tune8360_part/"; // placeholder
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
	// read data; mymeas.ReadFile(path, true, 0, path + "/cal_results.root") for an explizit output file
	mymeas.ReadFile(path, false, 0, path + "/cal_results.root");

	// only plot channels specified below. Leaving it empty will plot all channels
	mymeas.plot_active_channels={0,1,2,3,4,5,6,7};

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	//mymeas.CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 0, 0, false, false, 8);

	// cfd-stuff; here: get the cfd-times off all waveforms
	//float cfd_x = .3;
	//mymeas.GetTimingCFD(cfd_x, 110, 140, 0);

	// apply cut for time difference between two channels
	//mymeas.SkipEventsTimeDiffCut(10, 13, 1, 6, false);

	// print events above a threshold to identify interesting events
	// mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	//**********//
	// Plotting + Fitting
	//**********//

	// plot sums of all events per channel --> get parameters
	// Syntax: ...(bool doaverage, bool normalize, double shift, double sigma, bool doconv)
	//mymeas.PlotChannelSums(false);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 1.;	// lower integration window in ns rel. to max
	float intwindowplus = 1.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 130.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	float plotrangestart = 0; // decide on the
	float plotrangeend = 35;	// plotrange of x-axis
	float fitstart = 0;			// start of the fit (x-axis)
	float fitend = 35;			// end of the fit
	int channels_to_fit = 1; 	// numbers of channels to apply the fit to (counts like this: i=0;i<channels_to_fit;i++); thats why one should use the first channels for data and the later channels for triggering
	int which_fit = 3;			// decide on which fit-function you want; options: default (no value besides 1-6) - default SiPM fit function
								// 1 - landau gauss convolution for large number of photons
								// 2 - biased: if pedestal is biased because of peak finder algorithm
								// 3 - SiPM fit function with exponential delayed afterpulsing
								// 4 - ideal PMT fit function (5 is similar)
								// 6 - PMT fit function with biased pedestal
	
	// old timing histogramm (should not use for cfd)
	//mymeas.PrintTimeDist(110, 140, 105, 145, 200, 1, cfd_x);

	// plot results between t=110 ns and t=140 ns and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1);

	//plot constant fration descrimination
	// mymeas.PlotConstantFrationDescrimination(100, 150);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotrangestart, plotrangeend, 200, fitstart, fitend, channels_to_fit, which_fit);

	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins)
	// PrintChargeSpectrumPMT will apply a fit automatically from rangestart to rangeend (which are also the boundaries for the plot)
	//mymeas.PrintChargeSpectrumPMT(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 5, 300, 202);

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