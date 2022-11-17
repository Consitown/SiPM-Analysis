#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>

// you'll need to compile the analysis (ReadRun.cc) first! --> use command 'make -f makefile' in the directory where ReadRun.cc and ReadRun.h are

using namespace std;

void read_cfd_PS(int which = 0) // main
{
	string path;

	// edit for your work-directory
	path = "/mnt/d/Work_SHK_Bachelor/analysis_programm/measurements/";

	// WARNING: the name of the folder where the .bin-files are stored must be the same as the name of the bin files
	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "5time_testrun_PS_better/"; // PS run with the right trigger settings for CFD
		break;
	}
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);

	// read data; mymeas.ReadFile(path, true, 0, path + "/cfd_results.root") for an explizit output file
	mymeas.ReadFile(path, true, 0, path + "/cfd_results.root");

	// only plot channels specified below. Leaving it empty will plot all channels
	//mymeas.plot_active_channels={};

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	mymeas.GetTimingCFD(cfd_x, 110, 140, 0);

	//write some function for spectrum of time difference of channel 10 to 11,12,13 --> justify cuts (cuts should be able to select certain event locations)

	// apply cut for time difference between two channels
	//mymeas.SkipEventsTimeDiffCut(10, 13, 1, 6, false);

	// print events above a threshold to identify interesting events
	// mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	//**********//
	// Plotting
	//**********//

	// plot sums of all events per channel --> get parameters
	//mymeas.PlotChannelSums(false);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 3.;	// lower integration window in ns rel. to max
	float intwindowplus = 5.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	
	// old timing histogramm (should not use for cfd)
	//mymeas.PrintTimeDist(110, 140, 105, 145, 200, 1, cfd_x);

	// plot results between t=110 ns and t=140 ns and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1);

	//plot constant fration descrimination
	// mymeas.PlotConstantFrationDescrimination(100, 150);


	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -10, 350, 300, 0, 0, 0, 0);

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