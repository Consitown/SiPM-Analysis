#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TSystem.h>
#include <TCanvas.h>

// you'll need to compile the analysis first!

using namespace std;

void read_exampledata() // main
{
	int which = 0; //select meas

	string path;

	// edit for your fs
	path = "/mnt/d/Work_SHK_Bachelor/analysis-program/measurements/";

	switch (which) { //specify folder to run below, ALL bin files in this folder will be used
	case(0): {
		path += "exampledata/"; //
		break;
	}//
	default: {
		cout << "\nerror: path to data not specified" << endl; // default
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);

	// read data
	mymeas.ReadFile(path, true);

	// only plot channels specified below (can be uncommented to analyze all active channels in data)
	int channel_to_plot = 9;
	mymeas.plot_active_channels.push_back(channel_to_plot);
	mymeas.plot_active_channels.push_back(14);
	mymeas.plot_active_channels.push_back(15);

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	// print events above a threshold to identify interesting events
	mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	////plotting

	// plot sums of all events per channel
	mymeas.PlotChannelSums(true);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 3.;	// lower integration window in ns rel. to max
	float intwindowplus = 5.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)

	// plot all charge spectrum of channels
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -1e2, 2.e3, 500, 0, 0, 0);

	mymeas.PrintChargeSpectrumPMT(0, 0, findmaxfrom, findmaxto, -2e1, 1.8e2, 202);

	// timing of maximum
	mymeas.PrintTimeDist(findmaxfrom, findmaxto, findmaxfrom - 5, findmaxto + 5, 60);

	// plot waveforms of individual events
	int event1 = 68;
	int event2 = 79;

	//plot range
	double ymin = -5;
	double ymax = 25;

	// plot waveforms for certain events with integration window
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event1, ymin, ymax);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event2, ymin, ymax);

	// only write to root file
	gROOT->SetBatch(kTRUE);
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 10)) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	}
}