
#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

using namespace std;

void timing_example(int which = 0) // main
{
	string path;
	// edit for your fs
	path = "/mnt/c/SHiP/macros/constantin/";

	switch (which) { //specify folder to run below
	default: {
		path += "exampledata_timing/"; 
		// cosmics box example data (only 500 events)
		break;
	}
	}

	// initialize class
	ReadRun mymeas(0);
	
	// read data
	mymeas.ReadFile(path, true, 0, path + "/results.root");

	// apply baseline correction to ALL waveforms
	mymeas.CorrectBaselineMin(30, 0, 5., 400, 250, 0);

	// get timing for 30% CFD between t=100 ns and t=140 ns, 0 means no smoothing, false to search from max (true to search from start)
	mymeas.GetTimingCFD(0.3, 100, 140, 0, false);
	cout << "\n\n test : " << mymeas.timing_results[1][1] << "\n\n";

	// apply cut for time difference between two channels (ch14 and ch26, events with time differences <1 ns or >5ns will be cut)
	mymeas.SkipEventsTimeDiffCut(14, 26, 1, 5, false);

	// plot results between t=100 ns and t=140 ns and fit gauss
	mymeas.Print_GetTimingCFD(100, 140, 1);

	// plotting
	
	// old timing plot
	mymeas.PrintTimeDist(100, 140, 90, 150, 200, 1, .3);

	// plot sums of all events per channel
	mymeas.PlotChannelSums(true, false);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 15.;	// lower integration window in ns rel. to max
	float intwindowplus = 85.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume signal from muon arrives between here ...
	float findmaxto = 130.;		// ... and here (depends on trigger delay setting)

	// plot all channels
	// integral
	//mymeas.PrintChargeSpectrum_pars = { 1e4, 0.5, .15, 6., 6., 40., 5 };
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -100, 1900, 500, -100, 1900);
	// amplitude
	//mymeas.PrintChargeSpectrum_pars = { 1e4, 0.5, .15, .5, .5, 2., 2 };
	mymeas.PrintChargeSpectrum(0, 0, findmaxfrom, findmaxto, 0, 50, 500, 0, 50, -1);
	
	// plot waveforms of individual events
	// plot range
	double ymin = -10;
	double ymax = 120;
	// plot waveforms for certain events with integration window and timing info
	gROOT->SetBatch(kTRUE); // only write to root file
	for (int i = 1; i < mymeas.nevents; i += static_cast<int>(mymeas.nevents / 50)) mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, i, ymin, ymax);
	gROOT->SetBatch(kFALSE);
}
