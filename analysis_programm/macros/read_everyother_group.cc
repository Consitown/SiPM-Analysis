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

	// WARNING: the name of the folder where the .bin-files are stored must be the same as the name of the bin files
	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "19cosmics_pcbj_everyothergroup_vb42_PSleft/"; // only every other group of the SiPM's
        break;
	}
	case(1): {
		path += "22cosmics_pcbj_everyothergroup_vb42_PSmid/"; // has optical coupling
		break;
	}
	case(2): {
		path += "23cosmics_pcbj_everyothergroup_vb42_PSright/"; // PS positions are left (closer to door), mid, right (farther to door)
		break;
	}
	case(3): {
		path += "24cosmics_pcbj_everyothergroup_vb42_PSleft_new/"; // PS positions are left (closer to door), mid, right (farther to door)
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
	mymeas.ReadFile(path, false, 0, path + "everyother_results.root");

	// only plot channels specified below. Leaving it empty will plot all channels
	mymeas.plot_active_channels={0,1,2,3,4,5,6,7};

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	//mymeas.CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	//mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 0, 0, false, false, 8);

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Sytax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start) --> last argument is selecting inverse/normal cfd: true results in normal cfd
	mymeas.GetTimingCFD(cfd_x, 110, 150, 0, false); // this creates the timing_results matrix

	// create histogram of delta t: average cfd-time upper PMT's - average cfd-time lower PMT's
	// PMT's are at channels 10-13; 10 - right upper PMT, 11 - left upper PMT, 12 - right lower PMT, 13 - left lower PMT (perspective from door)
	// but waveforms of PMT come directly after wf's of SiPM's --> wf 0 - channel 0; wf 1 - ch 1; ... ; wf 8 - ch 10; wf 9 - ch 11; ...
	gStyle->SetOptStat("nemr");
	TString his_name(Form("delta_t at %0.2f", cfd_x));
	TString his_title(Form("Average time diff lower to upper PMT's at cfd of %0.2f", cfd_x));
	TH1* his = new TH1F(his_name, his_title, 100, -50, 50);
	TCanvas* hisc = new TCanvas(his_name, his_title, 1600, 1000);
	int counter;
	for(int i=0 ; i < mymeas.nwf; i += mymeas.nchannels){
		float delta_t = ((mymeas.timing_results[i+10][1] + mymeas.timing_results[i+11][1])/2) - ((mymeas.timing_results[i+9][1] + mymeas.timing_results[i+8][1])/2); // average cfd-delta_t
		his->Fill(delta_t);
	}
	his->GetXaxis()->SetTitle("time [ns]");
	his->GetYaxis()->SetTitle("#Entries");
	his->Draw();
	hisc->Update();
	mymeas.root_out->WriteObject(his, "delta_t");

	// apply cut for time difference between two channels
	//mymeas.SkipEventsTimeDiffCut(10, 13, 1, 6, false);

	// print events above a threshold to identify interesting events
	// mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	//**********//
	// Plotting + Fitting
	//**********//

	// plot sums of all events per channel --> get parameters
	// Syntax: ...(float ystart, float yend, bool doaverage, bool normalize, double shift, double sigma, bool doconv)
	//mymeas.PlotChannelSums(-1e4, 3.6e5, false);

	// investigate charge spectrum. For the integration values, look at the plots from PlotChannelSums. --> you can determine findmaxfrom and findmaxto
	float intwindowminus = 30.;	// lower integration window in ns rel. to max
	float intwindowplus = 100.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 130.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)
	float plotrangestart = 0; 	// decide on the
	float plotrangeend = 12000;	// plotrange of x-axis
	float fitstart = 25;		// start of the fit (x-axis)
	float fitend = 350;			// end of the fit
	int channels_to_fit = 0; 	// numbers of channels to apply the fit to (counts like this: i=0;i<channels_to_fit;i++); thats why one should use the first channels for data and the later channels for triggering
	int which_fit = 1;			// decide on which fit-function you want; options: default (no value besides 1-6) - default SiPM fit function
								// 1 - landau gauss convolution for large number of photons
								// 2 - biased: if pedestal is biased because of peak finder algorithm
								// 3 - SiPM fit function with exponential delayed afterpulsing
								// 4 - ideal PMT fit function (5 is similar)
								// 6 - PMT fit function with biased pedestal
	
	// old timing histogramm (should not use for cfd, so cfd_x should be 1)
	//mymeas.PrintTimeDist(110, 140, 105, 145, 200, 1, cfd_x);

	// plot results between t=110 ns and t=140 ns and fit gauss (thats what the 1 is for)
	//mymeas.Print_GetTimingCFD(110, 140, 1);

	//plot constant fration descrimination (old)
	// mymeas.PlotConstantFrationDescrimination(100, 150);

	// plot all pseudo charge spectrum of channels (real charge spectrum would be gained by multiplying every integrated value [x-axis] with 1/resistance_of_sipms)
	// Syntax: ...(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf)
	//mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, plotrangestart, plotrangeend, 200, fitstart, fitend, channels_to_fit, which_fit);

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