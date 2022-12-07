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
	int run = 0;

	// WARNING: the name of the folder where the .bin-files are stored must be the same as the name of the bin files
	switch (which) { //specify folders to run below, ALL bin files in this folder will be used.
	case(0): {
		path += "5time_testrun_PS_better/"; // PS run with the right trigger settings for CFD
		run = 5;
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

	// only plot channels specified below. Leaving it empty will plot all channels
	//mymeas.plot_active_channels={};

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	//mymeas.CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel)
	mymeas.CorrectBaselineMinSlopeRMS(100, true, 10, 10, 0, false, false, 9); //smoothing improves the result significantly

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Sytax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start) --> last argument is selecting inverse/normal cfd: true results in normal cfd
	bool modus = true;
	mymeas.GetTimingCFD(cfd_x, 110, 150, 0, modus); // this creates the timing_results matrix
	//if(modus){TString pdf_name(Form("normal"));}else{TString pdf_name(Form("inverse"));}

	// create histogram of delta t: average cfd-time upper PMT's - average cfd-time lower PMT's
	// PMT's are at channels 10-13; 10 - right upper PMT, 11 - left upper PMT, 12 - right lower PMT, 13 - left lower PMT (perspective from door)
	// but waveforms of PMT come directly after wf's of SiPM's --> wf 0 - channel 0; wf 1 - ch 1; ... ; wf 8 - ch 10; wf 9 - ch 11; ...
	// match channel number to channel index
	int ch_to_plot[4] = {10,11,12,13};
	int fir_ch = 0; int sec_ch = 0; int thi_ch = 0; int for_ch = 0;
	for (int i = 0; i < mymeas.active_channels.size(); i++) {
		if (mymeas.active_channels[i] == ch_to_plot[0]) fir_ch = i;
		if (mymeas.active_channels[i] == ch_to_plot[1]) sec_ch = i;
		if (mymeas.active_channels[i] == ch_to_plot[2]) thi_ch = i;
		if (mymeas.active_channels[i] == ch_to_plot[3]) for_ch = i;
	}
	gStyle->SetOptStat("nemr"); //draws a box with some histogram parameters
	TString his_name(Form("delta_t at %0.2f", cfd_x)); //the name of the histogram
	TString his_title(Form("Average time diff lower to upper PMT's at cfd of %0.2f", cfd_x)); //the title of the histogram
	TH1* his = new TH1F(his_name, his_title, 200, -10, 15); //new Histogram
	TCanvas* hisc = new TCanvas(his_name, his_title, 1600, 1000); //new canvas to save the histogram on
	for(int i=0 ; i < mymeas.nwf; i += mymeas.nchannels){ //loop through all the events
		float delta_t = ((mymeas.timing_results[i+thi_ch][1] + mymeas.timing_results[i+for_ch][1])/2) - ((mymeas.timing_results[i+fir_ch][1] + mymeas.timing_results[i+sec_ch][1])/2); // average cfd-delta_t
		his->Fill(delta_t); //fill the data into the histogram
	}
	his->Fit("gaus", "L","same"); //Gaussian fit
	his->GetXaxis()->SetTitle("time [ns]"); //naming the axes
	his->GetYaxis()->SetTitle("#Entries");
	his->Draw(); //don't know if this is really necessary or only for the graphic output
	hisc->Update(); //put the histogram on the canvas
	mymeas.root_out->WriteObject(his, "delta_t"); //save the histogram in a .root-file
	TString pdf_filename(Form("delta_t_%2.0f_normal_run%1d.pdf", 100*cfd_x, run));
	gErrorIgnoreLevel = kError; //suppress root terminal output
	hisc->Print(pdf_filename); //write the histogram to a .pdf-file (this makes saving in a root-file kinda redundant)
	gErrorIgnoreLevel = kUnset;
	//todo for inclusion in ReadRun: include the skip_events array, make the 4 PMT-channel a parameter array
	//make the plotrange + bin number a parameter, include checking whether GetTimingCFD was done prior

	// apply cut for time difference between two channels
	//mymeas.SkipEventsTimeDiffCut(10, 13, 1, 6, false);

	// print events above a threshold to identify interesting events
	// mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

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