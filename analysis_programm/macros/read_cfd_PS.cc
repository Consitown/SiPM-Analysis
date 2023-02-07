#include <iostream>
#include <fstream> //for opening and edidting .txt-files
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
	case(5): {
		path += "29_normal_PS_coinc_with_ortho_PS_C1_20k/";
		run = 29;
		break;
	}
	case(6): {
		path += "32_nor_PS_coinc_w_ort_PS_at_7toC1/";
		run = 32;
		break;
	}
	case(7): {
		path += "33_nor_PS_coinc_w_ort_PS_at_7toC2/";
		run = 33;
		break;
	}
	case(8): {
		path += "34_nor_PS_coinc_w_ort_PS_at_205toC2/";
		run = 34;
		break;
	}
	case(9): {
		path += "35_nor_PS_coinc_w_ort_PS_at_1925toleft/";
		run = 35;
		break;
	}
	case(10): {
		path += "36_normal_PS_coinc_with_ortho_PS_C0/";
		run = 36;
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
	int counter_ortho = 0;
	for (int i = 0; i < mymeas.skip_event.size(); i++) {if (mymeas.skip_event[i]) ++counter_ortho;} //a bit debugging
	cout << "Number of orthogonal events: " << counter_ortho << endl;
	mymeas.skip_event.flip();

	//skip weird events: -5 mV for identifying out burst/weird events (high frequency oscillations) in large PS (-7 only if signals inverted)
	vector<double> thresholds_weird = {-5, -5, -5, -5, -5, -5};
	mymeas.SkipEventsPerChannel(thresholds_weird, 100, 200, false);

	// cfd-stuff; here: get the cfd-times off all waveforms
	float cfd_x = .3;
	// Syntax: ...(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, bool doconv, bool use_spline) --> fifth argument is selecting inverse/normal cfd: true results in normal cfd
	mymeas.GetTimingCFD(cfd_x, 110, 150, 3, true, false, false); // this creates the timing_results matrix

	//plotting amplitude spectrum; mymeas.timimng_results[waveform][1] contains cfd-time
	int channel1 = 12; int channel2 = 13;
	// match channel number to channel index
	int ch_index1, ch_index2 = 0;
	for (int i = 0; i < mymeas.active_channels.size(); i++) {
		if (mymeas.active_channels[i] == channel1) ch_index1 = i;
		if (mymeas.active_channels[i] == channel2) ch_index2 = i;
	}

	gStyle->SetOptStat("nemr"); //draws a box with some histogram parameters

	TString his_name(Form("ch%2d-ch%2d_at_cfd_%0.2f", channel2, channel1, cfd_x));
	TCanvas* combi = new TCanvas("combi", "different filters", 600, 400);
	int nbins = 200; int min = -15; int max = 15;
	TH1* his1 = new TH1F("all_events", "Combined_his_"+his_name, nbins, min, max);
	TH1* his2 = new TH1F("ortho_skipped", "Combined_his_"+his_name, nbins, min, max);
	TH1* his3 = new TH1F("ortho_only", "Combined_his_"+his_name, nbins, min, max);

	for (int i=0 ; i < mymeas.nevents ; i++){ //loop through all the events
		his1->Fill(mymeas.timing_results[i*mymeas.nchannels+ch_index2][1]-mymeas.timing_results[i*mymeas.nchannels+ch_index1][1]); //all events
		if (mymeas.skip_event[i]) {
			his2->Fill(mymeas.timing_results[i*mymeas.nchannels+ch_index2][1]-mymeas.timing_results[i*mymeas.nchannels+ch_index1][1]); //ortho events skipped
		}
		else {
			his3->Fill(mymeas.timing_results[i*mymeas.nchannels+ch_index2][1]-mymeas.timing_results[i*mymeas.nchannels+ch_index1][1]); //ortho events
		}
	}

	his1->GetXaxis()->SetTitle("time [ns]");
	his1->GetYaxis()->SetTitle("#Entries");
	his1->Draw(); gPad->Update(); // draw the first histogram on canvas // force the update, so stat box gets drawn
	TPaveStats* t = (TPaveStats*)his1->FindObject("stats"); // find the stats, only works if the stat-box is drawn in first place
	t->SetY1NDC(.77); t->SetY2NDC(.92); // new y start position // new y end position
	his2->SetLineColor(kRed);
	his2->Draw("sames"); gPad->Update(); // draw the second histogram on canvas // force the update, so stat box gets drawn
	TPaveStats* s = (TPaveStats*)his2->FindObject("stats"); // find the stats, only works if the stat-box is drawn in first place
	s->SetY1NDC(.6); s->SetY2NDC(.75); // new y start position // new y end position
	his3->SetLineColor(kOrange);
	his3->Draw("sames"); gPad->Update(); // draw the third histogram on canvas // force the update, so stat box gets drawn
	TPaveStats* u = (TPaveStats*)his3->FindObject("stats"); // find the stats, only works if the stat-box is drawn in first place
	u->SetY1NDC(.43); u->SetY2NDC(.58); // new y start position // new y end position
	combi->Update(); //put everything on the canvas
	
	TString pdf_filename(Form("run%2d_%2d-%2d_combined.pdf", run, channel2, channel1));
	gErrorIgnoreLevel = kError; //suppress root terminal output
	combi->Print(pdf_filename); //write the histogram to a .pdf-file (this makes saving in a root-file kinda redundant)
	gErrorIgnoreLevel = kUnset;	

	// Positioncuts (trial)
	//mymeas.SkipEventsTimeDiffCut(10, 11, -2.5, 2.5, false);
	//mymeas.SkipEventsTimeDiffCut(12, 13, -2.5, 2.5, false);


	// prints some stats for events above a threshold into the terminal to identify interesting events
	// Syntax: ...(float threshold, bool max, bool greater, double from, double to, bool verbose)
	//mymeas.FractionEventsAboveThreshold(5, true, true, 200, 250, false);

	//**********//
	// Plotting
	//**********//

	// plot sums of all events per channel --> get parameters
	// Syntax: ...(bool doaverage, bool normalize, double shift, double sigma, bool doconv)
	//mymeas.PlotChannelSums(false);

	// Syntax: ...(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int do_fit, int nbins, float fitrangestart, float fitrangeend, string fitoption, bool set_errors)
	//from entering lab: upper right: 10, upper left: 11, lower right: 12, lower left: 13; channel 8 is upper orthogonal PMT, 9 is lower
	//mymeas.Print_GetTimingCFD_diff({12}, {13}, -15, 15, 0, 200, -8, 8, "RS", false);

	// Syntax: ... (float rangestart, float rangeend, int do_fit, int nbins, string fitoption)
	//mymeas.Print_GetTimingCFD(110,140,1,200,"S"); //for channelwise cfd

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