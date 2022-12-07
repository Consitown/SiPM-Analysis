/// \mainpage Introduction
/// Page is under construction!! \n \n 
/// 
/// You can find the relevant documentation at <a href="classReadRun.html"> ReadRun Class Reference</a> \n \n
///  
/// The code and instructions for getting started can be found at <https://github.com/cscharf-hub/wavecatcher-analysis>
/// 
/// Creators: \n
/// Christian Scharf \n 
/// Contributors: \n 
/// Doramas Jimeno Sanchez \n
/// Alessia Brignoli \n
/// Ben Skodda \n 
/// Christophe Mullesch \n
/// Constantin Eckardt \n
/// Alexander Vagts \n



#include "ReadRun.h"


ClassImp(ReadRun)

/// @brief Constructor of the class with arguments to filter noise events in the cosmics setup. Default values do nothing 
/// 
/// PMT_threshold -> set to 0 to do nothing. If a value is set, events where the maximum value is > PMT_threshold in \n
/// channel > 7 are removed from the analysis (for cosmics setup). Used to filter out events where the PMTs have triggered on picked-up radio frequency noise signals
ReadRun::ReadRun(double PMT_threshold, int channels_above_threshold) {

	cout << "\ninitializing ..." << endl;

	skip_event_threshold = PMT_threshold;
	skip_event_threshold_nch = channels_above_threshold;
	if (skip_event_threshold > 0) {
		cout << "\n removing events where channels 9-16 have entries exceeding +" << skip_event_threshold << " mV amplitude in at least " << channels_above_threshold << " PMTs\n\n";
	}
	else if (skip_event_threshold < 0) {
		cout << "\n removing events where channels 9-16 have entries below " << skip_event_threshold << " mV amplitude in at least " << channels_above_threshold << " PMTs\n\n";
	}
	nwf = 0;

	PrintChargeSpectrum_cnt = 0;
	PrintChargeSpectrumPMT_cnt = 0;
	PrintChargeSpectrumPMTthreshold_cnt = 0;
}

/// @brief Routine to read files created by the wavecatcher.
/// 
/// Reads the data and can already do simple data manipulation during reading: \n
///   - Change polarity of range of channels with parameters (explained below). \n 
///   - Shift all waveforms to a constant fraction such they all start at the same time (-> see ReadRun::Shift_WFs_in_file_loop) \n 
///   - Simple baseline correction by calling CorrectBaseline() before calling ReadFile(). \n 
/// Stores the data as ROOT histograms in TClonesArray rundata. \n \n 
/// 
/// Reader modified from WaveCatcher binary -> root converter by manu chauveau@cenbg.in2p3.fr \n 
/// See https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2
/// 
/// @param path Path to the data. The name of the folder needs to match the file names (data\data.bin_001 etc.).
/// @param change_polarity Set true to change polarity (sign) of certain channels (see "change_sign_from_to_ch_num" below).
/// @param change_sign_from_to_ch_num All channels \f$\geq\f$"change_sign_from_to_ch_num" will be inverted if "change_polarity" is true. \n 
/// If negative number all channels \f$\leq\f$abs("change_sign_from_to_ch_num") will be inverted if "change_polarity" is true.
/// @param out_file_name Name of the .root file which stores the results.
/// @param debug Set true to increase the verbosity.
/// @todo Remove "change_polarity" and make back-compatible. 
void ReadRun::ReadFile(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debug) {

	// save output path
	data_path = path;

	printf("+++ saving analysis results in '%s' ...\n\n", out_file_name.c_str());
	root_out = TFile::Open(out_file_name.c_str(), "recreate");

	// Wavecatcher hardware/software properties
	SP = 0.3125;					// ns per bin
	coef = 2.5 / (4096 * 10);		//?????
	binNumber = 1024;				//default: 1024, hard coded later on so it can be removed
	const int nChannelsWC = 64;			//max number of channels default: 32

	rundata = new TClonesArray("TH1F", 1e7); //raw data will be stored here as TH1F
	rundata->BypassStreamer(kFALSE);  //Doramas: I don't know why is it used, but it's better to use when working with TClonesArray
	TClonesArray& testrundata = *rundata;

	// verbosity
	bool debug_header = debug;
	bool debug_data = debug;

	unsigned short output_channel;
	unsigned int output_event;
	unsigned long long int output_tdc;
	unsigned short output_nbchannels;
	unsigned short read_channels;

	amplValuessum = new double* [nChannelsWC]; //sum of all wf for each channel
	for (int i = 0; i < nChannelsWC; i++) {//init
		amplValuessum[i] = new double[binNumber];
		for (int k = 0; k < binNumber; k++) amplValuessum[i][k] = 0.;
	}

	maxSumBin = new int[nChannelsWC];

	//Start reading the raw data from .bin files.
	stringstream inFileList = list_files(path.c_str(), ".bin"); //all *.bin* files in folder path
	int nitem = 1;
	string fileName;
	int file_counter = 0;
	int currentPrint = -1;
	int wfcounter = 0;
	int event_counter = 0;
	int removed_event_counter = 0;

	while (inFileList >> fileName) {
		// file loop

		fileName = path + fileName;
		ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);

		bool has_measurement = false;

		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}
		printf("+++ reading '%s' ...\n", fileName.c_str());

		// Header
		string header_line;
		// HEADER 1 //
		//
		// "=== DATA FILE SAVED WITH SOFTWARE VERSION: V?.??.? ==="
		//
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t header_version_first = header_line.find_last_of('V');
		size_t header_version_last = header_line.find_first_of(' ', header_version_first);
		string software_version = header_line.substr(header_version_first, header_version_last - header_version_first);
		if (debug_header) printf("    |- data version = '%s'\n", software_version.data());

		//if (software_version == "V2.9.13")
		//	;
		//else if (software_version == "V2.9.15")
		//	;
		//else if (debug_header) printf("*** unsupported data version\n");

		// HEADER 2 //
		// "=== WAVECATCHER SYSTEM OF TYPE ?? WITH ?? CHANNELS AND GAIN: ??? ==="
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 3 //
		// === Rate coincidence masks ... === Posttrig in ns for SamBlock ... ===
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 4 //
		// V2.9.13: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1
		// V2.9.15: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1 == MEASUREMENTS: 0 ===
		getline(input_file, header_line, '\n');

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t nsamples_first = 1 + header_line.find_last_of('[');
		size_t nsamples_last = header_line.find_first_of(']', nsamples_first);
		string nsamples_str = header_line.substr(nsamples_first, nsamples_last - nsamples_first);

		int nsamples = atoi(nsamples_str.data());
		if (debug_header) printf("    |- data sample  = %d\n", nsamples);

		size_t nchannels_first = 10 + header_line.find("ACQUIRED: ", nsamples_first);
		size_t nchannels_last = header_line.find_first_of(' ', nchannels_first);
		string nchannels_str = header_line.substr(nchannels_first, nchannels_last - nchannels_first);

		nchannels = atoi(nchannels_str.data());
		if (debug_header) printf("    |- nchannels    = %d\n", nchannels);

		if (software_version == "V2.9.15" || software_version == "V2.9.16" || software_version == "V2.10.1") {
			size_t has_measurement_first = 14 + header_line.find("MEASUREMENTS: ", nsamples_first);
			size_t has_measurement_last = header_line.find_first_of(' ', has_measurement_first);
			string has_measurement_str = header_line.substr(has_measurement_first, has_measurement_last - has_measurement_first);
			has_measurement = atoi(has_measurement_str.data());
		}
		else {
			//if (software_version == "V2.9.13") {
				// V2.9.13 has always measurement stored
				// (everything is set to 0 when disabled!)
			has_measurement = true;
		}

		if (debug_header) printf("    `- measurement  = %d\n", has_measurement);

		// end of header reader

		event_data an_event;

		while (input_file.read((char*)(&an_event), sizeof(an_event))) {
			//event loop
			int event_flag_cnt = 0;

			if (debug_data) printf("%03d has %d channels\n", an_event.EventNumber, an_event.nchannelstored);

			output_event = an_event.EventNumber;
			output_tdc = an_event.TDCsamIndex;
			output_nbchannels = an_event.nchannelstored;

			if (debug_data && output_event % 200 == 0) printf("EventNr: %d, nCh: %d\n", output_event, output_nbchannels);

			//cout << "EvN:" << an_event.EventNumber << " EpT:" << an_event.EpochTime << " Yr:" << an_event.Year << " TDCt:" << an_event.TDCsamIndex;

			// do analysis only for limited range of channels to reduce memory usage for large datasets with many channels and many events
			int start_at_ch = 0;
			if (start_read_at_channel < output_nbchannels && start_read_at_channel >= 0) start_at_ch = start_read_at_channel;
			int end_at_ch = output_nbchannels - 1;
			if (end_read_at_channel == -1 && start_read_at_channel != -1) end_read_at_channel = start_read_at_channel;
			else if (end_read_at_channel < output_nbchannels && end_read_at_channel >= 0) end_at_ch = end_read_at_channel;
			read_channels = end_at_ch - start_at_ch + 1;

			if (event_counter == 0) cout << "\nstart at ch " << start_at_ch << " end at ch " << end_at_ch << endl;

			for (int ch = 0; ch < output_nbchannels; ++ch) { // channel loop
				//
				channel_data_with_measurement a_channel_data;
				channel_data_without_measurement a_channel_data_without_measurement;

				if (has_measurement) {
					// read with 'channel_data_with_measurement' struct
					input_file.read((char*)(&a_channel_data), sizeof(channel_data_with_measurement));
				}
				else {
					// read with 'channel_data_without_measurement' struct
					input_file.read((char*)(&a_channel_data_without_measurement), sizeof(channel_data_without_measurement));

					// copy the content into 'channel_data_with_measurement' struct
					a_channel_data.channel = a_channel_data_without_measurement.channel;
					a_channel_data.EventIDsamIndex = a_channel_data_without_measurement.EventIDsamIndex;
					a_channel_data.FirstCellToPlotsamIndex = a_channel_data_without_measurement.FirstCellToPlotsamIndex;
					memcpy(a_channel_data.waveform, a_channel_data_without_measurement.waveform, 1024 * sizeof(short));
				}


				output_channel = a_channel_data.channel;
				if (debug_data) printf("- reading channel %d\n", output_channel);

				//if (event_counter < 10) cout << "\touter: " << ch;
				//---------------------------------------------------------------------------------------------------------------
				if (ch >= start_at_ch && ch <= end_at_ch) {
					//if (event_counter < 10) cout << "\tinner: " << ch;

					if (event_counter == 0) active_channels.push_back(static_cast<int>(output_channel));

					TString name(Form("channel_%02d, event %05d ", output_channel, an_event.EventNumber));
					TString title(Form("Channel %d, event %d raw data", output_channel, an_event.EventNumber));
					auto hCh = (TH1F*)testrundata.ConstructedAt(wfcounter);
					hCh->SetName(name.Data());
					hCh->SetTitle(title.Data());
					hCh->SetBins(binNumber, -0.5 * SP, 1023.5 * SP);

					int nshift = 0;
					if (Shift_WFs_in_file_loop) {
						// shift all waveforms to tWF_CF_bin
						float max = 0.;
						int nmax = 0;
						float cf = tWF_CF;
						int count_fall = 0;

						// do a mini-baseline correction (needed in case a voltage offset is set in the wavecathcer)
						short bsln = (a_channel_data.waveform[0] + a_channel_data.waveform[1] + a_channel_data.waveform[3]) / 3;
						for (int lll = 0; lll < 1024; lll++) a_channel_data.waveform[lll] -= bsln;

						float global_max = TMath::MaxElement(1024, a_channel_data.waveform);
						//nmax = TMath::LocMax(1024, a_channel_data.waveform);
						if (global_max < TMath::Abs(TMath::MinElement(1024, a_channel_data.waveform))) {
							global_max = TMath::Abs(TMath::MinElement(1024, a_channel_data.waveform));
							//nmax = TMath::LocMin(1024, a_channel_data.waveform);
						}


						for (int s = tWF_CF_lo; s < tWF_CF_hi; ++s) {
							if (max < TMath::Abs(a_channel_data.waveform[s])) {
								max = TMath::Abs(a_channel_data.waveform[s]);
								nmax = s;
							}

							// stop search if the current maximum is at least 0.5 * global maximum and if there are at least three consecutive bins where the waveform amplitude is decreasing
							if (max > .5 * global_max && s - nmax < 4 && TMath::Abs(a_channel_data.waveform[s + 2]) + TMath::Abs(a_channel_data.waveform[s + 3]) < TMath::Abs(a_channel_data.waveform[s]) + TMath::Abs(a_channel_data.waveform[s + 1])) count_fall++;
							else count_fall = 0;
							if (count_fall == 3) break;
						}
						cf *= max;
						//if (debug) cout << " || " << cf << ";" << max << ";" << nmax << ";";
						for (int s = nmax; s > tWF_CF_lo; --s) {
							if (cf >= TMath::Abs(a_channel_data.waveform[s])) {
								if (cf - TMath::Abs(a_channel_data.waveform[s]) < TMath::Abs(a_channel_data.waveform[s - 1]) - cf) nshift = tWF_CF_bin - s;
								else nshift = tWF_CF_bin - s - 1;
								break;
							}
						}
						//if (debug) cout << nshift;
					}

					float val = 0.;
					for (int s = 0; s < binNumber; ++s) {
						int shiftind = s - nshift;
						if (shiftind < 0) shiftind += 1023;
						else if (shiftind > 1023) shiftind -= 1023;
						val = a_channel_data.waveform[shiftind] * coef * 1000.;
						if (change_polarity && (output_channel >= change_sign_from_to_ch_num) || (change_sign_from_to_ch_num < 0 && output_channel <= abs(change_sign_from_to_ch_num))) {
							val *= -1.;
						}
						hCh->SetBinContent(s + 1, val);
						//hCh->SetBinError(s, 0.5); //DO NOT USE! Will double the memory usage

						// channel sums
						amplValuessum[ch][s] += static_cast<double>(val);
					}

					// baseline correction
					if (Using_BaselineCorrection_in_file_loop) {
						CorrectBaseline_function(hCh, tCutg, tCutEndg, wfcounter);
					}


					// skip events where there are large positive amplitudes in the PMT channels (real PMT photoelectron signals are negative, positive signals are pick up noise)
					for (int s = 0; s < binNumber; ++s) {
						if (skip_event_threshold != 0 && (skip_event_threshold > 0 && output_channel > 8 && val >= skip_event_threshold) || skip_event_threshold < 0 && output_channel > 8 && val <= skip_event_threshold) {
							event_flag_cnt++;
						}
					}

					//hCh->SetLineColor(ch + 1); // gets a bit too colorful
					//hCh->SetMarkerColor(ch + 1);

					wfcounter++;
				}//--------------------------------------------------------------------------------------------------------------

			} // for ch

			bool event_flag = false;
			if (event_flag_cnt >= skip_event_threshold_nch) event_flag = true;
			eventnr_storage.push_back(output_event);	//  Adds the current event number(the one from the WaveCatcher) to the storage vector
			skip_event.push_back(event_flag);
			event_counter++;
			if (event_flag) removed_event_counter++;
		} // while an_event

		input_file.close();
		file_counter++;
	} // for file_id

	// in case there are empty channels, nchannels is the number of channels which contain data
	nchannels = read_channels;

	// get bins where the sum spectrum has its maximum for runs with fixed trigger delay and fixed integration window relative to the max of the sum spectrum (not working for DC measurement)
	for (int ch = 0; ch < nchannels; ch++) {
		double max = 0.;
		for (int i = 0; i < binNumber; i++) {
			if (amplValuessum[ch][i] > max) {
				max = amplValuessum[ch][i];
				maxSumBin[ch] = i;
			}
		}
	}

	nevents = event_counter;
	nwf = wfcounter;

	if (skip_event_threshold != 0) cout << "\n\n total number of events: " << nevents << "\n number of removed events: " << removed_event_counter << "\n\n";
}

/// @brief Destructor
ReadRun::~ReadRun() {
	rundata->Clear("C");
	//delete[] maxSumBin;
	//delete baseline_correction_result;
	//plot_active_channels.clear();
	//root_out->Close();
	cout << "\ndeleting nothing currently..." << endl;
}

/// @brief Plot sums of all waveforms for each channel
/// @param doaverage If true it will apply smoothing to plots. \n 
/// Do not use without very good reason as it biases the results.
/// @param normalize If true will normalize the maximum to 1.
/// @param shift Shift histogram by "shift" ns
/// @param sigma Number of bins for running average/gauss sigma in ns for convolution.
/// @param doconv If false use running average (default). \n 
/// If true use gaussian smoothing (slower).
void ReadRun::PlotChannelSums(bool doaverage, bool normalize, double shift, double sigma, bool doconv) {

	double* xv = getx(shift);
	TMultiGraph* mgsums = new TMultiGraph();
	mgsums->SetTitle("channel sums; t [ns]; amplitude [arb.]");

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			double* yv = amplValuessum[i];
			if (doaverage) SmoothArray(yv, binNumber, sigma, doconv);

			TGraph* gr = new TGraph(binNumber, xv, yv);
			delete[] yv;
			if (normalize) gr->Scale(1. / TMath::MaxElement(gr->GetN(), gr->GetY()));

			TString name(Form("channel_%02d", active_channels[i]));
			TString title(Form("Channel %d", active_channels[i]));
			gr->SetName(name.Data());
			gr->SetTitle(title.Data());
			gr->SetLineColor(rcolor(i));
			gr->SetMarkerColor(rcolor(i));
			mgsums->Add(gr);
		}
	}
	delete[] xv;

	TCanvas* sumc = new TCanvas("Sums", "", 1600, 1000);
	mgsums->Draw("APL");
	mgsums->GetYaxis()->SetRangeUser(-1e4, 1e6);
	if (normalize) mgsums->GetYaxis()->SetRangeUser(-0.2, 1);
	sumc->BuildLegend(0.85, 0.70, .99, .95);
	root_out->WriteObject(mgsums, "channelsums");
	root_out->WriteObject(sumc, "channelsums_c");
}

/// @brief Averaging all waveforms (for testing, do not use for analysis!)
/// @param sigma Number of bins for running average/gauss sigma in ns for convolution.
/// @param doconv If false use running average (default). \n 
/// If true use gaussian smoothing (slower).
void ReadRun::SmoothAll(double sigma, bool doconv) { //deprecated since it can be done with baseline correction??
	// just for testing, not very efficient
	cout << "\nsmoothing wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his);
		SmoothArray(yvals, binNumber, sigma, doconv);
		for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

/// @brief derivative of all waveforms (for measurements w/o pole-zero cancellation)
///
/// Experimental!
void ReadRun::DerivativeAll() {
	// just for testing
	cout << "\nderivative of wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his);
		for (int i = 1; i < his->GetNbinsX() - 1; i++) his->SetBinContent(i, yvals[i + 1] - yvals[i]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

/// @brief Baseline correction constant window
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Uses the mean between t=0 and t="tCut" or between t="tCut" and t="tCutEnd" as offset. \n
/// Call method before ReadFile() if you want it to happen while reading. \n \n 
/// 
/// Most simple and fast method. \n 
/// Useful for measurements with very few background events/dark counts. \n 
/// Using a constant window for the baseline means there might be a background pulse which would lead to wrong correction. \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). 
/// 
/// @param tCut Time denoting the end or the beginning (if "tCutEnd" is set) of the integration window.
/// @param tCutEnd Time denoting the end of the integration window.
void ReadRun::CorrectBaseline(float tCut, float tCutEnd) {

	printf("\nUsing method CorrectBaseline\n");
	tCutg = tCut;
	tCutEndg = tCutEnd;
	if (nwf == 0) {
		Using_BaselineCorrection_in_file_loop = true;
	}
	else {
		printf("Baseline correction (%d waveforms) :: ", nwf);
		for (int j = 0; j < nwf; j++) {
			//TH1F* his = ((TH1F*)rundata->At(j));
			CorrectBaseline_function((TH1F*)rundata->At(j), tCut, tCutEnd, j);

			if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
		}
		cout << endl;
	}
}

/// @brief Helper function called by CorrectBaseline()
/// 
/// See CorrectBaseline()
/// 
void ReadRun::CorrectBaseline_function(TH1F* his, float tCut, float tCutEnd, int nwaveform) {
	int iCut, iCutEnd;
	float corr = 0;

	iCut = his->GetXaxis()->FindBin(tCut);

	if (tCutEnd <= 0) { //
		corr = his->Integral(1, iCut) / static_cast<float>(iCut);
	}
	else {
		iCutEnd = his->GetXaxis()->FindBin(tCutEnd);
		corr = his->Integral(iCut, iCutEnd) / static_cast<float>(iCutEnd - iCut + 1);
	}

	// write corrected values to histograms
	if (tCut >= 0) {
		for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, his->GetBinContent(i) - corr);
	}

	if (!Using_BaselineCorrection_in_file_loop) {
		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[nwaveform].push_back(corr);
		baseline_correction_result[nwaveform].push_back(0);
		baseline_correction_result[nwaveform].push_back(tCut);
		baseline_correction_result[nwaveform].push_back(tCutEnd);
	}
}

/// @brief Baseline correction using minimum slope of RMS method
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Determines the region of "nIntegrationWindow" bins where the root mean square of the slope of the smoothed waveform reaches its minimum. \n \n 
/// 
/// Slow but most versatile since it searches for the optimal baseline candidate region in a defined range. \n 
/// Will prefer constant sections of the waveform for the estimation of the baseline. \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n \n 
/// 
/// @param nIntegrationWindow Number of bins of integration window
/// @param doaverage If true will use averaging for more reliable slope. Use with care!
/// @param sigma Number of bins for running average/gauss sigma in ns for convolution. Use with care!
/// @param max_bin_for_baseline Maximum bin for search window.
/// @param start_at Minimum bin for search window.
/// @param search_min Experimental, use with care.
/// @param convolution If false use running average (default). \n 
/// If true use gaussian smoothing (slower).
/// @param skip_channel Skip a channel
/// @todo Work on "skip_channel" and remove "search_min"
void ReadRun::CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel) {

	const int binNumberSlope = binNumber - 1;
	double* slope = new double[binNumberSlope];
	skip_channel += 1;

	if (start_at > max_bin_for_baseline - nIntegrationWindow) start_at = 0;

	int min_distance_from_max = 25 + nIntegrationWindow;
	float corr = 0;
	float minchange = 1.e9;
	float minsum = 0;
	float minsum0 = 0;
	float minsqsum = 0;
	int iintwindowstart = 0;

	float sum = 0.;
	float sum0 = 0.;
	float sqsum = 0.;
	float change = 0.;
	float sign = 1.;

	int imax = 0;
	int search_before = 0;

	printf("\nBaseline correction (%d waveforms) :: ", nwf);

	for (int j = 0; j < nwf; j++) {
		corr = 0;
		minchange = 1.e9;
		iintwindowstart = 0;

		sum = 0.;
		sum0 = 0.;
		sqsum = 0.;
		change = 0.;
		sign = 1.;

		imax = 0;
		search_before = 0;

		if (j == 0 || j != skip_channel - 1 || j % skip_channel != 0) { //eventnr * nchannels + i
			TH1F* his = ((TH1F*)rundata->At(j));
			double* yvals = gety(his); //find faster way
			SmoothArray(yvals, binNumber, sigma, convolution); // smoothing important to suppress variations in slope due to noise so the method is more sensitve to excluding peaks

			//calculate slope
			for (int i = 0; i < binNumberSlope; i++) slope[i] = yvals[i + 1] - yvals[i];

			if (max_bin_for_baseline != 0 && max_bin_for_baseline > nIntegrationWindow) {
				search_before = max_bin_for_baseline - nIntegrationWindow - 1;
			}
			else {
				imax = his->GetMaximumBin();
				search_before = imax - min_distance_from_max;
			}

			for (int i = start_at; i < search_before; i += 3) { // currently in steps of 3 bins (~1 ns) to make it faster
				sum = 0.;
				sum0 = 0.;
				sqsum = 0.;
				change = 0.;
				sign = 1.;
				for (int k = i; k < nIntegrationWindow + i; k += 3) {
					sum += slope[k];
					sum0 += yvals[k] / 100; // completely random choice
					sqsum += (slope[k] * slope[k]);
				}
				if (sum0 < 0) sign = -1.;

				if (search_min) change = sqsum + sum * sum + sum0 * sign;
				else change = sqsum + sum * sum;

				if (change < minchange) {
					minchange = change;
					iintwindowstart = i;
					minsum = sum * sum;
					minsqsum = sqsum;
					minsum0 = sum0 * sign;
				}
			}

			corr = 0.;
			if (!doaverage) {
				corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow + 1);
			}
			else {
				for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
				corr /= static_cast<float>(nIntegrationWindow);
			}

			for (int i = 0; i < binNumber; i++) {
				if (!doaverage) his->SetBinContent(i, his->GetBinContent(i) - corr);
				else his->SetBinContent(i, yvals[i] - corr);
			}
			delete[] yvals; //delete slow
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(minchange);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart) * SP);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart + nIntegrationWindow) * SP);
		baseline_correction_result[j].push_back(minsum);
		baseline_correction_result[j].push_back(minsum0);
		baseline_correction_result[j].push_back(minsqsum);

		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
	delete[] slope;

	//TCanvas* sumc = new TCanvas("sumc", "sumc", 1600, 1000);
	//TH1F* hiss = new TH1F("sum", "sum", 1e4, -2, 2);
	//TH1F* hiss0 = new TH1F("sum0", "sum0", 1e4, -2, 2);
	//TH1F* hisssq = new TH1F("sqsum", "sqsum", 1e4, -2, 2);
	//for (int i = 0; i < nwf; i++) {
	//	hiss->Fill(baseline_correction_result[i][4]);
	//	hiss0->Fill(baseline_correction_result[i][5]);
	//	hisssq->Fill(baseline_correction_result[i][6]);
	//}
	//hiss0->SetLineColor(2); hisssq->SetLineColor(1);
	//hiss->Draw(); hiss0->Draw("same"); hisssq->Draw("same");
}

/// @brief Baseline correction using minimum mean in range for correction 
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Uses min(mean("nIntegrationWindow")) in range {"start_at", "max_bin_for_baseline"} \n \n 
/// 
/// Helpful for (groups of/irradiated) SiPMs with very high dark count rate where the signal never relaxes to the baseline: \f$\Rightarrow DCR \sim 1/\tau_{quenching}\f$. \n 
/// Estimates a baseline as the minimum level before the main peak. \n 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n 
/// For parameters see CorrectBaselineMinSlopeRMS() 
/// 
void ReadRun::CorrectBaselineMin(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool convolution, int skip_channel) {

	int binNumberSlope = binNumber - 1;
	skip_channel += 1;

	if (start_at > max_bin_for_baseline - nIntegrationWindow) start_at = 0;

	int min_distance_from_max = 25 + nIntegrationWindow;

	float corr = 0;
	float minchange = 1.e9;
	int iintwindowstart = 0;
	float sum0 = 0.;
	int imax = 0;
	int search_before = 0;

	printf("\nBaseline correction (%d waveforms) :: ", nwf);

	for (int j = 0; j < nwf; j++) {
		minchange = 1.e9;
		iintwindowstart = 0;

		if (j == 0 || j != skip_channel - 1 || j % skip_channel != 0) { //eventnr * nchannels + i
			TH1F* his = ((TH1F*)rundata->At(j));
			double* yvals = gety(his); //find faster way
			SmoothArray(yvals, binNumber, sigma, convolution); // smoothing

			if (max_bin_for_baseline != 0 && max_bin_for_baseline > nIntegrationWindow) {
				search_before = max_bin_for_baseline - nIntegrationWindow - 1;
			}
			else {
				imax = his->GetMaximumBin();
				search_before = imax - min_distance_from_max;
			}

			for (int i = start_at; i < search_before; i++) { // can also be done in coarser steps
				sum0 = 0.;
				for (int k = i; k < nIntegrationWindow + i; k += 2) { // can also be done in coarser steps
					sum0 += yvals[k];
				}

				if (sum0 < minchange) {
					minchange = sum0;
					iintwindowstart = i;
				}
			}

			corr = 0.;
			if (!doaverage) {
				corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow + 1);
			}
			else {
				for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
				corr /= static_cast<float>(nIntegrationWindow);
			}

			for (int i = 0; i < binNumber; i++) {
				if (!doaverage) his->SetBinContent(i, his->GetBinContent(i) - corr);
				else his->SetBinContent(i, yvals[i] - corr);
			}
			delete[] yvals; //delete slow
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(minchange);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart) * SP);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart + nIntegrationWindow) * SP);

		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

/// @brief Determine the timing of the maximum peak with constant fraction discrimination
/// 
/// Determines timing in the time interval ["start_at_t", "end_at_t"] with CFD for fraction of maximum "cf_r". \n 
/// Stores timing information for all channels and all events in ReadRun::timing_results. \n 
/// Per event results will be visualized in PrintChargeSpectrumWF(). \n 
/// Cumulative per channel results can be visualized with Print_GetTimingCFD(). \n 
/// Results can be used for time difference cuts with SkipEventsTimeDiffCut().
/// 
/// @param cf_r Fraction of maximum for CFD.
/// @param start_at_t Time in ns to start searching.
/// @param end_at_t Time in ns to end searching.
/// @param sigma Should be set to 0. \n 
/// If >0 will use running average for CFD. This will bias the results! Do not use (or use very carefully, only for noisy data)!
/// @param find_CF_from_start If true will start search from "start_at_t" \n 
/// If false searches backwards from time of maximum (default setting).
void ReadRun::GetTimingCFD(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start) {

	int start_at = static_cast<int>(start_at_t / SP);
	int end_at = static_cast<int>(end_at_t / SP);
	int n_range = end_at - start_at;

	printf("\nGet timing at CF=%.2f between %.2f ns and %.2f ns (%d waveforms) :: ", cf_r, start_at_t, end_at_t, nwf);

	for (int j = 0; j < nwf; j++) {

		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his, start_at, end_at); // get range where to search for CFD for timing

		if (sigma > 0.) SmoothArray(yvals, n_range, sigma); // smoothing to suppress noise, will also change timing so use with care!

		float max = 0.;
		int n_max = 0;
		int i = 0;
		for (i = 0; i < n_range; i++) {
			if (yvals[i] > max) {
				max = yvals[i];
				n_max = i;
			}
		}

		float cf = cf_r * max;
		if (!find_CF_from_start) {
			i = n_max;
			while (yvals[i] > cf && i >= 0) i--;
			i++;
		}
		else {
			i = 0;
			while (yvals[i] < cf && i <= n_max) i++;
			i--;
		}
		float interpol_bin = .0;
		interpol_bin = LinearInterpolation(cf, static_cast<float>(i), static_cast<float>(i+1), yvals[i], yvals[i+1]);

		timing_results.push_back(vector<float>());
		timing_results[j].push_back(interpol_bin); // return the bin we looked for
		timing_results[j].push_back(static_cast<float>(interpol_bin + static_cast<float>(start_at)) * SP); // return the cfd-time we looked for
		timing_results[j].push_back(static_cast<float>(start_at) * SP); // return starting time
		timing_results[j].push_back(static_cast<float>(end_at) * SP); // return the end time
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
	cout << endl;
}

/// @brief Skip events where the time difference between two channels is outside of specified range
/// 
/// Skip events where the time difference between channel "first_channel" and channel "second_channel" 
/// is less than "time_diff_min" or more than "time_diff_max" \n 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param first_channel_abs First channel.
/// @param second_channel_abs Second channel.
/// @param time_diff_min Skip events for \f$\Delta t<t_{diff,min}\f$
/// @param time_diff_max Skip events for \f$\Delta t>t_{diff,max}\f$
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsTimeDiffCut(int first_channel_abs, int second_channel_abs, double time_diff_min, double time_diff_max, bool verbose) {

	printf("\n\n Removing events if the event-wise time difference between the main peaks in ch%d and ch%d is <%.2f ns or >%.2f ns\n\n", first_channel_abs, second_channel_abs, time_diff_min, time_diff_max);
	int counter = 0;

	// match channel number to channel index
	int first_channel = 0;
	int second_channel = 0;
	for (int i = 0; i < active_channels.size(); i++) {
		if (active_channels[i] == first_channel_abs) first_channel = i;
		if (active_channels[i] == second_channel_abs) second_channel = i;
	}

	// call GetTimingCFD() in case it was not initialized
	if (timing_results.size() == 0) GetTimingCFD();

	// loop through events, calculate timing difference between channels and compare with cuts
	for (int j = 0; j < nwf; j += nchannels) {
		if (!skip_event[floor(j / nchannels)]) {
			float time_diff = fabs(timing_results[j + first_channel][1] - timing_results[j + second_channel][1]);

			if (j <= timing_results.size() && (time_diff < time_diff_min || time_diff > time_diff_max)) {
				int currevent = eventnr_storage[floor(j / nchannels)];
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannels:\t" << first_channel_abs << " & " << second_channel_abs << "\ttime diff:\t" << time_diff;
				skip_event[floor(j / nchannels)] = true;
				counter++;
			}
		}
	}
	cout << "\n\n\t" << counter << " events will be cut out of " << nevents << "\n\n";
}

/// @brief Find events with max/min above/below a certain threshold
/// 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param threshold Threshold in mV.
/// @param max If true uses max, else uses min.
/// @param greater If true looks for events with max/min>threshold, else looks for events with max/min<threshold.
/// @param from Start search at "from" in ns.
/// @param to End search at "to" in ns.
/// @param verbose Set true for extra verbosity.
void ReadRun::FractionEventsAboveThreshold(float threshold, bool max, bool greater, double from, double to, bool verbose) {

	int occurences = 0;
	int occurences2ch = 0;
	int o2ch = 0;
	int currchannel = 0;
	int currevent = 0;
	int lastevent = 0;
	if (plot_active_channels.empty()) plot_active_channels = active_channels;
	vector<int> counter_abovethr(plot_active_channels.size());	// DORAMAS: It stores a counter of events above threshold for each channel that will be plotted

	cout << "\n\n ------> ";
	if (max) cout << "max";
	else cout << "min";

	if (greater) cout << " > ";
	else cout << " < ";
	cout << threshold << " mV:\n";

	for (int j = 0; j < nwf; j++) {
		auto his = (TH1F*)((TH1F*)rundata->At(j))->Clone(); // use Clone() to not change ranges of original histogram

		// set range where to search for amplitudes above threshold
		if (from >= 0 && to > 0) {
			his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
		}

		currchannel = j - nchannels * floor(j / nchannels);
		if (find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[currchannel]) != plot_active_channels.end()) {
			if ((max && greater && his->GetMaximum() > threshold) || (max && !greater && his->GetMaximum() < threshold) || (!max && greater && his->GetMinimum() > threshold) || (!max && !greater && his->GetMinimum() < threshold)) {
				currevent = eventnr_storage[floor(j / nchannels)];

				if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel];

				// We must use 'distance' to make sure the position in 'counter_above' matches with the corresponding channel's position at 'plot_active_channels'
				counter_abovethr[distance(plot_active_channels.begin(), find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[currchannel]))] += 1;
				// This is to detect events w/ at least two channels above threshold
				if (lastevent != currevent) occurences += 1;
				if (lastevent == currevent && o2ch != occurences) {
					occurences2ch += 1;
					o2ch = occurences;
				}
				lastevent = currevent;
			}
		}
	}

	if (verbose) cout << endl;

	//  Loop to show the fraction of events above threshold for each channel that will be plotted
	for (int i = 0; i < plot_active_channels.size(); i++) {
		cout << "\nfraction of events in channel " << plot_active_channels[i] << " above threshold: " << 100. * static_cast<float>(counter_abovethr[i]) / static_cast<float>(nevents) << "%\n";
	}
	//
	cout << "\nfraction of events w/ at least 2 channels above threshold: " << 100. * static_cast<float>(occurences2ch) / static_cast<float>(nevents) << "%\n";
	cout << "\tfor a total of " << nevents << " events\n" << endl;
}

/// @brief Skip events above/below indiviual thresholds per channel
/// 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). Negative thresholds mean events below threshold will be cut.
/// @param rangestart Range start in ns
/// @param rangeend Range end in ns
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsPerChannel(vector<double> thresholds, double rangestart, double rangeend, bool verbose) { // merge with IntegralFilter()?

	cout << "\n\n Removing events with individual amplitude threshold per channel!!!\n\n";
	int counter = 0;

	for (int j = 0; j < nwf; j++) {
		if (!skip_event[floor(j / nchannels)]) {
			auto his = (TH1F*)((TH1F*)rundata->At(j))->Clone(); // use Clone() to not change ranges of original histogram
			if (rangestart != 0 && rangeend != 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(rangeend));

			int currchannel = j - nchannels * floor(j / nchannels);
			if (currchannel <= thresholds.size() && thresholds[currchannel] != 0 && !skip_event[floor(j / nchannels)] && ((thresholds[currchannel] > 0 && his->GetMaximum() > thresholds[currchannel]) || (thresholds[currchannel] < 0 && his->GetMinimum() < thresholds[currchannel]))) {

				int currevent = eventnr_storage[floor(j / nchannels)];
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\tthreshold\t" << thresholds[currchannel];
				skip_event[floor(j / nchannels)] = true;
				counter++;
			}
		}
	}

	cout << "\n\n\t" << counter << " events will be cut out of " << nevents << "\n\n";
}

/// @brief Skip events with threshold on integral
/// 
/// Compare with SkipEventsPerChannel()
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// Please note that the old version was missing start and end, so please add these parameters to older scripts if they throw errors.
/// 
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). Negative thresholds mean events below threshold will be cut.
/// @param highlow  Vector should contain a bool for each active channel. True means events with integrals above threshold will be cut, false means below threshold.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param use_AND_condition If set true it will overrule previously applied cuts and unskip events that pass integral criterium.
/// @param verbose Set true for extra verbosity.
void ReadRun::IntegralFilter(vector<double> thresholds, vector<bool> highlow, float windowlow, float windowhi, float start, float end, bool use_AND_condition, bool verbose) {

	cout << "\n\n Removing events with individual integral threshold per channel!!!\n\n";
	int counter = 0;
	float integral = 0;

	for (int j = 0; j < nwf; j++) {
		int currevent_counter = floor(j / nchannels);

		if (use_AND_condition || !skip_event[currevent_counter]) {
			auto his = (TH1F*)((TH1F*)rundata->At(j))->Clone(); // use Clone() to not change ranges of original histogram

			int currchannel = j - nchannels * currevent_counter;
			integral = GetPeakIntegral(his, windowlow, windowhi, start, end, currchannel);

			if (currchannel <= thresholds.size() && thresholds[currchannel] != 0 && !skip_event[currevent_counter] && ((highlow[currchannel] && integral > thresholds[currchannel]) || (!highlow[currchannel] && integral < thresholds[currchannel]))) {
				int currevent = eventnr_storage[currevent_counter];
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\tthreshold\t" << thresholds[currchannel];
				skip_event[currevent_counter] = true;
				while (floor((j + 1) / nchannels) == currevent_counter) j++;
				counter++;
			}
			else if (use_AND_condition && thresholds[currchannel] != 0 && skip_event[currevent_counter] && !highlow[currchannel]) {
				int currevent = eventnr_storage[currevent_counter];
				skip_event[currevent_counter] = false;
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\thas been flagged good by integral";
			}
		}
	}
	cout << "\n\n\t" << counter << " events will be cut out of " << nevents << endl;
}

/// @brief Prints a list of all skipped events into the terminal for diagnostics
void ReadRun::PrintSkippedEvents() {
	int counter = 0;
	for (int j = 0; j < skip_event.size(); j++) {
		if (skip_event[j]) {
			int currevent = eventnr_storage[j];
			cout << "\nevent:\t" << currevent;
			counter++;
		}
	}
	cout << "\n\ntotal number of skipped events:\t" << counter << "\tout of:\t" << nevents << endl;
}

// functions for charge spectrum

/// @brief Determine indices for integration window for peaks
/// 
/// Default usage: Find maximum in range ("start", "end") and return bin numbers for [0] the max, [1] t_max - "windowlow", and [2] t_max + "windowhi" \n 
/// If ("start" < 0 || "end" < 0) doesn't return max and integration window is fixed t(max(sum_spectrum[channel])) +/- "windowhi"/"windowlow" \n 
/// If ("windowlow" == "start" && "windowwhi" == "end") doesn't return max and sets fixed integration window from "start" until "end" for all channels.
/// 
/// @param his Histgrom to integrate.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel Channel index in case the integration should be around the maximum of the sum of all waveforms
/// @return Array { max, \f$n_{t,start}\f$, \f$n_{t,end}\f$ }
int* ReadRun::GetIntWindow(TH1F* his, float windowlow, float windowhi, float start, float end, int channel) {

	int istart, iend;
	int* foundindices = new int[3];//
	foundindices[0] = 0;

	if (start < 0 || end < 0) {									// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(maxSumBin[channel]) - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(maxSumBin[channel]) + windowhi);
	}
	else if (windowlow == start && windowhi == end) {				// fixed integration window for all channels
		foundindices[1] = his->GetXaxis()->FindBin(windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(windowhi);
	}
	else {															// fixed integration window relative to maximum of each individual waveform
		istart = his->GetXaxis()->FindBin(start);
		iend = his->GetXaxis()->FindBin(end);
		foundindices[0] = istart;

		if (istart<1 || iend>his->GetNbinsX()) {
			cout << "\nError: start or end out of range" << endl;
			return 0;
		}

		float max = -9.e99;
		float val = 0;
		for (int i = istart; i < iend; i++) {
			val = his->GetBinContent(i);
			if (val > max) {
				max = val;
				foundindices[0] = i;
			}
		}

		foundindices[1] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(foundindices[0]) - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(his->GetXaxis()->GetBinCenter(foundindices[0]) + windowhi);
	}
	return foundindices;
}

/// @brief Calculate the integral around a peak with several options explained in GetIntWindow().
/// @param his Histogram with peak.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel_index Channel index in case the integration should be around the maximum of the sum of all waveforms.
/// @return Integral/amplitude.
float ReadRun::GetPeakIntegral(TH1F* his, float windowlow, float windowhi, float start, float end, int channel_index) {
	int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, channel_index);	// find integration window
	string integral_option(""); // For amplitude -> unit[mV].
	if (windowind[1] != windowind[2]) integral_option = "width"; // 'width' (bin width) for integral -> unit[mV x ns].
	float integral = his->Integral(windowind[1], windowind[2], integral_option.c_str());
	delete[] windowind;
	return integral;
}

/// @brief Plot waveforms of all channels for a given event number and add the determined integration windows to the plot
/// 
/// See GetIntWindow() for explanation of parameters.
/// Will also add CFD timing if GetTimingCFD() was called before.
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param eventnr Event number
/// @param ymin Y axis minimum
/// @param ymax Y axis maximum
void ReadRun::PrintChargeSpectrumWF(float windowlow, float windowhi, float start, float end, int eventnr, float ymin, float ymax) {

	gStyle->SetOptStat(0);

	TString name(Form("waveforms_event__%05d", eventnr));
	TCanvas* intwinc = new TCanvas(name.Data(), name.Data(), 1600, 1000);
	SplitCanvas(intwinc);
	int event_index = GetEventIndex(eventnr);

	int current_canvas = 0;
	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his = Getwf(i, event_index);
			int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, i);
			// create lines to indicate the integration window
			TLine* low = new TLine(his->GetXaxis()->GetBinCenter(windowind[1]), -5, his->GetXaxis()->GetBinCenter(windowind[1]), 10);
			low->SetLineColor(2);
			TLine* hi = new TLine(his->GetXaxis()->GetBinCenter(windowind[2]), -2, his->GetXaxis()->GetBinCenter(windowind[2]), 3);
			hi->SetLineColor(2);
			TLine* zero = new TLine(0, 0, 320, 0); // draw line at x=0 to check if baseline correction worked
			zero->SetLineColor(1);
			delete[] windowind;

			// draw to canvas
			intwinc->cd(current_canvas);
			his->Draw();
			if (ymin != 0. && ymax != 0.) his->GetYaxis()->SetRangeUser(ymin, ymax); //for better comparison fix y range
			low->Draw("same");
			hi->Draw("same");
			zero->Draw("same");

			// draw baseline parameters
			if (baseline_correction_result.size() > event_index * nchannels + i) {
				TLine* baselinel = new TLine(baseline_correction_result[event_index * nchannels + i][2], -1, baseline_correction_result[event_index * nchannels + i][2], 1);
				baselinel->SetLineColor(6);
				baselinel->SetLineWidth(2);
				TLine* baselineh = new TLine(baseline_correction_result[event_index * nchannels + i][3], -1, baseline_correction_result[event_index * nchannels + i][3], 1);
				baselineh->SetLineColor(6);
				baselineh->SetLineWidth(2);
				TLine* baseline = new TLine(baseline_correction_result[event_index * nchannels + i][2], 0, baseline_correction_result[event_index * nchannels + i][3], 0);
				baseline->SetLineColor(6);
				TLine* correction_value = new TLine(baseline_correction_result[event_index * nchannels + i][2], baseline_correction_result[event_index * nchannels + i][0], baseline_correction_result[event_index * nchannels + i][3], baseline_correction_result[event_index * nchannels + i][0]);
				correction_value->SetLineColor(2);

				baselinel->Draw("same");
				baselineh->Draw("same");
				baseline->Draw("same");
				correction_value->Draw("same");
			}

			if (timing_results.size() > event_index * nchannels + i) {
				TLine* timing = new TLine(timing_results[event_index * nchannels + i][1], -10, timing_results[event_index * nchannels + i][1], 100);
				timing->SetLineColor(9);
				timing->SetLineWidth(2);
				timing->Draw("same");
			}
		}
	}
	intwinc->Update();

	root_out->WriteObject(intwinc, name.Data());
}

/// @brief Histogram of the "charge" spectrum for one channel
/// 
/// See PrintChargeSpectrum() for parameters.
/// 
/// @return Histogram for one channel.
TH1F* ReadRun::ChargeSpectrum(int channel_index, float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {

	TString name(Form("channel__%02d", active_channels[channel_index]));
	TH1F* h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	for (int j = 0; j < nevents; j++) {
		if (!skip_event[j]) {
			TH1F* his = ((TH1F*)rundata->At(j * nchannels + channel_index));
			h1->Fill(GetPeakIntegral(his, windowlow, windowhi, start, end, channel_index)); // fill charge spectrum
		}
	}
	return h1;
}

/// @brief Plots the "charge" spectrums of all channels
/// 
/// Integrate all pulses in range ("start", "end") from t_max - "windowlow" to t_max + "windowhi" for a given channel 
/// and return the charge histogram with x range ("rangestart", "rangeend") and the number of bins "nbins". \n 
/// It is not really charge, but either amplitude (mV) or integral (mV x ns).
/// See ChargeSpectrum() and GetIntWindow().
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param rangestart Plot x range start
/// @param rangeend Plot x range end
/// @param nbins Number of bins of histogram
/// @param fitrangestart Fit range start
/// @param fitrangeend Fit range end
/// @param max_channel_nr_to_fit Fit only channels with index > "max_channel_nr_to_fit". Set to -1 to skip fitting.
/// @param which_fitf Choose fit function: \n 
/// 1 - landau gauss convolution for large number of photons \n 
/// 2 - if pedestal is biased because of peak finder algorithm \n 
/// 3 - SiPM fit function with exponential delayed afterpulsing \n 
/// 4 - ideal PMT fit function \n 
/// 5 - PMT fit function \n 
/// 6 - PMT fit function with biased pedestal \n 
/// 7 - default SiPM fit function + dark count spectrum (for lots of false triggers) \n 
/// else - default SiPM fit function \n 
void ReadRun::PrintChargeSpectrum(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf) {
	// print ReadRun::ChargeSpectrum for all channels optimized for SiPM signals
	PrintChargeSpectrum_cnt++;

	gStyle->SetOptStat("ne");
	gStyle->SetOptFit(1111);

	if (fitrangestart == 0.) fitrangestart = rangestart;
	if (fitrangeend == 0.) fitrangeend = rangeend;

	string ctitle("\"charge\" spectra" + to_string(PrintChargeSpectrum_cnt));
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 1600, 1000);
	SplitCanvas(chargec);

	cout << "\n\nThere is data recorded in " << active_channels.size() << " channels \n\n\n";
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, start, end, rangestart, rangeend, nbins);
			his->GetYaxis()->SetTitle("#Entries");
			if (windowlow + windowhi > 0.) his->GetXaxis()->SetTitle("integral in mV#timesns");
			else his->GetXaxis()->SetTitle("amplitude in mV");

			TString name(Form("ChargeSpectrum channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name.Data());

			chargec->cd(current_canvas);

			//do the fit
			if (which_fitf == 1) { // landau gauss convolution for large number of photons
				Fitf_langaus fitf;
				TF1* f = new TF1("fitf_langaus", fitf, fitrangestart, fitrangeend, 4); f->SetLineColor(3);

				f->SetParName(0, "Width");				f->SetParameter(0, 35);
				f->SetParName(1, "MPV");				f->SetParameter(1, 1000);
				f->SetParName(2, "Area");			    f->SetParameter(2, 10000);
				f->SetParName(3, "#sigma_{Gauss}");		f->SetParameter(3, 100);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else if (which_fitf == 2) { // if pedestal is biased because of peak finder algorithm
				Fitf_biased fitf;
				TF1* f = new TF1("fitf_biased", fitf, fitrangestart, fitrangeend, 9); f->SetLineColor(3);

				f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
				f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
				f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
				f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4); //f->SetParLimits(4, 1.e-9, 1.e3);
				f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 10.);
				f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
				f->SetParName(7, "norm_{0}");			f->SetParameter(7, 0.7);
				f->SetParName(8, "x_{0}");				f->SetParameter(8, 5.);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}

				// get number of excess events in the pedestal in the fit region. To get the absolute number of excess events the full pedestal needs to be inside of the fit range (fitrangestart, fitrangeend)
				//double excessEventsInPedestal = f->Integral(fitrangestart, fitrangeend)/.3125;
				//f->SetParameter(7, 1.);
				//excessEventsInPedestal -= f->Integral(fitrangestart, fitrangeend)/.3125;
				//cout << "\nNumber of excess events in the pedestal within the fit range:\t" << excessEventsInPedestal << "\n\n";
			}
			else if (which_fitf == 3) { // SiPM fit function with exponential delayed afterpulsing
				Fitf_full fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 9); f->SetLineColor(3);

				f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
				f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
				f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
				f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
				f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
				f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
				f->SetParName(7, "#alpha");				f->SetParameter(7, .1); //f->FixParameter(7, .2);
				f->SetParName(8, "#beta");				f->SetParameter(8, 80.); //f->FixParameter(8, 80);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else if (which_fitf == 4) { // ideal PMT fit function
				Fitf_PMT_ideal fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 4); f->SetLineColor(3);

				f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
				f->SetParName(1, "#mu");				f->SetParameter(1, 1.);
				f->SetParName(2, "#sigma");				f->SetParameter(2, 5.);
				f->SetParName(3, "gain");				f->SetParameter(3, 10.);

				f->SetLineColor(2);
				f->SetNpx(1000);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else if (which_fitf == 5) { // PMT fit function
				Fitf_PMT fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 8); f->SetLineColor(3);

				f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
				f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-99, 4.e-1); //probability for type II BG
				f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-99, 5.e-2); //coefficient of exponential decrease of typ II BG
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
				f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);
				f->SetParName(5, "#mu");				f->SetParameter(5, 1.);
				f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
				f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);

				f->SetLineColor(2);
				f->SetNpx(1000);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else if (which_fitf == 6) { // PMT fit function with biased pedestal
				Fitf_PMT_pedestal fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 9); f->SetLineColor(3);

				f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
				f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
				f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
				f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
				f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
				f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
				f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);	f->SetParLimits(7, 1.e-9, 1.e9);
				f->SetParName(8, "A_{0}");				f->SetParameter(8, 1.);		f->SetParLimits(8, 1.e-9, 1.e1);

				f->SetLineColor(2);
				f->SetNpx(1000);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else if (which_fitf == 7) { // default SiPM fit function + dark count spectrum (for lots of false triggers)
				Fitf_plus_DC fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 9); f->SetLineColor(3);

				f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
				f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
				f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
				f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
				f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
				f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
				f->SetParName(7, "#mu_darkcount");		f->SetParameter(7, .1);		f->SetParLimits(7, 1.e-9, 1.);
				f->SetParName(8, "N_{0}_darkcount");	f->SetParameter(8, .05);	f->SetParLimits(8, 1.e-9, .3);

				f->SetLineColor(2);
				f->SetNpx(1000);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			else { // default SiPM fit function
				Fitf fitf;
				TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 7); f->SetLineColor(3);

				f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
				f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
				f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
				f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
				f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
				f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
				f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);

				if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < PrintChargeSpectrum_pars.size(); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

				if (i < max_channel_nr_to_fit) {
					cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
					TFitResultPtr fresults = his->Fit(f, "LRS");
					fit_results.push_back(fresults);
				}
			}
			his->Draw();
		}
	}

	chargec->Update();
	root_out->WriteObject(chargec, "ChargeSpectra");
}

/// @brief "Charge" spectrum optimized for PMT signals
/// 
/// Just for plotting. To analyze the data use PrintChargeSpectrum() with Fitf_PMT_pedestal() for low number of photons and Fitf_langaus() for >10-15 photons. \n 
/// See PrintChargeSpectrum() for parameters.
/// 
void ReadRun::PrintChargeSpectrumPMT(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {

	PrintChargeSpectrumPMT_cnt++;

	string ctitle("charge spectra PMT" + to_string(PrintChargeSpectrumPMT_cnt));
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 1600, 1000);
	SplitCanvas(chargec);

	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, start, end, rangestart, rangeend, nbins);
			if (windowlow + windowhi > 0.) his->GetXaxis()->SetTitle("integral in mV#timesns");
			else his->GetXaxis()->SetTitle("amplitude in mV");
			chargec->cd(current_canvas);

			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("integral in mV#timesns");
			his->Draw();
			stringstream allname; allname << his->GetEntries() << " entries";
			his->SetTitle(allname.str().c_str());

			TString name(Form("ChargeSpectrumPMT channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name.Data());

			// fit sum of two gaussians
			// estimate starting values of fit by fitting only one gauss
			auto gauss = new TF1("gauss", "gaus", rangestart, rangeend);
			TFitResultPtr fres_est = his->Fit(gauss, "QRSL");
			// now use these fit results to fit the sum of two gauss
			auto two_gauss = new TF1("two gaussians", "gaus(0)+gaus(3)", rangestart, rangeend); two_gauss->SetTitle("Sum of two gauss");
			two_gauss->SetParameters(fres_est->Parameter(0) * .95, fres_est->Parameter(1) * .95, fres_est->Parameter(2) * .95, fres_est->Parameter(0) * .3, fres_est->Parameter(1) * 1.05, fres_est->Parameter(2) * .85); // factors are pretty much random
			two_gauss->SetParName(0, "A_{pedestal}");
			two_gauss->SetParName(1, "#mu_{pedestal}");
			two_gauss->SetParName(2, "#sigma_{pedestal}");
			two_gauss->SetParName(3, "A_{SPE}");
			two_gauss->SetParName(4, "#mu_{SPE}");
			two_gauss->SetParName(5, "#sigma_{SPE}");

			if (!PrintChargeSpectrumPMT_pars.empty()) {
				for (int j = 0; j < PrintChargeSpectrumPMT_pars.size(); j++) two_gauss->SetParameter(j, PrintChargeSpectrumPMT_pars[j]);
			}

			//two_gauss->SetLineColor(4);
			TFitResultPtr fresults = his->Fit(two_gauss, "RSL");
			fit_results.push_back(fresults);

			two_gauss->Draw("same");

			auto pedestal = new TF1("pedestal", "gaus", rangestart, rangeend); pedestal->SetTitle("pedestal");
			pedestal->SetParameters(fresults->Parameter(0), fresults->Parameter(1), fresults->Parameter(2));
			pedestal->SetLineColor(3);
			pedestal->Draw("same");

			auto pepeak = new TF1("pepeak", "gaus", rangestart, rangeend); pepeak->SetTitle("pepeak");
			pepeak->SetParameters(fresults->Parameter(3), fresults->Parameter(4), fresults->Parameter(5));
			pepeak->SetLineColor(4);
			pepeak->Draw("same");

			gPad->BuildLegend();
		}
	}

	chargec->Update();
	root_out->WriteObject(chargec, "ChargeSpectraPMT");
}

/// @brief Print "charge" spectrum with highlight a threshold.
/// 
/// Can also be used to determine the dark count rate of SiPMs.
/// See PrintChargeSpectrum() for parameters.
/// 
/// @param threshold Threshold
/// @param calculate_SiPM_DCR Set true to calculate the fraction above and below threshold and the rate (SiPM DCR).
void ReadRun::PrintChargeSpectrumPMTthreshold(float windowlow, float windowhi, float rangestart, float rangeend, int nbins, double threshold, bool calculate_SiPM_DCR) {

	PrintChargeSpectrumPMTthreshold_cnt++;

	gStyle->SetOptStat(0); // 11 is title + entries

	// show fraction of events above 0.5 pe charge = pedestal + gain/2
	// dark count rate for SiPMs (currently only automated for fit function Fitf)
	// need to call the SiPM fitfunction before this one for this functionality
	bool use_fit_result_for_threshold = false;
	if (threshold == 999) {
		calculate_SiPM_DCR = true;
		use_fit_result_for_threshold = true;
	}

	string unit(" mV");
	string title("amplitude in mV"); // amplitude spectrum not good for fitting, will be biased
	if (windowlow + windowhi > 0.) {
		unit = " mV#timesns";
		title = "integral in mV#timesns";
	}
	string ctitle("charge spectra PMT threshold" + to_string(PrintChargeSpectrumPMTthreshold_cnt));
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 1600, 1000);
	SplitCanvas(chargec);

	int current_canvas = 0;
	double threshold_bin_center = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {

			if (use_fit_result_for_threshold) threshold = fit_results[current_canvas]->Parameter(6) + fit_results[current_canvas]->Parameter(5) / 2.;

			chargec->cd(++current_canvas);

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, rangestart, rangeend, rangestart, rangeend, nbins);
			his->GetXaxis()->SetTitle(title.c_str());
			his->GetYaxis()->SetTitle("#Entries");
			his->Draw();
			stringstream allname; allname << his->GetEntries() << " entries";
			his->SetTitle(allname.str().c_str());

			auto his_lo = (TH1F*)his->Clone();
			his_lo->GetXaxis()->SetRange(his_lo->GetXaxis()->FindBin(rangestart), his_lo->GetXaxis()->FindBin(threshold));
			his_lo->SetLineColor(2);
			his_lo->SetFillColor(2);
			his_lo->Draw("LF2 same");

			stringstream lonamefrac;
			stringstream lonamerate;
			lonamefrac << 100. * his->Integral(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(threshold)) / his->GetEntries() << "% <= " << threshold << unit;
			lonamerate << "<0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(threshold)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << lonamerate.str().c_str() << endl;
			cout << "\n" << lonamefrac.str().c_str() << endl;
			his_lo->SetTitle(lonamerate.str().c_str());
			if (!calculate_SiPM_DCR) his_lo->SetTitle(lonamefrac.str().c_str());

			auto his_hi = (TH1F*)his->Clone();
			his_hi->GetXaxis()->SetRange(his_hi->GetXaxis()->FindBin(threshold), his_lo->GetXaxis()->FindBin(rangeend));
			his_hi->SetLineColor(1);
			his_hi->SetFillColor(1);
			his_hi->Draw("LF2 same");

			stringstream hinamefrac;
			stringstream hinamerate;
			hinamefrac << 100. * his->Integral(his->GetXaxis()->FindBin(threshold) + 1, his->GetXaxis()->FindBin(rangeend)) / his->GetEntries() << "% > " << threshold << unit;
			hinamerate << ">0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(threshold) + 1, his->GetXaxis()->FindBin(rangeend)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << hinamerate.str().c_str() << endl;
			cout << "\n" << hinamefrac.str().c_str() << endl;
			his_hi->SetTitle(hinamerate.str().c_str());
			if (!calculate_SiPM_DCR) his_hi->SetTitle(hinamefrac.str().c_str());

			gPad->BuildLegend();

			threshold_bin_center = his->GetXaxis()->GetBinCenter(his->GetXaxis()->FindBin(threshold) + 1);
			cout << "\n PMT charge spectrum is counting events above threshold from bin center >= " << threshold_bin_center << unit << " for a threshold setting of " << threshold << unit << "\n\n";
		}
	}

	chargec->Update();
	root_out->WriteObject(chargec, "ChargeSpectraPMTthreshold");
}

/// @brief Calculate SiPM dark count rate
/// 
/// See PrintChargeSpectrum() for parameters.
/// 
/// @param threshold Threshold
void ReadRun::PrintDCR(float windowlow, float windowhi, float rangestart, float rangeend, double threshold) {

	string unit(" mV");
	if (windowlow + windowhi > 0.) unit = " mV*ns";

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, rangestart, rangeend, rangestart, rangeend, 500);

			stringstream lonamerate;
			lonamerate << "<0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(threshold)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << lonamerate.str().c_str() << endl;

			stringstream hinamerate;
			hinamerate << ">0.5 pe=" << threshold << unit << " -> " << his->Integral(his->GetXaxis()->FindBin(threshold) + 1, his->GetXaxis()->FindBin(rangeend)) / his->GetEntries() / (1.e-3 * (rangeend - rangestart)) << " MHz";
			cout << "\n" << hinamerate.str().c_str() << endl;
		}
	}
}

/// @brief Time distribution of maximum, CFD, or 10% - 90% rise time in a certain time window
/// 
/// See PrintTimeDist() for parameters.
/// 
TH1F* ReadRun::TimeDist(int channel_index, float from, float to, float rangestart, float rangeend, int nbins, int which, float cf_r) {

	TString name(Form("timedist_ch%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	for (int j = 0; j < nevents; j++) {
		if (!skip_event[j]) {
			auto his = (TH1F*)((TH1F*)rundata->At(j * nchannels + channel_index))->Clone();
			if (from >= 0 && to > 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));

			int from_n = his->GetXaxis()->FindBin(from);
			int to_n = his->GetXaxis()->FindBin(to);

			if (which == 0) { // time of maximum 
				h1->Fill(his->GetXaxis()->GetBinCenter(his->GetMaximumBin()));
			}
			else if (which == 1) { // time of 50% CFD
				float max = his->GetMaximum();
				int max_n = his->GetMaximumBin();

				do {
					max_n--;
				} while (his->GetBinContent(max_n) >= cf_r * max && max_n > from_n);
				max_n++;

				h1->Fill(LinearInterpolation(cf_r * max, his->GetXaxis()->GetBinCenter(max_n - 1), his->GetXaxis()->GetBinCenter(max_n), his->GetBinContent(max_n - 1), his->GetBinContent(max_n)));
			}
			else { // 10%-90% rise time
				// todo -> search backwards from maximum
				double max = his->GetMaximum();
				int n10 = -1;
				int n90 = -1;
				do {
					from_n++;
					if (n10 == -1 && his->GetBinContent(from_n) >= .1 * max && his->GetBinContent(from_n - 1) <= .1 * max) n10 = from_n;
					if (n90 == -1 && his->GetBinContent(from_n) >= .9 * max && his->GetBinContent(from_n - 1) <= .9 * max) n90 = from_n;
				} while (his->GetBinContent(from_n) <= max && from_n < to_n);

				float t10 = LinearInterpolation(.1 * max, his->GetXaxis()->GetBinCenter(n10 - 1), his->GetXaxis()->GetBinCenter(n10), his->GetBinContent(n10 - 1), his->GetBinContent(n10));
				float t90 = LinearInterpolation(.9 * max, his->GetXaxis()->GetBinCenter(n90 - 1), his->GetXaxis()->GetBinCenter(n90), his->GetBinContent(n90 - 1), his->GetBinContent(n90));

				h1->Fill(t90 - t10);
			}
		}
	}
	if (which == 1) h1->Fit("gaus", "L", "same");
	return h1;
}

/// @brief Time distribution of maximum, CFD, or 10% - 90% rise time in a certain time window
/// 
/// Find peak time for a given channel in time window ["from", "to"] and return the peak time histogram 
/// with x range ["rangestart", "rangeend"] and the number of bins "nbins". \n 
/// Plots TimeDist() for all channels. \n \n  
/// 
/// For CFD it is advised to use GetTimingCFD() with Print_GetTimingCFD() instead of this function.
/// 
/// @param from Start of time interval in ns
/// @param to End of time interval in ns
/// @param rangestart Start of x range of histogram
/// @param rangeend End of x range of histogram
/// @param nbins Number of bins of histogram
/// @param which Options: \n 
/// 0 - Gives time of maximum. \n 
/// 1 - Gives constant fraction discrimination with fraction "cf_r" of maximum, searching backwards from the maximum until "from". \n 
/// else - gives the 10% - 90% rise time.
/// @param cf_r Fraction of max for CFD.
void ReadRun::PrintTimeDist(float from, float to, float rangestart, float rangeend, int nbins, int which, float cf_r) {
	// print ReadRun::TimeDist for all channels
	// see TimeDist for parameter explanations

	gStyle->SetOptStat(1111); // 11 is title + entries

	TCanvas* time_dist_c = new TCanvas("timing of maximum", "timing of maximum", 600, 400);
	SplitCanvas(time_dist_c);

	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;
			time_dist_c->cd(current_canvas);

			TH1F* his;
			his = TimeDist(i, from, to, rangestart, rangeend, nbins, which, cf_r);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("time [ns]");
			his->Draw();
			stringstream name; name << "t_{max} for " << from << "<t<" << to << " ns";
			his->SetTitle(name.str().c_str());

			TString name_save(Form("TimeDist channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name_save.Data());
		}
	}

	time_dist_c->Update();
	root_out->WriteObject(time_dist_c, "TimeDist");
}

/// @brief Finds maximum amplitude for a given channel in time window ["from", "to"] and creates 3d map of waveforms ordered by maxima
/// 
/// Use PrintMaxDist() to plot all channels.
/// Use only for small datasets as it will contain all of the data
/// 
/// @return TGraph2D of all histograms ordered by maximum amplitude
TGraph2D* ReadRun::MaxDist(int channel_index, float from, float to) {
	// find maximum amplitude for a given channel in time window [from, to] and return 3d histogram with the number of bins nbinsy,z

	TString name(Form("maxdist_ch%02d", active_channels[channel_index]));
	TGraph2D* g3d = new TGraph2D((1024 + 2) * nevents);
	g3d->SetTitle("waveforms; t [ns]; max. amplitude [mv]; amplitude [mV]");

	for (int j = 0; j < nevents; j++) {
		if (!skip_event[j]) {
			auto his = (TH1F*)((TH1F*)rundata->At(j * nchannels + channel_index))->Clone();
			if (from >= 0 && to > 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
			double max = his->GetMaximum();
			for (int i = 0; i < 1024; i++) g3d->SetPoint(j * 1024 + i, his->GetXaxis()->GetBinCenter(i), max, his->GetBinContent(i));
		}
	}
	root_out->WriteObject(g3d, name.Data());
	return g3d;
}

/// @brief Finds maximum amplitude for a given channel in time window ["from", "to"] and creates 3d map of waveforms ordered by maxima
/// 
/// Prints MaxDist() for all channels.
/// Use only for small datasets as it will contain all of the data
/// 
/// @param from From
/// @param to To 
void ReadRun::PrintMaxDist(float from, float to) {

	TCanvas* max_dist_c = new TCanvas("wf grouped by maximum", "wf grouped by maximum", 1600, 1000);
	SplitCanvas(max_dist_c);

	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;
			auto g3d = MaxDist(i, from, to);
			max_dist_c->cd(current_canvas);
			g3d->Draw();
		}
	}
	max_dist_c->Update();
	root_out->WriteObject(max_dist_c, "MaxDist");
}

/// @brief Plot results of GetTimingCFD()
/// 
/// See Print_GetTimingCFD() for parameters.
/// 
/// @return Timing histogram for one channel
TH1F* ReadRun::His_GetTimingCFD(int channel_index, float rangestart, float rangeend) {

	int nbins = static_cast<int>((rangeend - rangestart) / SP);

	TString name(Form("GetTimingCFD_ch%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) if (!skip_event[j]) h1->Fill(timing_results[j * nchannels + channel_index][1]);
	return h1;
}

/// @brief Plot results of GetTimingCFD()
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1 fits a gaussian. \n
/// Else do not fit. \n 
/// Fit results per channel are stored in ReadRun::timing_fit_results
void ReadRun::Print_GetTimingCFD(float rangestart, float rangeend, int do_fit) {

	// call GetTimingCFD() in case it was not initialized
	if (timing_results.size() == 0) GetTimingCFD();

	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	TCanvas* timing_cfd_c = new TCanvas("timing of cfd", "timing of cfd", 600, 400);
	SplitCanvas(timing_cfd_c);
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;
			timing_cfd_c->cd(current_canvas);

			TH1F* his;
			his = His_GetTimingCFD(i, rangestart, rangeend);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("time [ns]");
			his->Draw();

			if (do_fit == 1) {
				TFitResultPtr fresults = his->Fit("gaus", "LS", "same");
				timing_fit_results.push_back(fresults);
			}

			TString name_save(Form("Timing_cfd_channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name_save.Data());
		}
	}

	timing_cfd_c->Update();
	root_out->WriteObject(timing_cfd_c, "TimingCFD");
}

/// @brief Plot timining difference between two channels
/// 
/// See Print_GetTimingCFD_diff() for parameters.
/// 
/// @return Histogram with event-wise timing differences between two channels
TH1F* ReadRun::His_GetTimingCFD_diff(int channel1, int channel2, float rangestart, float rangeend) {

	int nbins = static_cast<int>((rangeend - rangestart) / SP);

	auto chin1 = find(active_channels.begin(), active_channels.end(), channel1);
	int channel_index1;
	if (chin1 != active_channels.end())	channel_index1 = chin1 - active_channels.begin();
	else cout << "\n\n ERROR: channel1 does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";

	auto chin2 = find(active_channels.begin(), active_channels.end(), channel2);
	int channel_index2;
	if (chin2 != active_channels.end())	channel_index2 = chin2 - active_channels.begin();
	else cout << "\n\n ERROR: channel2 does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";

	TString name(Form("GetTimingCFD_diff_ch%02d_ch%02d", channel1, channel2));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) if (!skip_event[j]) h1->Fill(timing_results[j * nchannels + channel_index2][1] - timing_results[j * nchannels + channel_index1][1]);
	return h1;
}

/// @brief Plot timining difference between two channels
/// 
/// Plots the difference between the peak times between two channels for each event. It does \f$\Delta t = t_{second} - t_{first}\f$. 
/// 
/// @param channel1 First channel number (wavecatcher channel number). 
/// @param channel2 Second channel number to compare. 
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1 fits a gaussian. \n
/// Else do not fit. \n 
void ReadRun::Print_GetTimingCFD_diff(int channel1, int channel2, float rangestart, float rangeend, int do_fit) {

	// call GetTimingCFD() in case it was not initialized
	if (timing_results.size() == 0) GetTimingCFD();

	gStyle->SetOptStat(1111); 
	gStyle->SetOptFit(111);

	TCanvas* timing_cfd_d_c = new TCanvas("timing of cfd diff", "timing of cfd diff", 600, 400);

	TH1F* his;
	his = His_GetTimingCFD_diff(channel1, channel2, rangestart, rangeend);
	his->GetYaxis()->SetTitle("#Entries");
	his->GetXaxis()->SetTitle("time [ns]");
	his->Draw();

	if (do_fit == 1) {
		TFitResultPtr fresults = his->Fit("gaus", "LS", "same");
		timing_fit_results.push_back(fresults);
	}

	TString name_save(Form("Timing_cfd_diff_channel_%02d_channel_%02d", channel1, channel2));
	root_out->WriteObject(his, name_save.Data());

	timing_cfd_d_c->Update();
	root_out->WriteObject(timing_cfd_d_c, "TimingCFD_diff");
}

// helper functions

/// @brief Helper. Creates a list of .bin data files in data folder to be read in
/// @param dirname Directory
/// @param ext File extension
/// @return Stringstream of line separated file names
stringstream ReadRun::list_files(const char* dirname, const char* ext) {

	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	if (files) {
		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				ss << fname.Data() << "\n";
				//cout << fname.Data() << "\n";
			}
		}
		TIter next2(files);
		while ((file = (TSystemFile*)next2())) {
			fname = file->GetName();
			if (!file->IsDirectory() && !fname.EndsWith(ext) && fname.Contains(ext)) {
				ss << fname.Data() << "\n";
				//cout << fname.Data() << "\n";
			}
		}
	}
	return ss;
}

/// @brief Helper that returns the waveform histogram for a certain channel number and a certain event number
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @param color Choose color of histogram
/// @return Waveform histogram 
TH1F* ReadRun::Getwf(int channelnr, int eventnr, int color) {
	TH1F* his;
	his = (TH1F*)rundata->At(eventnr * nchannels + channelnr);
	his->SetLineColor(color);
	his->SetMarkerColor(color);
	return his;
}

/// @brief Get array of x axis (time) for standard wavecatcher settings
/// @param shift Offset
/// @return Time array
double* ReadRun::getx(double shift) {
	double* xvals = new double[1024];
	for (int i = 0; i < 1024; i++) {
		xvals[i] = static_cast<double>(SP) * static_cast<double>(i) + shift;
	}
	return xvals;
}

/// @brief Get array of y values for a certain waveform
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @return Y values of waveform
double* ReadRun::gety(int channelnr, int eventnr) {
	TH1F* his = Getwf(channelnr, eventnr);
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Get array of y values for a certain waveform
/// @param his Waveform histogram
/// @return Y values of waveform
double* ReadRun::gety(TH1F* his) {
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Get truncated array of y values for a certain waveform
/// @param his Waveform histogram
/// @param start_at Truncate from index
/// @param end_at Truncate to index
/// @return Truncated Y values of waveform
double* ReadRun::gety(TH1F* his, int start_at, int end_at) {
	if (start_at < 0 || start_at >= his->GetNbinsX() || end_at >= his->GetNbinsX() || end_at - start_at < 1) {
		cout << "\nError: ReadRun::gety out of range" << endl;
		return 0;
	}
	const int n_bins_new = end_at - start_at;
	double* yvals = new double[n_bins_new];
	for (int i = start_at; i < end_at; i++) {
		yvals[i - start_at] = his->GetBinContent(i);
	}
	return yvals;
}

/// @brief Translate a random number into a useful root color https://root.cern.ch/doc/master/classTColor.html
/// @param i Index of your plotting loop that is to be translated into a useful ROOT color index
/// @return ROOT color index
int ReadRun::rcolor(int i) {
	int nclrs = 16;
	int rclrs[nclrs] = { 1, 2, 3, 4, 5, 6, 7, 13, 28, 30, 34, 38, 40, 31, 46, 49 };
	return rclrs[i - static_cast<int>(floor(i / nclrs)) * nclrs];
}

/// @brief Returns index of a certain event number (if data files are read in parallel threads)
/// @param eventnr Event number as stored in the data.
/// @return Corresponding event number in the internal data structure.
int ReadRun::GetEventIndex(int eventnr) {
	if (eventnr <= 0) eventnr = 1; // first event is 1
	if (eventnr > nevents) eventnr = nevents;
	return distance(eventnr_storage.begin(), find(eventnr_storage.begin(), eventnr_storage.end(), eventnr));
}

/// @brief Helper to split canvas accordind to the number of channels to be plotted
/// @param c Canvas to be splitted
void ReadRun::SplitCanvas(TCanvas*& c) {
	if (plot_active_channels.empty()) c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else if (plot_active_channels.size() > 1) c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
}

/// @brief Simple linear interpolation for x
/// @param ym Y value for evaluation
/// @param x1 X1 
/// @param x2 X2
/// @param y1 Y1
/// @param y2 Y2
/// @return x value at "ym"
float ReadRun::LinearInterpolation(float ym, float x1, float x2, float y1, float y2) {
	return x1 + (ym - y1) * (x2 - x1) / (y2 - y1);
}

/// @brief Helper to perform convolution of two 1D arrays
/// 
/// Used for smoothing etc.
/// 
/// @param[in,out] result Array containing convolution result
/// @param first First array for convolution
/// @param second Second array for convolution 
/// @param size1 Size of first array (ideally size1 should be equal to size2)
/// @param size2 Size of second array
void ReadRun::Convolute(double*& result, double* first, double* second, int size1, int size2) {

	// uncomment to use sum instead of FFT
	// faster if size1<size2
	//for (int i = 0; i < size2; i++) {
	//	result[i] = 0.;
	//	for (int j = 0; j < TMath::Min(size1, i); j++) {
	//		result[i] += first[j] * second[i - j];
	//	}
	//}

	double* refirst = new double[size1];
	double* imfirst = new double[size1];
	double* resecond = new double[size1];
	double* imsecond = new double[size1];
	double* reres = new double[size1];
	double* imres = new double[size1];

	TVirtualFFT* fftfirst = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftfirst->SetPoints(first);
	fftfirst->Transform();
	fftfirst->GetPointsComplex(refirst, imfirst);
	delete fftfirst;

	TVirtualFFT* fftsecond = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftsecond->SetPoints(second);
	fftsecond->Transform();
	fftsecond->GetPointsComplex(resecond, imsecond);
	delete fftsecond;

	TComplex cofirst;
	TComplex cosecond;
	TComplex cores;

	for (int i = 0; i < size1; i++) {
		cofirst(refirst[i], imfirst[i]);
		cosecond(resecond[i], imsecond[i]);

		cores = cofirst * cosecond / static_cast<double>(size1);

		reres[i] = cores.Re();
		imres[i] = cores.Im();
	}

	//cout << "performing IFFT ... ";
	TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &size1, "C2R ES");
	fft_back->SetPointsComplex(reres, imres);
	fft_back->Transform();
	fft_back->GetPoints(result);
	delete fft_back;
	delete[] imres; delete[] reres; delete[] refirst; delete[] imfirst; delete[] resecond; delete[] imsecond;
}

/// @brief Apply smoothing array of double with length nbins
/// 
/// Very inefficient, try not to use.
/// 
/// @param[in,out] ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Number of bins for running average/gauss sigma in ns for convolution.
/// @param doconv If false use running average (default). \n 
/// If true use gaussian smoothing (slower).
void ReadRun::SmoothArray(double*& ar, int nbins, double sigma, bool doconv) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	if (doconv) {
		// convolution with gauss with sigma (experimental, not yet tested)
		double* gauss = new double[nbins];

		double sum = 0.;

		for (int i = 0; i < nbins; i++) {
			gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i) * SP - 5 * sigma/*-> offset?*/), 2.) / (2. * sigma * sigma)) / (sigma * 2.506628);
			sum += gauss[i];
		}

		for (int i = 0; i < nbins; i++) {
			gauss[i] /= sum;
		}

		Convolute(ar, artmp, gauss, nbins, nbins);
		delete[] gauss;
	}
	else {
		// calculate running average from -sigma until +sigma
		for (int i = 0; i < nbins; i++) {
			double mean1 = 0.;
			int nmn = 0;
			for (int k = -1 * static_cast<int>(floor(sigma)); k <= static_cast<int>(ceil(sigma)); k++) {
				if (i + k >= 0 && i + k < nbins) {
					mean1 += artmp[i + k];
					nmn++;
				}
			}
			if (nmn != 0.) {
				ar[i] = mean1 / static_cast<double>(nmn);
			}
		}
	}
	delete[] artmp;
}

/// @brief Print fourier transform of waveform
/// @param eventnr Event number
/// @param xmin Range min
/// @param xmax Range max
/// @param multiplier Multiplier for resolution of plot
void ReadRun::PrintFFTWF(int eventnr, float xmin, float xmax, int multiplier) {
	// plot waveforms of all channels for a given event number eventnr and add the determined integration windwos to the plot
	TString name(Form("fft_waveforms_event__%04d", eventnr));
	TCanvas* fftc = new TCanvas(name.Data(), name.Data(), 1600, 1000);
	SplitCanvas(fftc);

	TString imname(Form("fft_im_waveforms_event__%04d", eventnr));
	TCanvas* imfftc = new TCanvas(imname.Data(), imname.Data(), 1600, 1000);
	SplitCanvas(imfftc);

	int size = 1024 * multiplier;
	double* xvals = new double[size];
	for (int i = 0; i < size; i++) {
		xvals[i] = static_cast<double>(i) / (SP * static_cast<double>(size));
	}

	double* refft = new double[size];
	double* imfft = new double[size];

	double* yvals = new double[size];

	for (int i = 0; i < nchannels; i++) {
		TH1F* his = Getwf(i, GetEventIndex(eventnr));

		for (int j = 0; j < size; j++) {
			if (j < 1024) yvals[j] = his->GetBinContent(j + 1);
			else yvals[j] = 0.;
		}

		TVirtualFFT* ffthis = TVirtualFFT::FFT(1, &size, "R2C ES");
		ffthis->SetPoints(yvals);
		ffthis->Transform();
		ffthis->GetPointsComplex(refft, imfft);

		TGraph* re = new TGraph(size, xvals, refft);
		TGraph* im = new TGraph(size, xvals, imfft);

		// draw to canvas
		fftc->cd(i + 1);
		stringstream renamess; renamess << "Channel " << active_channels[i] << ", event " << eventnr << ", Re(FFT(data))";
		re->SetTitle(renamess.str().c_str());
		re->Draw("AL");
		gPad->Modified();
		if (xmin != 0. || xmax != 0.) re->GetXaxis()->SetLimits(xmin, xmax);
		//re->GetYaxis()->SetRangeUser(-1*size, size);

		imfftc->cd(i + 1);
		stringstream imnamess; imnamess << "Channel " << active_channels[i] << ", event " << eventnr << ", Im(FFT(data))";
		im->SetTitle(imnamess.str().c_str());
		im->Draw("AL");
		gPad->Modified();
		if (xmin != 0. || xmax != 0.) im->GetXaxis()->SetLimits(xmin, xmax);
		//im->GetYaxis()->SetRangeUser(-1*size, size);

		delete ffthis;
	}
	fftc->Update();
	imfftc->Update();

	root_out->WriteObject(fftc, name.Data());
	root_out->WriteObject(imfftc, imname.Data());

	delete[] yvals;
	delete[] refft;
	delete[] imfft;
	delete[] xvals;
}