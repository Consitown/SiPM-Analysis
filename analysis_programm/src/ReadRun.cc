/// \mainpage Introduction
/// This page serves as documentation of the waveform analysis framework ```wavecatcher-analysis``` for WaveCatcher setups 
/// in the Experimental Elementary Particle Physics Group at the Institute of Physics at Humboldt University of Berlin. \n \n \n
/// 
/// You can find the documentation of the functions and variables at <a href="classReadRun.html"> ReadRun Class Reference</a>
/// 
/// Installation instructions can be found at <https://cscharf-hub.github.io/wavecatcher-analysis/>
/// 
/// The source code is available at <https://github.com/cscharf-hub/wavecatcher-analysis>
/// 
/// Development and maintenance: \n
/// Christian Scharf \n 
/// Contributors: \n 
/// Doramas Jimeno Sanchez \n
/// Alessia Brignoli \n
/// Ben Skodda \n 
/// Christophe Mullesch \n
/// Constantin Eckardt \n
/// Alexander Vagts \n

#include "ReadRun.h"

/// @brief Constructor of the class
ReadRun::ReadRun(int no_of_bin_files_to_read) {

	cout << "\ninitializing ..." << endl;
	if (no_of_bin_files_to_read > 0) cout << "will read <=" << no_of_bin_files_to_read << " bin files" << endl;
	ROOT::EnableImplicitMT();
	TH1::AddDirectory(kFALSE);

	nwf = 0;
	PrintChargeSpectrum_cnt = 0;
	PlotChannelAverages_cnt = 0;
	NoOfBinFilesToRead = no_of_bin_files_to_read;

	root_out = new TFile();	// init results file
}

/// @brief Routine to read files created by the wavecatcher.
/// 
/// Reads the data and can already do simple data manipulation during reading: \n
///   - Change polarity of range of channels with parameters (explained below). \n 
///   - Shift all waveforms to a constant fraction such they all start at the same time (-> see ReadRun::Shift_WFs_in_file_loop ) \n 
///   - Simple baseline correction by calling CorrectBaseline() before calling ReadFile(). \n 
/// Stores the data as ROOT histograms in ```TClonesArray rundata```. \n \n 
/// 
/// Reader modified from WaveCatcher binary -> root converter by manu chauveau@cenbg.in2p3.fr \n 
/// See https://owncloud.lal.in2p3.fr/public.php?service=files&t=56e4a2c53a991cb08f73d03f1ce58ba2
/// 
/// @param path Path to the data. All files in this folder containing ```.bin``` in the file name will be read in.
/// @param change_polarity Set ```true``` to change polarity (sign) of certain channels (see ```change_sign_from_to_ch_num``` below).
/// @param change_sign_from_to_ch_num All channels \f$ \geq \f$ ```change_sign_from_to_ch_num``` will be inverted if ```change_polarity``` is ```true```. \n 
/// If negative number all channels \f$ \leq \f$ ```abs(change_sign_from_to_ch_num)``` will be inverted if ```change_polarity``` is ```true```.
/// @param out_file_name Name of the ```.root``` file which stores the results, e. g. ```results.root```.
/// @param debug Set ```true``` to increase the verbosity.
/// @todo Remove ```change_polarity``` and make back-compatible. 
void ReadRun::ReadFile(string path, bool change_polarity, int change_sign_from_to_ch_num, string out_file_name, bool debug) {

	// save output path
	data_path = path;

	printf("+++ saving analysis results in '%s' ...\n\n", out_file_name.c_str());
	root_out = TFile::Open(out_file_name.c_str(), "recreate");

	/// Wavecatcher hardware max. number of channels (reduce if not using the 64 channel crate)
	const int nChannelsWC = 64;

	rundata = new TClonesArray("TH1F", maxNWF); //raw data will be stored here as TH1F
	rundata->BypassStreamer(kFALSE);			//Doramas: I don't know why is it used, but it's better to use when working with TClonesArray
	TClonesArray& testrundata = *rundata;

	// verbosity
	bool debug_header = debug;
	bool debug_data = debug;

	unsigned short output_channel;
	unsigned int output_event;
	//unsigned long long int output_tdc;
	unsigned short output_nbchannels;
	unsigned short read_channels = 0;

	amplValuessum = new double* [nChannelsWC]; //sum of all wf for each channel
	for (int i = 0; i < nChannelsWC; i++) {//init
		amplValuessum[i] = new double[binNumber];
		for (int k = 0; k < binNumber; k++) amplValuessum[i][k] = 0.;
	}

	maxSumBin = new int[nChannelsWC];

	//Start reading the raw data from .bin files.
	stringstream inFileList;
	inFileList << ListFiles(path.c_str(), ".bin"); //all *.bin* files in folder path
	string fileName;
	int file_counter = 0;
	int wfcounter = 0;
	int event_counter = 0;

	while (inFileList >> fileName) {
		// file loop
		if (NoOfBinFilesToRead > 0 && file_counter >= NoOfBinFilesToRead) break;

		fileName = path + fileName;
		ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);

		bool has_measurement = false;

		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}

		if (file_counter < 10 || file_counter % 10 == 0 || debug) printf("+++ reading '%s' ...\n", fileName.c_str());

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
			if (debug_data) printf("%03d has %d channels\n", an_event.EventNumber, an_event.nchannelstored);

			output_event = an_event.EventNumber;
			//output_tdc = an_event.TDCsamIndex;
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
					memcpy(a_channel_data.waveform, a_channel_data_without_measurement.waveform, binNumber * sizeof(short));
				}


				output_channel = a_channel_data.channel;
				if (debug_data) printf("- reading channel %d\n", output_channel);

				//---------------------------------------------------------------------------------------------------------------
				if (ch >= start_at_ch && ch <= end_at_ch) {

					if (event_counter == 0) active_channels.push_back(static_cast<int>(output_channel));

					TString name(Form("channel%02d_event%05d", output_channel, an_event.EventNumber));
					TString title(Form("Channel %d, event %d;time [ns];signal [mV]", output_channel, an_event.EventNumber));
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

						// do a mini-baseline correction (needed in case a voltage offset is set in the wavecatcher)
						short bsln = (a_channel_data.waveform[0] + a_channel_data.waveform[1] + a_channel_data.waveform[3]) / 3;
						for (int lll = 0; lll < binNumber; lll++) a_channel_data.waveform[lll] -= bsln;

						float global_max = TMath::MaxElement(1024, a_channel_data.waveform);
						//nmax = TMath::LocMax(binNumber, a_channel_data.waveform);
						if (global_max < TMath::Abs(TMath::MinElement(binNumber, a_channel_data.waveform))) {
							global_max = TMath::Abs(TMath::MinElement(binNumber, a_channel_data.waveform));
							//nmax = TMath::LocMin(binNumber, a_channel_data.waveform);
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
						if (change_polarity && ((output_channel >= change_sign_from_to_ch_num) || (change_sign_from_to_ch_num < 0 && output_channel <= abs(change_sign_from_to_ch_num)))) {
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

					wfcounter++;
				}//--------------------------------------------------------------------------------------------------------------

			} // for ch

			skip_event.push_back(false);
			if (!discard_original_eventnr) eventnr_storage.push_back(output_event); // Stores the current WaveCatcher event number
			else eventnr_storage.push_back(event_counter);
			event_counter++;
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
}

/// @brief Destructor
ReadRun::~ReadRun() {
	//rundata->Clear("C");
	//if (maxSumBin) delete[] maxSumBin;
	//plot_active_channels.clear();
	if (root_out->IsOpen()) root_out->Close();
	cout << "\ndeleting nothing currently..." << endl;
}

/// @brief Plot sums of all raw waveforms for each channel
/// 
/// To plot the average waveforms after baseline correction etc. use PlotChannelAverages().
/// 
/// \image html PlotChannelSums.png "This function plots the sum of all waveforms for each channel, without any corrections. Channel 9 was measured with an offset as visible here. Code in example." width=75%
/// 
/// @param smooth If true it will apply smoothing to plots. \n 
/// Do not use without very good reason as it biases the results.
/// @param normalize If true will normalize the maximum to 1.
/// @param shift Shift histogram by "shift" ns
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution.
/// @param smooth_method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks! Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::PlotChannelSums(bool smooth, bool normalize, double shift, double sigma, int smooth_method) {

	double* xv = getx(shift);
	TMultiGraph* mgsums = new TMultiGraph();
	mgsums->SetTitle("channel sums; t [ns]; amplitude [mV]");
	if (normalize) mgsums->SetTitle("channel sums; t [ns]; amplitude [arb.]");

	double max = 0., min = 0.;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			double* yv = amplValuessum[i];
			if (smooth) SmoothArray(yv, binNumber, sigma, smooth_method);

			TGraph* gr = new TGraph(binNumber, xv, yv);
			delete[] yv;

			double tmp_min = TMath::MinElement(gr->GetN(), gr->GetY());
			if (tmp_min < min) min = tmp_min;
			double tmp_max = TMath::MaxElement(gr->GetN(), gr->GetY());
			if (tmp_max > max) max = tmp_max;
			if (normalize) {
				for (int j = 0; j < gr->GetN(); j++) gr->SetPointY(j, gr->GetPointY(j) / tmp_max);
			}

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

	TCanvas* sumc = new TCanvas("Sums", "", 600, 400);
	mgsums->Draw("AL");
	if (normalize) mgsums->GetYaxis()->SetRangeUser(-0.2, 1);
	else mgsums->GetYaxis()->SetRangeUser(min, max);
	sumc->BuildLegend(0.85, 0.70, .99, .95);
	root_out->WriteObject(mgsums, "channelsums");
	root_out->WriteObject(sumc, "channelsums_c");
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Plot averages only of the good, corrected waveforms for each channel
/// 
/// Similar to PlotChannelSums(), but will average all non-skipped waveforms. \n
/// Can be used to inspect average waveforms after baseline correction etc. has been applied. 
/// To do so, call function after calling correction and event filter functions.
/// 
/// \image html PlotChannelAverages.png "This function plots the sum of all non-skipped waveforms for each channel, with corrections. Compare with the result for PlotChannelSums(). Code in example." width=75%
/// 
/// @param normalize If true will normalize the maximum to 1.
void ReadRun::PlotChannelAverages(bool normalize) {
	
	double* xv = getx();
	TMultiGraph* mgav = new TMultiGraph();
	mgav->SetTitle("channel averages; t [ns]; amplitude [mV]");
	if (normalize) mgav->SetTitle("channel averages; t[ns]; amplitude[arb.]");

	double max = 0., min = 0.;
		
	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			
			double* yv = new double[binNumber];
			for (int k = 0; k < binNumber; k++) yv[k] = 0.;

			for (int j = 0; j < nevents; j++) {
				if (!skip_event[j]) {
					double* y_tmp = gety(i, j);
					for (int k = 0; k < binNumber; k++) yv[k] += y_tmp[k];
					delete[] y_tmp;
				}
			}

			double norm = static_cast<double>(Nevents_good());
			for (int k = 0; k < binNumber; k++) yv[k] /= norm;

			TGraph* gr = new TGraph(binNumber, xv, yv);
			delete[] yv;

			double tmp_min = TMath::MinElement(gr->GetN(), gr->GetY());
			if (tmp_min < min) min = tmp_min;
			double tmp_max = TMath::MaxElement(gr->GetN(), gr->GetY());
			if (tmp_max > max) max = tmp_max;
			if (normalize) {
				for (int j = 0; j < gr->GetN(); j++) gr->SetPointY(j, gr->GetPointY(j) / tmp_max);
			}

			TString name(Form("channel_%02d", active_channels[i]));
			TString title(Form("Channel %d", active_channels[i]));
			gr->SetName(name.Data());
			gr->SetTitle(title.Data());
			gr->SetLineColor(rcolor(i));
			gr->SetMarkerColor(rcolor(i));
			mgav->Add(gr);
		}
	}
	delete[] xv;

	string cname("Averages_" + to_string(PlotChannelAverages_cnt++));
	TCanvas* avc = new TCanvas(cname.c_str(), cname.c_str(), 600, 400);
	mgav->Draw("AL");
	if (normalize) mgav->GetYaxis()->SetRangeUser(-0.2, 1);
	else mgav->GetYaxis()->SetRangeUser(min, max);
	avc->BuildLegend(0.85, 0.70, .99, .95);
	root_out->WriteObject(mgav, ("channelaverages" + to_string(PlotChannelAverages_cnt)).c_str());
	root_out->WriteObject(avc, ("channelaverages_c" + to_string(PlotChannelAverages_cnt)).c_str());
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Smoothing all waveforms which are not skipped (for testing, careful when using for analysis!)
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution.
/// @param method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
void ReadRun::SmoothAll(double sigma, int method) {
	// just for testing, not very efficient
	cout << "\nsmoothing wfs";
	for (int j = 0; j < nwf; j++) {
		if (!skip_event[GetCurrentEvent(j)]) {
			TH1F* his = Getwf(j);
			double* yvals = gety(his);
			SmoothArray(yvals, binNumber, sigma, method, SP);
			for (int i = 1; i <= his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i - 1]);
			delete[] yvals;
		}
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

/// @brief Filter all waveforms
/// 
/// Experimental. See FilterArray().
/// 
/// @param sigma1 First.
/// @param sigma2 Second.
/// @param factor Factor for negative part (0 < factor < 1).
void ReadRun::FilterAll(double sigma1, double sigma2, double factor) {
	cout << "\nfiltering wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = Getwf(j);
		double* yvals = gety(his);
		FilterArray(yvals, binNumber, sigma1, sigma2, factor);
		for (int i = 1; i <= his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i - 1]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

/// @brief This function shifts all waveforms to the average signal starting times for each channel.
/// 
/// The signal starting times are determined with constant fraction discrimination. 
/// **Before** calling this function, please call GetTimingCFD() with suitable parameters for your data. \n
/// Additionally, PrintChargeSpectrumWF() should be called **before** calling this function since the timing reference (blue line) won't be shifted.
/// 
void ReadRun::ShiftAllToAverageCF() {
	cout << "\nshifting all WFs to the average CF time for each channel.\n";
	
	//call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();
	
	double* timing_mean = new double[nchannels];
	for (int i = 0; i < nchannels; i++) timing_mean[i] = 0.;
	
	for (int j = 0; j < nwf; j++) {
		if (!skip_event[GetCurrentEvent(j)]) timing_mean[GetCurrentChannel(j)] += timing_results[j][0];
	}

	double norm = static_cast<double>(Nevents_good());

	int* timing_mean_n = new int[nchannels];
	for (int i = 0; i < nchannels; i++) {
		timing_mean_n[i] = static_cast<int>(round(timing_mean[i] / norm ));
	}
	delete[] timing_mean;

	for (int j = 0; j < nwf; j++) {
		if (!skip_event[GetCurrentEvent(j)]) {
			TH1F* his = Getwf(j);
			double* yvals = gety(his);
			int shift = static_cast<int>(timing_results[j][0]) - timing_mean_n[GetCurrentChannel(j)];

			for (int i = 0; i < his->GetNbinsX(); i++) {
				int icycle = 0;
				if (i + shift >= his->GetNbinsX()) icycle = -1 * his->GetNbinsX();
				else if (i + shift < 0) icycle = his->GetNbinsX();
				his->SetBinContent(i + 1, yvals[i + shift + icycle]);
			}
			delete[] yvals;
		}
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
	delete[] timing_mean_n;
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
			CorrectBaseline_function(Getwf(j), tCut, tCutEnd, j);

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
		for (int i = 1; i <= his->GetNbinsX(); i++) his->SetBinContent(i, his->GetBinContent(i) - corr);
	}

	if (!Using_BaselineCorrection_in_file_loop) {
		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[nwaveform].push_back(corr);
		baseline_correction_result[nwaveform].push_back(0);
		baseline_correction_result[nwaveform].push_back(tCut);
		baseline_correction_result[nwaveform].push_back(tCutEnd);
	}
}

/// @brief Baseline correction method searching for non-monotonic, rather constant regions
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Determines the region of "nIntegrationWindow" bins where the squared sum plus the square of the sum of the slope of the smoothed waveform reaches its minimum: \n \n
/// \f$\mathbf{min}\left( \sum \left(\Delta y_i \right)^2 + \left(\sum \Delta y_i \right)^2 \right) \f$ \n \n
/// 
/// Here, \f$\sum \left(\Delta y_i \right)^2 \to 0\f$ if the region is constant and \f$\left( \sum \Delta y_i \right)^2 \to 0\f$ if the region is constant or 
/// oscillating around a constant value. The second term penalizes regions with a small monotonic slope (e. g. tails). \n \n 
/// 
/// Slow but most versatile since it searches for the optimal baseline candidate region in a defined range. \n 
/// Will prefer constant sections of the waveform for the estimation of the baseline. \n \n 
/// 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n \n 
/// 
/// @param nIntegrationWindow Number of bins of integration window
/// @param smooth If true will use averaging for more reliable slope. Use with care!
/// @param sigma NNumber of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution. Use with care!
/// @param max_bin_for_baseline Maximum bin for search window.
/// @param start_at Minimum bin for search window.
/// @param search_min Experimental, use with care.
/// @param smooth_method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
/// @param skip_channel Skip a channel
/// @todo Work on "skip_channel" and remove "search_min"
void ReadRun::CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool smooth, double sigma, int max_bin_for_baseline, int start_at, bool search_min, int smooth_method, int skip_channel) {

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

		imax = 0;
		search_before = 0;

		if (j == 0 || j != skip_channel - 1 || j % skip_channel != 0) { //eventnr * nchannels + i
			TH1F* his = Getwf(j);
			double* yvals = gety(his); //find faster way
			if (sigma > 0) SmoothArray(yvals, binNumber, sigma, smooth_method); // smoothing important to suppress variations in slope due to noise so the method is more sensitive to excluding peaks

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
			if (!smooth) {
				corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow + 1);
			}
			else {
				for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
				corr /= static_cast<float>(nIntegrationWindow);
			}

			for (int i = 0; i < binNumber; i++) {
				if (!smooth) his->SetBinContent(i + 1, his->GetBinContent(i + 1) - corr);
				else his->SetBinContent(i + 1, yvals[i] - corr);
			}
			delete[] yvals;
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
}
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Baseline correction using minimum sum (\f$\propto\f$ mean) in range for correction 
/// 
/// Corrects the baseline (DC offset) of all waveforms. \n 
/// Uses \f$\mathbf{min}\left( \sum y_i \right) \f$ over "nIntegrationWindow" in range {"start_at", "max_bin_for_baseline"}. \n
/// Make sure the search range is shortly before the triggered signal is expected to arrive. \n \n 
/// 
/// Helpful for (groups of/irradiated) SiPMs with very high dark count rate where the signal voltage rarely relaxes back to the baseline before the next signal arrives: \n
/// \f$ \Rightarrow DCR \sim 1/t_{signal} \f$ \n \n
/// 
/// Estimates a baseline as the minimum level before the main peak. \n 
/// Stores results for all channels and all events in ReadRun::baseline_correction_result. \n 
/// Results will be visualized for each event in PrintChargeSpectrumWF(). \n 
/// For parameters see CorrectBaselineMinSlopeRMS() 
/// 
void ReadRun::CorrectBaselineMin(int nIntegrationWindow, double sigma, int max_bin_for_baseline, int start_at, int smooth_method, int skip_channel) {
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
			TH1F* his = Getwf(j);
			double* yvals = gety(his); //find faster way
			if (sigma > 0) SmoothArray(yvals, binNumber, sigma, smooth_method); // smoothing

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

			corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow + 1);

			for (int i = 0; i < binNumber; i++) {
				his->SetBinContent(i + 1, his->GetBinContent(i + 1) - corr);
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
/// @example timing_example.cc

/// @brief Determine the timing of the maximum peak with constant fraction discrimination
/// 
/// Determines timing in the time interval ["start_at_t", "end_at_t"] with CFD for fraction of maximum "cf_r". \n 
/// Stores timing information for all channels and all events in ReadRun::timing_results. \n 
/// Per event results will be visualized in PrintChargeSpectrumWF(). \n 
/// Cumulative per channel results can be visualized with Print_GetTimingCFD(). \n 
/// Results can be used for time difference cuts with SkipEventsTimeDiffCut().
/// 
/// @param cf_r Fraction of maximum for CFD. \n
/// If set > 1 it will be used as a fixed threshold in mV, so it will no longer use CFD.
/// @param start_at_t Time in ns to start searching.
/// @param end_at_t Time in ns to end searching.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution.
/// This will bias the results! Do not use (or use very carefully, only for noisy data)! Set to 0 if you do not want to use smoothing.
/// @param find_CF_from_start If true will start search from "start_at_t" to find the first arriving photon (default setting). \n 
/// If false search backwards from the time of the maximum.  
/// @param smooth_method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
/// @param use_spline If false will use linear interpolation between the two bins closest to cf_r. \n
/// If true will use a 5th order spline and bisection method for interpolation. 
/// Performs a bit better in most cases compared to only linear interpolation. 
void ReadRun::GetTimingCFD(float cf_r, float start_at_t, float end_at_t, double sigma, bool find_CF_from_start, int smooth_method, bool use_spline) {

	int start_at = static_cast<int>(floor(start_at_t / SP));
	int end_at = static_cast<int>(ceil(end_at_t / SP));
	int n_range = end_at - start_at;

	if (cf_r > 0 && cf_r <= 1) printf("\nGet timing at CF=%.2f between %.2f ns and %.2f ns (%d waveforms) :: ", cf_r, start_at_t, end_at_t, nwf);
	else printf("\nGet timing at threshold=%.2f between %.2f ns and %.2f ns (%d waveforms) :: ", cf_r, start_at_t, end_at_t, nwf);

	for (int j = 0; j < nwf; j++) {

		TH1F* his = Getwf(j);
		double* yvals = gety(his, start_at, end_at); // get range where to search for CFD for timing

		if (sigma > 0.) SmoothArray(yvals, n_range, sigma, smooth_method); // smoothing to suppress noise, will also change timing so use with care!

		float max = 0.;
		int n_max = 0;
		int i = 0;
		for (i = 0; i < n_range; i++) {
			if (yvals[i] > max) {
				max = yvals[i];
				n_max = i;
			}
		}

		float cf = cf_r;
		if (cf_r > 0 && cf_r <= 1) cf *= max;

		if (!find_CF_from_start) {
			i = n_max;
			while (yvals[i] > cf && i >= 0) i--;
			i++;
		}
		else {
			i = 0;
			while (yvals[i] < cf && i < n_max) i++;
			i--;
		}

		// do interpolation for cf
		float interpol_bin = .0;
		interpol_bin = LinearInterpolation(cf, static_cast<float>(i), static_cast<float>(i + 1), yvals[i], yvals[i + 1]);
		// go to center of bin
		interpol_bin += .5;

		if (use_spline) { // use spline interpolation with tolerance epsilon*bin_size
			double epsilon = 1e-4;
			double x_low = interpol_bin - .5;
			double x_high = interpol_bin + .5;

			double* xvals = new double[n_range];
			for (int k = 0; k < n_range; k++) xvals[k] = static_cast<double>(k) + .5;

			TSpline5* wfspl = 0;
			wfspl = new TSpline5("wf_spline", xvals, yvals, n_range, "b1e1b2e2", 0., 0., 0., 0.);
			delete[] xvals;
			
			// using bisection method: halving search window until cf is less than epsilon bins from spline value
			while (x_high - x_low > epsilon) {
				double x_mid = (x_low + x_high) / 2;
				double f_mid = wfspl->Eval(x_mid);
				if (f_mid == cf) break;

				if (f_mid > cf) x_high = x_mid;
				else x_low = x_mid;
			}
			interpol_bin = (x_low + x_high) / 2;
			delete wfspl;
		}

		timing_results.push_back(vector<float>());
		timing_results[j].push_back(interpol_bin);											// the bin we looked for
		timing_results[j].push_back((interpol_bin + static_cast<float>(start_at)) * SP);	// cfd-time we looked for
		timing_results[j].push_back(max);													// maximum value
		timing_results[j].push_back(n_max);													// bin of maximum
		timing_results[j].push_back(cf);													// constant fraction
		timing_results[j].push_back(static_cast<float>(start_at) * SP);						// starting time
		timing_results[j].push_back(static_cast<float>(end_at) * SP);						// end time
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
	cout << endl;
}
/// @example timing_example.cc

/// @brief Skip events where the time difference between two channels is outside of specified range
/// 
/// Skip events where the time difference between channel "first_channel" and channel "second_channel" 
/// is less than "time_diff_min" or more than "time_diff_max" \n 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param first_channel_abs First channel.
/// @param second_channel_abs Second channel.
/// @param time_diff_min Skip events for \f$ \Delta t<t_{diff,min} \f$
/// @param time_diff_max Skip events for \f$ \Delta t>t_{diff,max} \f$
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsTimeDiffCut(int first_channel_abs, int second_channel_abs, double time_diff_min, double time_diff_max, bool verbose) {

	printf("\n\n Removing events if the event-wise time difference between the main peaks in ch%d and ch%d is <%.2f ns or >%.2f ns\n\n", first_channel_abs, second_channel_abs, time_diff_min, time_diff_max);
	int counter = 0;
	int first_channel = GetChannelIndex(first_channel_abs);
	int second_channel = GetChannelIndex(second_channel_abs);

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	// loop through events, calculate timing difference between channels and compare with cuts
	for (int j = 0; j < nwf; j += nchannels) {
		if (!skip_event[GetCurrentEvent(j)]) {
			float time_diff = timing_results[j + second_channel][1] - timing_results[j + first_channel][1];

			if (j <= static_cast<int>(timing_results.size()) && (time_diff < time_diff_min || time_diff > time_diff_max)) {
				int currevent = eventnr_storage[GetCurrentEvent(j)];
				if (verbose) cout << "\nevent:\t" << currevent << "\tchannels:\t" << first_channel_abs << " & " << second_channel_abs << "\ttime diff:\t" << time_diff;
				skip_event[GetCurrentEvent(j)] = true;
				counter++;
			}
		}
	}
	cout << "\n\n\t" << counter << " events will be cut out of " << nevents << "\n\n";
}
/// @example timing_example.cc

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
	vector<int> counter_abovethr(static_cast<int>(plot_active_channels.size()));	// DORAMAS: It stores a counter of events above threshold for each channel that will be plotted

	cout << "\n\n ------> ";
	if (max) cout << "max";
	else cout << "min";

	if (greater) cout << " > ";
	else cout << " < ";
	cout << threshold << " mV:\n";

	for (int j = 0; j < nwf; j++) {
		auto his = (TH1F*)(Getwf(j))->Clone(); // use Clone() to not change ranges of original histogram

		// set range where to search for amplitudes above threshold
		if (from >= 0 && to > 0) {
			his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
		}

		currchannel = j - nchannels * GetCurrentEvent(j);
		if (find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[currchannel]) != plot_active_channels.end()) {
			if ((max && greater && his->GetMaximum() > threshold) || (max && !greater && his->GetMaximum() < threshold) || (!max && greater && his->GetMinimum() > threshold) || (!max && !greater && his->GetMinimum() < threshold)) {
				currevent = eventnr_storage[GetCurrentEvent(j)];

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
		delete his;
	}

	if (verbose) cout << endl;

	//  Loop to show the fraction of events above threshold for each channel that will be plotted
	for (int i = 0; i < static_cast<int>(plot_active_channels.size()); i++) {
		cout << "\nfraction of events in channel " << plot_active_channels[i] << " above threshold: " << 100. * static_cast<float>(counter_abovethr[i]) / static_cast<float>(nevents) << "%\n";
	}
	//
	cout << "\nfraction of events w/ at least 2 channels above threshold: " << 100. * static_cast<float>(occurences2ch) / static_cast<float>(nevents) << "%\n";
	cout << "\tfor a total of " << nevents << " events\n" << endl;
}

/// @brief Skip events above/below individual thresholds per channel
/// 
/// Needs to be called before the charge spectrum etc functions. \n 
/// Baseline correction should be called before this function.
/// 
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). Negative thresholds mean events below threshold will be cut. If only one value is given this value will be used for all active channels.
/// @param rangestart Range start in ns
/// @param rangeend Range end in ns
/// @param verbose Set true for extra verbosity.
void ReadRun::SkipEventsPerChannel(vector<double> thresholds, double rangestart, double rangeend, bool verbose) { // merge with IntegralFilter()?

	if (thresholds.empty()) cout << "\nError: thresholds is empty";
	while (thresholds.size() <= plot_active_channels.size()) thresholds.push_back(thresholds[0]);

	cout << "\n\n Removing events with individual amplitude threshold per channel!!!\n\n";
	int counter = 0;

	for (int j = 0; j < nwf; j++) {
		if (!skip_event[GetCurrentEvent(j)]) {
			int currchannel = j - nchannels * GetCurrentEvent(j);
			if (currchannel < static_cast<int>(thresholds.size())) {
				auto his = (TH1F*)(Getwf(j))->Clone(); // use Clone() to not change ranges of original histogram
				if (rangestart != 0 && rangeend != 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(rangeend));

				if (thresholds[currchannel] != 0. && !skip_event[GetCurrentEvent(j)] && ((thresholds[currchannel] > 0 && his->GetMaximum() > thresholds[currchannel]) || (thresholds[currchannel] < 0 && his->GetMinimum() < thresholds[currchannel]))) {

					int currevent = eventnr_storage[GetCurrentEvent(j)];
					if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\tthreshold\t" << thresholds[currchannel];
					skip_event[GetCurrentEvent(j)] = true;
					counter++;
				}
				delete his;
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
/// @param thresholds Vector should contain a threshold for each active channel saved in the data, in ascending order (ch0, ch1 ...). Negative thresholds mean events below threshold will be cut. A threshold of 0 means the channel will not be evaluated. If only one value is given this value will be used for all active channels.
/// @param highlow  Vector should contain a bool for each active channel. True means events with integrals above threshold will be cut, false means below threshold. If only one value is given this value will be used for all active channels.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param use_AND_condition If set true it will overrule previously applied cuts and un-skip events that pass integral criterium.
/// @param verbose Set true for extra verbosity.
void ReadRun::IntegralFilter(vector<double> thresholds, vector<bool> highlow, float windowlow, float windowhi, float start, float end, bool use_AND_condition, bool verbose) {

	if (thresholds.empty() || highlow.empty()) cout << "\nError: thresholds or highlow are empty";
	while (thresholds.size() <= plot_active_channels.size()) { 
		thresholds.push_back(thresholds[0]);
		highlow.push_back(highlow[0]);
	}

	cout << "\n\nRemoving events with individual integral threshold per channel :: ";
	int counter = 0;
	float integral = 0;
	int currevent_counter = 0;
	int currchannel = 0;
	int currevent = 0;

	for (int j = 0; j < nwf; j++) {
		currevent_counter = GetCurrentEvent(j);

		if (use_AND_condition || !skip_event[currevent_counter]) {
			currchannel = GetCurrentChannel(j);// -nchannels * GetCurrentEvent(j);

			if (currchannel < static_cast<int>(thresholds.size())) {
				auto his = (TH1F*)(Getwf(j))->Clone(); // use Clone() to not change ranges of original histogram
				integral = GetPeakIntegral(his, windowlow, windowhi, start, end, currchannel);

				if (thresholds[currchannel] != 0 && !skip_event[currevent_counter] && 
					((highlow[currchannel] && integral > thresholds[currchannel]) || (!highlow[currchannel] && integral < thresholds[currchannel]))) {
					currevent = eventnr_storage[currevent_counter];
					if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\tthreshold\t" << thresholds[currchannel] << "\tintegral:\t" << integral;
					skip_event[currevent_counter] = true;
					while (floor((j + 1) / nchannels) == currevent_counter) j++;
					counter++;
				}
				else if (use_AND_condition && thresholds[currchannel] != 0 && skip_event[currevent_counter] && !highlow[currchannel]) {
					currevent = eventnr_storage[currevent_counter];
					skip_event[currevent_counter] = false;
					if (verbose) cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel] << "\thas been flagged good by integral:\t" << integral;
				}
				delete his;
			}
			if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
		}
	}

	cout << "\n\n\t" << counter << " events will be cut out of " << nevents << endl;
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Prints a list of all skipped events into the terminal for diagnostics
void ReadRun::PrintSkippedEvents() {
	int counter = 0;
	for (int j = 0; j < static_cast<int>(skip_event.size()); j++) {
		if (skip_event[j]) {
			int currevent = eventnr_storage[j];
			cout << "\nevent:\t" << currevent;
			counter++;
		}
	}
	cout << "\n\ntotal number of skipped events:\t" << counter << "\tout of:\t" << nevents << endl;
}

/// @brief Sets skip_event flag to false for all events, removing any previous cuts
void ReadRun::UnskipAll() {
	for (int j = 0; j < static_cast<int>(skip_event.size()); j++) skip_event[j] = false;
	cout << "\n\nAll event cuts were removed" << endl;
}

/// @brief Number of good events that are not skipped
int ReadRun::Nevents_good() {
	int nevents_good = 0;
	for (int i = 0; i < nevents; i++) if (!skip_event[i]) nevents_good++;
	return nevents_good;
}

// functions for charge spectrum

/// @brief Determine indices for integration window for peaks
/// 
/// Default usage: Find maximum in range ("start", "end") and return bin numbers for [0] the max, [1] t_max - "windowlow", and [2] t_max + "windowhi" \n 
/// If ("start" < 0 || "end" < 0) doesn't return max and integration window is fixed t(max(sum_spectrum[channel])) +/- "windowhi"/"windowlow" \n 
/// If ("windowlow" == "start" && "windowhi" == "end") doesn't return max and sets fixed integration window from "start" until "end" for all channels.
/// 
/// @param his Histogram to integrate.
/// @param windowlow Integration time left to the maximum of the peak.
/// @param windowhi Integration time right to the maximum of the peak.
/// @param start Range for finding maximum.
/// @param end Range for finding maximum.
/// @param channel Channel index in case the integration should be around the maximum of the sum of all waveforms
/// @return Array { max, \f$ n_{t,start} \f$ , \f$ n_{t,end} \f$ }
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
/// 
/// Calculated over an integer number of bins.
/// 
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
/// See GetIntWindow() for explanation of parameters. \n
/// Will also add CFD timing if GetTimingCFD() was called before. 
/// 
/// image html PrintChargeSpectrumWF.png "Waveforms in all channels for a single event. Code in example."
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param eventnr Event number
/// @param ymin Y axis minimum
/// @param ymax Y axis maximum
/// @param xmin X axis maximum
/// @param xmax X axis maximum
void ReadRun::PrintChargeSpectrumWF(float windowlow, float windowhi, float start, float end, int eventnr, float ymin, float ymax, float xmin, float xmax) {

	gStyle->SetOptStat(0);

	TString name(Form("waveforms_event__%05d", eventnr));
	TCanvas* intwinc = new TCanvas(name.Data(), name.Data(), 600, 400);
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
			// formatting
			gPad->SetTopMargin(.01);
			int last_canvas = nchannels;
			if (!plot_active_channels.empty()) last_canvas = static_cast<int>(plot_active_channels.size());
			if (current_canvas == 1 && last_canvas < 4) gPad->SetLeftMargin(.15);
			if (current_canvas % 4 == 0 || current_canvas == last_canvas) gPad->SetRightMargin(.01);
			his->Draw("HIST");
			his->SetStats(0); 

			low->Draw("same");
			hi->Draw("same");
			zero->Draw("same");

			// draw baseline parameters
			if (static_cast<int>(baseline_correction_result.size()) > event_index * nchannels + i) {
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

			if (static_cast<int>(timing_results.size()) > event_index * nchannels + i) {
				TLine* timing = new TLine(timing_results[event_index * nchannels + i][1], -10, timing_results[event_index * nchannels + i][1], 100);
				timing->SetLineColor(9);
				timing->SetLineWidth(2);
				timing->Draw("same");
			}

			if (ymin != 0. && ymax != 0.) his->GetYaxis()->SetRangeUser(ymin, ymax); // fix y range for better comparison 
			if (xmin != 0. && xmax != 0.) his->GetXaxis()->SetRangeUser(xmin, xmax);
		}
	}
	intwinc->Update();

	root_out->WriteObject(intwinc, name.Data());
}
/// @example timing_example.cc

/// @brief Returns array with the individual "charge"/amplitude for all events of one channel
/// 
/// See SaveChargeLists() and GetIntWindow().
/// 
/// @param channel_index Index of the channel.
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param negative_vals If true will save negative values. If false will set negative values to 0.
float* ReadRun::ChargeList(int channel_index, float windowlow, float windowhi, float start, float end, bool negative_vals) {
	float* charge_list = new float[nevents];
	for (int j = 0; j < nevents; j++) {
		TH1F* his = Getwf(j * nchannels + channel_index);
		charge_list[j] = GetPeakIntegral(his, windowlow, windowhi, start, end, channel_index);
		if (!negative_vals && charge_list[j] < 0.) charge_list[j] = 0.;
	}
	return charge_list;
}

/// @brief Saves TGraphs to root file with the individual "charge"/amplitude for all events and all channels
/// 
/// Event with a skip_event flag will be removed. Call before filtering or call UnskipAll() before to get all events. \n
/// See GetIntWindow().
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns.
/// @param negative_vals If true will save negative values. If false will set negative values to 0.
void ReadRun::SaveChargeLists(float windowlow, float windowhi, float start, float end, bool negative_vals) {
	float* event_list = new float[nevents];
	for (int i = 0; i < nevents; i++) event_list[i] = static_cast<float>(i);

	TMultiGraph* charge_list_mg = new TMultiGraph();
	if (windowlow + windowhi > 0.) charge_list_mg->SetTitle("event-wise integrals; Event number; integral [mV#timesns]");
	else charge_list_mg->SetTitle("event-wise amplitudes; Event number; amplitude [mV]");

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			TString name(Form("charge_list_ch_%02d", active_channels[i]));
			float* charge_list = ChargeList(i, windowlow, windowhi, start, end, negative_vals);
			TGraph* charge_list_graph = new TGraph(nevents, event_list, charge_list);
			charge_list_graph->SetLineWidth(0);
			charge_list_graph->SetMarkerStyle(2);
			charge_list_graph->SetMarkerColor(rcolor(i));
			charge_list_graph->SetTitle(name.Data());

			//remove skipped events
			for (int j = 0; j < nevents; j++) {
				if (skip_event[j]) charge_list_graph->RemovePoint(j);
			}

			charge_list_mg->Add(charge_list_graph);
			root_out->WriteObject(charge_list_graph, name.Data());
			delete[] charge_list;
		}
	}
	root_out->WriteObject(charge_list_mg, "all_charge_lists");
	delete[] event_list;
}

/// @brief Plot correlation of integrals/amplitudes between two channels
/// 
/// See GetIntWindow() and PrintChargeSpectrum() for parameters.
/// 
/// \image html ChargeCorrelation.png "Integrals of the signals of two channels for each event plottet against each other. A clear correlation is visible: On average, if channel 2 records a large signal also channel 3 records a large signal. Around (0,0) the correlation of the pedestals can be seen (no signal in both channels). Code in example." width=75% 
/// 
/// @param windowlow Integrate from "windowlow" ns from max...
/// @param windowhi ...to "windowhi" ns from max.
/// @param start Find max from "start" in ns...
/// @param end ...to "end" in ns. 
/// @param rangestart Plot x & y range start
/// @param rangeend Plot x & y range end
/// @param nbins Number of x & y bins of the histogram
/// @param channel1 First channel number to compare
/// @param channel2 second channel number to compare
/// @param ignore_skipped_events Set true to plot only events which passed filtering, else all events will be plotted
void ReadRun::ChargeCorrelation(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, int channel1, int channel2, bool ignore_skipped_events) {
	gStyle->SetOptStat(1111);
	stringstream name;
	name << "charge_correlation_ch" << channel1 << "_ch" << channel2;
	stringstream title;
	if (windowlow + windowhi > 0.) title << ";integral ch" << channel1 << " in mV#timesns;integral ch" << channel2 << " in mV#timesns;Entries";
	else title << ";amplitude ch" << channel1 << " in mV;amplitude ch" << channel2 << " in mV;Entries";

	auto charge_corr_canvas = new TCanvas(name.str().c_str(), "canvas", 600, 600);
	charge_corr_canvas->SetRightMargin(0.15);

	float* charge1 = ChargeList(GetChannelIndex(channel1), windowlow, windowhi, start, end);
	float* charge2 = ChargeList(GetChannelIndex(channel2), windowlow, windowhi, start, end);

	auto charge_corr = new TH2F(name.str().c_str(), title.str().c_str(), nbins, rangestart, rangeend, nbins, rangestart, rangeend);
	for (int i = 0; i < nevents; i++) {
		if (!ignore_skipped_events || !skip_event[i]) charge_corr->Fill(charge1[i], charge2[i]);
	}
	
	charge_corr->Draw("colz");
	root_out->WriteObject(charge_corr, name.str().c_str());

	charge_corr_canvas->Update();
	charge_corr_canvas->SetGrid();
	// move stat box out of the way (causing problems since May 23?)
	//TPaveStats* stat_box = (TPaveStats*)charge_corr->FindObject("stats"); 
	//stat_box->SetX1NDC(0.6); 
	//stat_box->SetX2NDC(0.85);
	charge_corr->SetStats(0);
	charge_corr_canvas->Modified();
	name << "_c";
	root_out->WriteObject(charge_corr_canvas, name.str().c_str());
	delete[] charge1;
	delete[] charge2;
}
/// @example timing_example.cc

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
			TH1F* his = Getwf(j * nchannels + channel_index);
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
/// \image html PrintChargeSpectrum.png "Simple example of the integrated signals for a single channel. The resulting integrated signals are fitted with a Landau-Gauss convolution (-> energy deposition of a minimum ionizing particle in a thin absorber). Code in example." width=75%
/// 
/// \image html cosmics-fit-example.png "Integrated signals for a SiPM in blue and a fit with the fit function Fitf for SiPMs missing after-pulses and dark counts in red. The spectrum is not well fit by the simplified model. You can try different fit models in the example cosmics-fit.ipynb." width=50%
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
/// @param max_channel_nr_to_fit Fit only channels with index < "max_channel_nr_to_fit". Set to -1 to skip fitting.
/// @param which_fitf Choose fit function: \n 
/// 0 - do not fit \n
/// 1 - landau gauss convolution for large number of photons \n 
/// 2 - if pedestal is biased because of peak finder algorithm \n 
/// 3 - SiPM fit function with exponential delayed after pulsing \n 
/// 4 - ideal PMT fit function \n 
/// 5 - PMT fit function \n 
/// 6 - PMT fit function with biased pedestal \n 
/// 7 - default SiPM fit function + dark count spectrum (for lots of false triggers) \n 
/// else - default SiPM fit function \n 
void ReadRun::PrintChargeSpectrum(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit, int which_fitf) {
	// print ReadRun::ChargeSpectrum for all channels optimized for SiPM signals

	gStyle->SetOptStat("ne");
	gStyle->SetOptFit(1111);

	if (fitrangestart == 0.) fitrangestart = rangestart;
	if (fitrangeend == 0.) fitrangeend = rangeend;

	string ctitle("\"charge\" spectra" + to_string(PrintChargeSpectrum_cnt++));
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
	SplitCanvas(chargec);

	cout << "\n\nThere is data recorded in " << active_channels.size() << " channels \n\n\n";
	int current_canvas = 0;

	float default_rangestart = -2000;
	float default_rangeend = 30000;
	int default_nbins = static_cast<int>((default_rangeend - default_rangestart) * nbins / (rangeend - rangestart));

	TH1F* his;
	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;
			
			his = ChargeSpectrum(i, windowlow, windowhi, start, end, default_rangestart, default_rangeend, default_nbins);
			his->GetYaxis()->SetTitle("#Entries");
			if (windowlow + windowhi > 0.) his->GetXaxis()->SetTitle("integral in mV#timesns");
			else his->GetXaxis()->SetTitle("amplitude in mV");

			//store the mean integral of each channel --> used for correction factors of phi_ew analysis
			mean_integral.push_back(his->GetMean());

			chargec->cd(current_canvas);

			//fitting
			if (i < max_channel_nr_to_fit) {
				if (which_fitf == 0) {}
				else if (which_fitf == 1) { // landau gauss convolution for large number of photons
					Fitf_langaus fitf;
					int n_par = 4;
					TF1* f = new TF1("fitf_langaus", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "Width");				f->SetParameter(0, 35);
					f->SetParName(1, "MPV");				f->SetParameter(1, 1000);
					f->SetParName(2, "Area");			    f->SetParameter(2, 10000);
					f->SetParName(3, "#sigma_{Gauss}");		f->SetParameter(3, 100);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(4, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 2) { // if pedestal is biased because of peak finder algorithm
					Fitf_biased fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf_biased", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4); //f->SetParLimits(4, 1.e-9, 1.e3);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 10.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
					f->SetParName(7, "norm_{0}");			f->SetParameter(7, 0.7);
					f->SetParName(8, "x_{0}");				f->SetParameter(8, 5.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

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
				else if (which_fitf == 3) { // SiPM fit function with exponential delayed after pulsing
					Fitf_full fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
					f->SetParName(7, "#alpha");				f->SetParameter(7, .1); //f->FixParameter(7, .2);
					f->SetParName(8, "#beta");				f->SetParameter(8, 80.); //f->FixParameter(8, 80);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 4) { // ideal PMT fit function
					Fitf_PMT_ideal fitf;
					int n_par = 4;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 1.);
					f->SetParName(2, "#sigma");				f->SetParameter(2, 5.);
					f->SetParName(3, "gain");				f->SetParameter(3, 10.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 5) { // PMT fit function
					Fitf_PMT fitf;
					int n_par = 8;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-99, 4.e-1); //probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-99, 5.e-2); //coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 6) { // PMT fit function with biased pedestal
					Fitf_PMT_pedestal fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "Q_{1}");				f->SetParameter(7, 10.);	f->SetParLimits(7, 1.e-9, 1.e9);
					f->SetParName(8, "A_{0}");				f->SetParameter(8, 1.);		f->SetParLimits(8, 1.e-9, 1.e1);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else if (which_fitf == 7) { // default SiPM fit function + dark count spectrum (for lots of false triggers)
					Fitf_plus_DC fitf;
					int n_par = 9;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "A");					f->SetParameter(0, his->Integral());
					f->SetParName(1, "w");					f->SetParameter(1, .05);	f->SetParLimits(1, 1.e-9, 4.e-1);	//probability for type II BG
					f->SetParName(2, "#alpha");				f->SetParameter(2, .05);	f->SetParLimits(2, 1.e-9, 5.e-2);	//coefficient of exponential decrease of typ II BG
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 5.);		f->SetParLimits(3, 1.e-9, 1.e3);
					f->SetParName(4, "Q_{0}");				f->SetParameter(4, 0.);		f->SetParLimits(4, -1.e2, 1.e2);
					f->SetParName(5, "#mu");				f->SetParameter(5, 1.);		f->SetParLimits(5, 1.e-9, 1.e2);
					f->SetParName(6, "#sigma_{1}");			f->SetParameter(6, 5.);		f->SetParLimits(6, 1.e-9, 1.e3);
					f->SetParName(7, "#mu_darkcount");		f->SetParameter(7, .1);		f->SetParLimits(7, 1.e-9, 1.);
					f->SetParName(8, "N_{0}_darkcount");	f->SetParameter(8, .05);	f->SetParLimits(8, 1.e-9, .3);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
				else { // default SiPM fit function
					Fitf fitf;
					int n_par = 7;
					TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, n_par); f->SetLineColor(3); f->SetNpx(1000);

					f->SetParName(0, "N_{0}");				f->SetParameter(0, his->Integral());
					f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
					f->SetParName(2, "#lambda");			f->SetParameter(2, .04);
					f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);
					f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);
					f->SetParName(5, "Gain");				f->SetParameter(5, 30.); //f->FixParameter(5, 40.);
					f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);

					if (!PrintChargeSpectrum_pars.empty()) for (int j = 0; j < min(n_par, static_cast<int>(PrintChargeSpectrum_pars.size())); j++) f->SetParameter(j, PrintChargeSpectrum_pars[j]);

					if (i < max_channel_nr_to_fit) {
						cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
						TFitResultPtr fresults = his->Fit(f, "LRS");
						fit_results.push_back(fresults);
					}
				}
			}
			TString name(Form("ChargeSpectrum channel_%02d_%d", active_channels[i], PrintChargeSpectrum_cnt));
			root_out->WriteObject(his, name.Data());
			his->Draw();
		}
	}

	SetRangeCanvas(chargec, rangestart, rangeend);
	root_out->WriteObject(chargec, ("ChargeSpectra" + to_string(PrintChargeSpectrum_cnt)).c_str());
}
/// @example timing_example.cc
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Calculate (SiPM) dark count rate
/// 
/// See PrintChargeSpectrum() for parameters.
/// 
/// @param threshold 1.5 photoelectron threshold
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
			auto his = (TH1F*)(Getwf(j * nchannels + channel_index))->Clone();
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
			delete his;
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
/// @example read_exampledata.cc
/// @example read_exampledata.py

/// @brief Finds maximum amplitude for a given channel in time window ["from", "to"] and creates 3d map of waveforms ordered by maxima
/// 
/// Use PrintMaxDist() to plot all channels.
/// Use only for small datasets as it will contain all of the data
/// 
/// @return TGraph2D of all histograms ordered by maximum amplitude
TGraph2D* ReadRun::MaxDist(int channel_index, float from, float to) {
	// find maximum amplitude for a given channel in time window [from, to] and return 3d histogram with the number of bins nbinsy,z

	TString name(Form("maxdist_ch%02d", active_channels[channel_index]));
	TGraph2D* g3d = new TGraph2D((binNumber + 2) * nevents);
	g3d->SetTitle("waveforms; t [ns]; max. amplitude [mv]; amplitude [mV]");

	for (int j = 0; j < nevents; j++) {
		if (!skip_event[j]) {
			auto his = (TH1F*)(Getwf(j * nchannels + channel_index))->Clone();
			if (from >= 0 && to > 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
			double max = his->GetMaximum();
			for (int i = 0; i < binNumber; i++) g3d->SetPoint(j * binNumber + i, his->GetXaxis()->GetBinCenter(i), max, his->GetBinContent(i));
			delete his;
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

	TCanvas* max_dist_c = new TCanvas("wf grouped by maximum", "wf grouped by maximum", 600, 400);
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
TH1F* ReadRun::His_GetTimingCFD(int channel_index, float rangestart, float rangeend, int nbins) {

	if (nbins == -999) nbins = static_cast<int>((rangeend - rangestart) / SP);

	TString name(Form("GetTimingCFD_ch%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) if (!skip_event[j]) h1->Fill(timing_results[j * nchannels + channel_index][1]);
	return h1;
}

/// @brief Plot results of GetTimingCFD()
/// 
/// \image html Print_GetTimingCFD.png "Beginning of the signals for all good events determined with constant fraction discrimination. The red lines are gauss functions fitted to the distrutions. Code in example." width=75%
/// 
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1: fits a gaussian. \n
/// Else do not fit. \n 
/// Fit results per channel are stored in ReadRun::timing_fit_results.
/// @param nbins Number of bins for histogram. Will use 320 MHz sampling rate for binning if nbins = -999.
/// @param fitoption ROOT fit option, default is "S".
/// @param set_errors Assign errors to the bins. Will assign errors of 1 to empty bins and \f$ \sqrt(N) \f$ if they are not empty. 
/// Can improve the chi^2 fit.
void ReadRun::Print_GetTimingCFD(float rangestart, float rangeend, int do_fit, int nbins, string fitoption, bool set_errors) {

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

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
			his = His_GetTimingCFD(i, rangestart, rangeend, nbins);
			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("time [ns]");

			if (set_errors) {
				for (int i = 1; i <= his->GetNbinsX(); i++) {
					if (his->GetBinContent(i) < 2) his->SetBinError(i, 1);
					else his->SetBinError(i, sqrt(his->GetBinContent(i)));
				}
			}

			his->Draw();

			if (do_fit == 1) {
				TFitResultPtr fresults = his->Fit("gaus", fitoption.c_str(), "same");
				timing_fit_results.push_back(fresults);
			}

			TString name_save(Form("Timing_cfd_channel_%02d", active_channels[i]));
			root_out->WriteObject(his, name_save.Data());
		}
	}

	timing_cfd_c->Update();
	root_out->WriteObject(timing_cfd_c, "TimingCFD");
}
/// @example timing_example.cc

/// @brief Plot timing difference between the mean timings of two channel ranges
/// 
/// See Print_GetTimingCFD_diff() for parameters.
/// 
/// @return Histogram with event-wise timing differences between two channel ranges
TH1F* ReadRun::His_GetTimingCFD_diff(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int nbins) {

	if (nbins == -999) nbins = static_cast<int>((rangeend - rangestart) / SP);

	stringstream name;
	name << "GetTimingCFD_diff <";

	// find channel indices and assemble title
	int counter = 0;
	for (int& entry : channels2) {
		if (counter > 0) name << "&";
		auto chin2 = find(active_channels.begin(), active_channels.end(), entry);
		if (chin2 != active_channels.end()) {
			name << "ch" << entry;
			entry = chin2 - active_channels.begin();
		}
		else cout << "\n\n ERROR: channels2 = " << entry << " does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";
		counter++;
	}
	name << ">-<";
	counter = 0;

	for (int& entry : channels1) {
		if (counter > 0) name << "&";
		auto chin1 = find(active_channels.begin(), active_channels.end(), entry);
		if (chin1 != active_channels.end()) {
			name << "ch" << entry;
			entry = chin1 - active_channels.begin();
		}
		else cout << "\n\n ERROR: channels1 = " << entry << " does not exist in data. Check parameters for Print_GetTimingCFD_diff()\n\n";
		counter++;
	}
	name << ">";

	// fill histogram
	auto h1 = new TH1F(name.str().c_str(), name.str().c_str(), nbins, rangestart, rangeend);
	for (int j = 0; j < nevents; j++) {
		if (!skip_event[j]) {
			float mean1 = 0., mean2 = 0., cnt1 = 0., cnt2 = 0.;
			for (int i : channels1) {
				mean1 += timing_results[j * nchannels + i][1];
				cnt1 += 1.;
			}
			for (int i : channels2) {
				mean2 += timing_results[j * nchannels + i][1];
				cnt2 += 1.;
			}
			h1->Fill(mean2 / cnt2 - mean1 / cnt1);
		}
	}

	return h1;
}

/// @brief Plot timing difference between the mean timings of two channel ranges
/// 
/// Plots the difference between the peak times between the mean times of two ranges of channels for each event. \n
/// It calculates \f$ \Delta t = <t_{second,i}> - <t_{first,i}> \f$ . \n \n \n
/// 
/// The vectors of channels to compare are added with curly brackets:
/// > mymeas.Print_GetTimingCFD_diff({ 26, 14 }, { 19 }, 0, 20, 2, 200); \n
/// would plot \f$ \Delta t = t_{ch19} - (t_{ch26} + t_{ch14})/2 \f$ from 0 ns to 20 ns with 200 bins (100 ps bin width). Another example is given in the plot below.
/// 
/// \image html Print_GetTimingCFD_diff.png "Event-wise time differences of the start of the signals of two channels. Code in example." width=75%
/// 
/// @param channels1 Vector of first channel numbers (wavecatcher channel numbers). 
/// @param channels2 Vector of second channel numbers to compare. 
/// @param rangestart Start of x range for plot in ns.
/// @param rangeend End of x range for plot in ns.
/// @param do_fit If 1: Fit a gaussian. \n
/// If 2: Fit a gaussian-exponential convolution to account for different arrival times of photons due to different possible light paths in the scintillator/light guide \n
/// and/or delay due to self-absorption and reemission of photons in the scintillator. \n
/// To be used for long light paths in the scintillator. See https://doi.org/10.1016/S0029-554X(79)90170-8 . \n
/// This option only works for sufficient asymmetry \f$\tau > \sigma/2\f$. Otherwise, the exponential decay time becomes too small to be fitted. In this case please use option 1.\n
/// If 3: Fits the sum of two gaussians where the second gauss serves as a rough background estimate. Background means events that should have been filtered out. \n
/// Else: Do not fit. \n 
/// @param nbins Number of bins for histogram.
/// @param fitrangestart Start of fitting range.
/// @param fitrangeend End of fitting range.
/// @param fitoption ROOT fitting option. Default is "RS" (chi^2). You can try to use the likelihood method with "LRS" if the data is very clean. 
/// Outliers will influence the results for the likelihood method so it is advisable to limit the fit range to exclude outliers for "LRS".
/// @param set_errors Assign errors to the bins. Will assign errors of 1 to empty bins and \f$ \sqrt(N) \f$ if they are not empty. 
/// Can improve the chi^2 fit.
void ReadRun::Print_GetTimingCFD_diff(vector<int> channels1, vector<int> channels2, float rangestart, float rangeend, int do_fit, int nbins, float fitrangestart, float fitrangeend, string fitoption, bool set_errors) {

	// call GetTimingCFD() in case it was not initialized
	if (static_cast<int>(timing_results.size()) == 0) GetTimingCFD();

	if (fitrangestart == -999) {
		fitrangestart = rangestart;
		fitrangeend = rangeend;
	}

	//gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	TCanvas* timing_cfd_d_c = new TCanvas("timing of cfd diff", "timing of cfd diff", 600, 400);

	TH1F* his;
	his = His_GetTimingCFD_diff(channels1, channels2, rangestart, rangeend, nbins);
	his->GetYaxis()->SetTitle("#Entries");
	his->GetXaxis()->SetTitle("time [ns]");

	if (set_errors) {
		for (int i = 1; i <= his->GetNbinsX(); i++) {
			if (his->GetBinContent(i) < 2) his->SetBinError(i, 1);
			else his->SetBinError(i, sqrt(his->GetBinContent(i)));
		}
	}

	his->Draw();

	double skewness = his->GetSkewness();

	if (do_fit == 1 || (do_fit == 2 && abs(skewness) < .15)) {
		// gauss (default)
		TFitResultPtr fresults = his->Fit("gaus", fitoption.c_str(), "same", fitrangestart, fitrangeend);
		timing_fit_results.push_back(fresults);
		if (do_fit == 2) cout << "\nWARNING: Print_GetTimingCFD_diff\nFITTING GAUSS INSTEAD OF GAUSS x EXP CONVOLUTION BC SYMMETRY\n";
	}
	else if (do_fit == 2) {
		// gauss x exp convolution (effective delay from random light path and/or self-absorption and reemission)
		string gxe = "[3]/(2*TMath::Abs([0]))*TMath::Exp(([1]*[1]+2*[2]*[0]-2*[0]*x)/(2*[0]*[0]))*TMath::Erfc(([1]*[1]+[0]*([2]-x))/(1.4142*TMath::Abs([0])*[1]))";
		auto expgconv = new TF1("exp x gauss convolution", gxe.c_str(), fitrangestart, fitrangeend);
		expgconv->SetNpx(5000);

		// this parameter describes the sigma from different light paths 
		// and/or the effective decay time constant for self-absorption and reemission
		expgconv->SetParName(0, "#tau_{eff}");		expgconv->SetParameter(0, skewness);
		if (skewness>0) expgconv->SetParLimits(0, .15, 5.);
		else expgconv->SetParLimits(0, -5., -.15);
		//expgconv->FixParameter(0, 1.55);

		expgconv->SetParName(1, "#sigma_{gaus}");		expgconv->SetParameter(1, his->GetStdDev());
		expgconv->SetParLimits(1, 1e-1, 7.);	//expgconv->FixParameter(1, .7);

		expgconv->SetParName(2, "t_{0}");		expgconv->SetParameter(2, his->GetMean());
		expgconv->SetParLimits(2, fitrangestart, fitrangeend);	//expgconv->FixParameter(2, 6.6);

		expgconv->SetParName(3, "norm");		expgconv->SetParameter(3, his->Integral("width"));
		expgconv->SetParLimits(3, 1., 1e8);		//expgconv->FixParameter(3, 105.5);

		TFitResultPtr fresults = his->Fit(expgconv, "SR", "same");
		timing_fit_results.push_back(fresults);

		// for the phi_ew-analysis: print out the time value of the maximum of the best fit --> used to determine timing cuts
		float t_of_maximum = expgconv->GetMaximumX(-5, 5);
		cout << "Maximum of the fit is at t=" << t_of_maximum << " ns" << endl;

		TLine* mean = new TLine(expgconv->GetParameter(2), 1e-2, expgconv->GetParameter(2), his->GetMaximum());
		mean->SetLineColor(1); mean->SetLineWidth(2);
		mean->Draw("same");
	}
	else if (do_fit == 3) {
		// sum of two gaussians (one as background estimate)
		auto two_gauss = new TF1("two gaussians", "gaus(0)+gaus(3)", rangestart, rangeend);
		two_gauss->SetTitle("Sum of two gauss");
		float posmax = his->GetXaxis()->GetBinCenter(his->GetMaximumBin());
		two_gauss->SetParameters(his->Integral("width"), posmax, 0.35, his->Integral("width") / 30, posmax, 2);
		two_gauss->SetParName(0, "norm_{peak}");		two_gauss->SetParName(1, "#mu_{peak}");			two_gauss->SetParName(2, "#sigma_{peak}");			two_gauss->SetParLimits(2, 1e-9, 1e2);
		two_gauss->SetParName(3, "norm_{background}");	two_gauss->SetParName(4, "#mu_{background}");	two_gauss->SetParName(5, "#sigma_{background}");	two_gauss->SetParLimits(5, 1e-9, 1e2);
		TFitResultPtr fresults = his->Fit(two_gauss, fitoption.c_str(), "same", fitrangestart, fitrangeend);
		timing_fit_results.push_back(fresults);
	}

	root_out->WriteObject(his, his->GetTitle());
	timing_cfd_d_c->Update();
	root_out->WriteObject(timing_cfd_d_c, "TimingCFD_diff");
}
/// @example timing_example.cc






/// @brief Helper. Creates a list of .bin data files in data folder to be read in
/// @param dirname Directory
/// @param ext File extension
/// @return String of line separated file names
string ReadRun::ListFiles(const char* dirname, const char* ext) {

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
			}
		}
		TIter next2(files);
		while ((file = (TSystemFile*)next2())) {
			fname = file->GetName();
			if (!file->IsDirectory() && !fname.EndsWith(ext) && fname.Contains(ext)) {
				ss << fname.Data() << "\n";
			}
		}
	}
	return ss.str();
}

/// @brief Helper that returns the waveform histogram for a certain waveform number number 
/// @param wf_nr Waveform number
/// @return Waveform histogram 
TH1F* ReadRun::Getwf(int wf_nr) {
	return (TH1F*)rundata->At(wf_nr);
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
	double* xvals = new double[binNumber];
	for (int i = 0; i < binNumber; i++) {
		xvals[i] = static_cast<double>(SP) * static_cast<double>(i) + shift;
	}
	return xvals;
}

/// @brief Get array of y values for a certain waveform
/// @param waveform_index Waveform index
/// @return Y values of waveform
double* ReadRun::gety(int waveform_index) {
	TH1F* his = Getwf(waveform_index);
	double* yvals = new double[binNumber];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i + 1);
	}
	return yvals;
}

/// @brief Get array of y values for a certain waveform
/// @param channelnr Channel number index (not the actual channel number)
/// @param eventnr Event number
/// @return Y values of waveform
double* ReadRun::gety(int channelnr, int eventnr) {
	TH1F* his = Getwf(channelnr, eventnr);
	double* yvals = new double[binNumber];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i + 1);
	}
	return yvals;
}

/// @brief Get array of y values for a certain waveform
/// @param his Waveform histogram
/// @return Y values of waveform
double* ReadRun::gety(TH1F* his) {
	double* yvals = new double[binNumber];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i + 1);
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
		yvals[i - start_at] = his->GetBinContent(i + 1);
	}
	return yvals;
}

/// @brief Translate a random number into a useful root color https://root.cern.ch/doc/master/classTColor.html
/// @param i Index of your plotting loop that is to be translated into a useful ROOT color index
/// @return ROOT color index
int ReadRun::rcolor(unsigned int i) {
	const int nclrs = 17;
	int rclrs[nclrs] = { 1, 2, 3, 4, 6, 8, 9, 13, 20, 28, 30, 34, 38, 40, 31, 46, 49 };
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

/// @brief  Match channel number (wavecatcher input channel) to channel index
/// @param channel_number Number of the channel as defined in the wavecatcher software
/// @return Corresponding index for this channel
int ReadRun::GetChannelIndex(int channel_number) {
	int channel_index = -1;
	for (int i = 0; i < static_cast<int>(active_channels.size()); i++) {
		if (active_channels[i] == channel_number) channel_index = i;
	}
	if (channel_index == -1) {
		cout << "\n\n\tERROR: channel " << channel_number << " does not exist in data. Will continue with first channel\n\n";
		channel_index = 0;
	}
	return channel_index;
}

/// @brief Get the current channel index for a certain waveform index
/// @param waveform_index 
/// @return Current channel index
int ReadRun::GetCurrentChannel(int waveform_index) {
	return (waveform_index - nchannels * floor(waveform_index / nchannels));
}

/// @brief Get the current event index for a certain waveform index
/// @param waveform_index 
/// @return Current event index
int ReadRun::GetCurrentEvent(int waveform_index) {
	return floor(waveform_index / nchannels);
}

/// @brief Helper to split canvas according to the number of channels to be plotted
/// @param c Canvas to be split
void ReadRun::SplitCanvas(TCanvas*& c) {
	// cross check if user input exists in data
	vector<int> rmv;
	for (int i = 0; i < static_cast<int>(plot_active_channels.size()); i++) {
		if (find(active_channels.begin(), active_channels.end(), plot_active_channels[i]) == active_channels.end()) {
			cout << "\n\n\n ------------ WARNING ------------\n";
			cout << "YOUR SELECTED CHANNEL " << plot_active_channels[i] << " DOES NOT EXIST IN DATA\n";
			cout << "PLEASE CHANGE plot_active_channels\n\n\n";
			rmv.push_back(plot_active_channels[i]);
		}
	}

	for (int i = 0; i < static_cast<int>(rmv.size()); i++) {
		auto it = find(plot_active_channels.begin(), plot_active_channels.end(), rmv[i]);
		if (it != plot_active_channels.end()) plot_active_channels.erase(it);
	}

	if (plot_active_channels.empty()) c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else if (static_cast<int>(plot_active_channels.size()) > 1) c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
}

/// @brief Set consistent x-axis and y-axis range for all TH1 histograms on a canvas. 
/// 
/// Will only change the axes where ```min != max```.
/// 
/// @param c Canvas Canvas
/// @param x_range_min X-axis minimum
/// @param x_range_max X-axis maximum
/// @param y_range_min Y-axis minimum
/// @param y_range_max Y-axis maximum
void ReadRun::SetRangeCanvas(TCanvas*& c, double x_range_min, double x_range_max, double y_range_min, double y_range_max) {

	// Lambda function to set the axis ranges
	auto setAxisRanges = [&](TH1 *his) {
		if (x_range_min != x_range_max) his->GetXaxis()->SetRangeUser(x_range_min, x_range_max);
		if (y_range_min != y_range_max) his->GetYaxis()->SetRangeUser(y_range_min, y_range_max);
	};

	// Get the list of pads on the canvas and loop over pads
	TList* pads = c->GetListOfPrimitives();
	TIter nextPad(pads);
	TObject* object;
	while ((object = nextPad())) {
		if (object->InheritsFrom(TPad::Class())) {
			TPad* pad = static_cast<TPad*>(object);

			// Set the axis ranges for all plots on the current pad
			pad->cd();
			TList* primitives = pad->GetListOfPrimitives();
			TIter nextPrimitive(primitives);
			TObject* primitive;
			while ((primitive = nextPrimitive())) {
				if (primitive->InheritsFrom(TH1::Class())) {
					//TH1* his = static_cast<TH1*>(primitive);
					setAxisRanges(static_cast<TH1*>(primitive));
				}
			}
		}
		else if (object->InheritsFrom(TH1::Class())) {
			//TH1* his = static_cast<TH1*>(object);
			setAxisRanges(static_cast<TH1*>(object));
		}
	}

	c->Modified();
	c->Update();
}

/// @brief Simple linear interpolation for x
/// 
/// Will return the mean of x1 and x2 if y1=y2.
/// 
/// @param ym Y value for evaluation
/// @param x1 X1 
/// @param x2 X2
/// @param y1 Y1
/// @param y2 Y2
/// @return x value at "ym"
float ReadRun::LinearInterpolation(float ym, float x1, float x2, float y1, float y2) {
	if (y1 == y2) return (x1 + x2) / 2.;
	else if (y1 > ym) { 
		cout << "\nError in LinearInterpolation:\value ym=" << ym << " out of range" << endl;
		return 1;
	}
	else return x1 + (ym - y1) * (x2 - x1) / (y2 - y1);
}

/// @brief Helper to perform convolution of two 1D arrays
/// 
/// Used for smoothing etc.
/// 
/// @param[in,out] result Array containing convolution result
/// @param first First array for convolution
/// @param second Second array for convolution 
/// @param size Size of arrays
void ReadRun::Convolute(double*& result, double* first, double* second, int size) {

	// FFT real -> im
	auto fftfirst = [&size](double*& orig, double*& re, double*& im) {
		TVirtualFFT* fft = TVirtualFFT::FFT(1, &size, "R2C ES");
		fft->SetPoints(orig);
		fft->Transform();
		fft->GetPointsComplex(re, im);
		delete fft;
	};

	double* refirst = new double[size];
	double* imfirst = new double[size];
	double* resecond = new double[size];
	double* imsecond = new double[size];
	double* reres = new double[size];
	double* imres = new double[size];

	fftfirst(first, refirst, imfirst);
	fftfirst(second, resecond, imsecond);

	TComplex cofirst;
	TComplex cosecond;
	TComplex cores;

	for (int i = 0; i < size; i++) {
		cofirst(refirst[i], imfirst[i]);
		cosecond(resecond[i], imsecond[i]);

		cores = cofirst * cosecond / static_cast<double>(size);

		reres[i] = cores.Re();
		imres[i] = cores.Im();
	}


	//cout << "performing IFFT ... ";
	TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &size, "C2R ES");
	fft_back->SetPointsComplex(reres, imres);
	fft_back->Transform();
	fft_back->GetPoints(result);
	delete fft_back;
	delete[] imres; delete[] reres; delete[] refirst; delete[] imfirst; delete[] resecond; delete[] imsecond;
}

/// @brief Apply smoothing array of double with length nbins
/// 
/// Use with care. Method=2 is preferred. \n \n
///
/// Please note that if you want to use gaussian smoothing for data with a binning different from 0.3125 ns/bin 
/// you need to set the variable bin_size to the new bin size.
/// 
/// \image html use_functions_wo_measurement.png "Gaussian smoothing of a simple array with 15 entries. Code in example." width=50%
/// 
/// @param[in,out] ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Number of bins before and after central bin for running average OR gauss sigma in ns for gauss kernel and convolution (see parameter bin_size).
/// @param method If 0 use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1 use 5 sigma gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else use 3 sigma gaussian kernel smoothing. Preferred method, fast.
/// @param bin_size Bin width of the array to smooth for gauss sigma. Default is .3125 for wavecatcher sampling rate. 
/// Set to 1 to change sigma unit to number of bins.
void ReadRun::SmoothArray(double*& ar, int nbins, double sigma, int method, double bin_size) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	if (method == 0) {
		// calculate running average from -sigma until +sigma (sigma = number of bins)
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
	else if (method == 1) {
		// convolution with gauss clipped at +-5 sigma (very inefficient and slow)
		double* gauss = new double[nbins];

		double sum = 0.;

		for (int i = 0; i < nbins; i++) {
			if (static_cast<double>(i) * bin_size < 5 * sigma) gauss[i] = 
				TMath::Exp(-1. * TMath::Power((static_cast<double>(i) * bin_size - 5 * sigma), 2.) / (2. * sigma * sigma)) / (sigma * 2.506628);
			else gauss[i] = 0.;
			sum += gauss[i];
		}

		for (int i = 0; i < nbins; i++) {
			gauss[i] /= sum;
		}

		Convolute(ar, artmp, gauss, nbins);
		delete[] gauss;
	}
	else {
		// gauss kernel 3*sigma
		int nbins_3sigma = static_cast<int>(ceil(6. * sigma / bin_size));
		if (nbins_3sigma % 2 == 0) nbins_3sigma++;
		if (nbins_3sigma > 1) {
			double* gauss = new double[nbins_3sigma];
			double gauss_offset = floor(static_cast<double>(nbins_3sigma) / 2.) * bin_size;
			double denom = 2. * sigma * sigma;
			for (int i = 0; i < nbins_3sigma; i++) {
				gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i)) * bin_size - gauss_offset, 2.) / denom);
			}

			double res = 0;
			double norm = 0;
			for (int i = 0; i < nbins; i++) {
				res = 0.;
				norm = 0.;
				for (int j = max(0, nbins_3sigma / 2 - i); j < min(nbins - i + nbins_3sigma / 2, nbins_3sigma); j++) {
					res += gauss[j] * artmp[i + j - nbins_3sigma / 2];
					norm += gauss[j];
				}
				if (norm != 0.) ar[i] = res / norm;
			}
			delete[] gauss;
		}
	}
	delete[] artmp;
}
/// @example use_functions_wo_measurement.cc

/// @brief Apply filter for array of double with length nbins
/// 
/// Experimental, can be used to highlight peaks and suppress long tails (suppresses low and high frequencies). Use Filter_test.ipynb to test parameters.
/// 
/// @param[in,out] ar Array to be filtered.
/// @param nbins Number of bins of input.
/// @param sigma1 First.
/// @param sigma2 Second.
/// @param factor Factor for negative part (<=1).
/// @param bin_size Bin width. Default is .3125. Set to 1 to get sigma in units of bins.
void ReadRun::FilterArray(double*& ar, int nbins, double sigma1, double sigma2, double factor, double bin_size) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	// shifted difference of two gauss functions (~smoothed derivative)
	int nbins_2sigma = static_cast<int>(ceil((2. * sigma1 + 3. * sigma2) / bin_size));
	double* sdog = new double[nbins_2sigma];

	double denom1 = 2. * sigma1 * sigma1;
	double denom2 = 2. * sigma2 * sigma2;
	for (int i = 0; i < nbins_2sigma; i++) {
		sdog[i] = TMath::Exp(-1. * TMath::Power(static_cast<double>(i) * bin_size - 3. * sigma2, 2.) / denom1) -
			factor * TMath::Exp(-1. * TMath::Power(static_cast<double>(i) * bin_size - 2. * sigma2, 2.) / denom2);
	}

	double res = 0;
	double norm = 0;
	for (int i = 0; i < nbins; i++) {
		res = 0.;
		norm = 0.;
		for (int j = max(0, nbins_2sigma / 2 - i); j < min(nbins - i + nbins_2sigma / 2, nbins_2sigma); j++) {
			res += sdog[j] * artmp[i + j - nbins_2sigma / 2];
			if (sdog[j] > 0.) norm += sdog[j];
		}
		if (norm != 0.) ar[i] = res / norm;
	}
	delete[] sdog;
}