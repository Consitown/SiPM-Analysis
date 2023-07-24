#include "CosmicsBox.h"

/// @brief Angular distribution of passing particles in the cosmics setup 
/// 
/// Plots distribution of the event-dependent angle phi_ew in a histogram. \n 
/// The angle is calculated using vectorial addition and the lightyield of each channel in each event. \n
/// For more info, see master thesis of Alexander Vagts. \n
/// Currently uses the uncorrected formulas.
/// 
/// \image html phi-ew.png "Example fit of the reconstructed angles for a certain event selection in one corner of the cosmics box. Code in example phi-ew.ipynb." width=50%
/// 
/// @param phi_chx Vector of the associated angles of each channel. Order is important. Angle of the first channel first. \n
/// E.g. if channel 0 is used, the first angle would be from channel 0
/// @param ly_C0 Vector of the lightyield of all channels at C0 position for the correction. Order must be the same as in phi_chx
/// @param SiPMchannels all SiPM channels which should be included in analysis, e.g. {0, 2, 4, 6}
/// @param windowmin left edge of integration window for lightyield
/// @param windowmax right edge of integration window for lightyield
/// @param maxfrom searches for peak from this to
/// @param maxto this
/// @param nbins Number of bins in histogram
/// @param corr selection bool for corrected or uncorrected spectra \n
/// If true - corrected spectra \n
/// If false - uncorrected spectra
/// @param periodic If true, will print all phi_ew shifted by +/- 360 deg (so normal phi_ew distri * 3) and fit a periodic gauss
/// @return Phi_ew spectrum
void CosmicsBox::Print_Phi_ew(vector<int> phi_chx, vector<float> ly_C0, vector<int> SiPMchannels, float windowmin, float windowmax, float maxfrom, float maxto, int nbins, bool corr, bool periodic) {

	// plotting uncorrected/corrected phi_ew - spectra ; first: map the phi_i to the channels, like Alex did, but with new channel positions
	// match channel number to channel index (still very specific)
	int sipmnum = SiPMchannels.size();
	vector<int> ch_index; for (int i = 0; i < sipmnum; i++) ch_index.push_back(0); // initialize the ch_index vector
	for (int i = 0; i < static_cast<int>(active_channels.size()); i++) {
		for (int j = 0; j < sipmnum; j++) {
			if (active_channels[i] == SiPMchannels[j]) ch_index[j] = i;
		}
	}

	// compute correction factor
	vector<float> ly_corr;
	if (corr) {
		float ly_av = 0;
		// average lightyield of all channels * sipmnum
		for (int i = 0; i < sipmnum; i++) ly_av += ly_C0[i];
		// divide ly of one channel by (1/sipmnum)*sipmnum*average ly of all channels
		for (int i = 0; i < sipmnum; i++) ly_corr.push_back(sipmnum * ly_C0[i] / ly_av);
	}
	else for (int i = 0; i < sipmnum; i++) ly_corr.push_back(1); // no correction

	// print correction factors
	for (int i = 0; i < sipmnum; i++) cout << "Correction factor for channel " << SiPMchannels[i] << ":" << ly_corr[i] << endl;

	// initialize canvas + histograms
	gStyle->SetOptStat("nemr"); gStyle->SetOptFit(1111); //draws a box with some histogram parameters
	TString his_name(Form("#phi_ew-spectrum_from_ch%d_to_ch%d", SiPMchannels.front(), SiPMchannels.back()));
	TCanvas* hisc = new TCanvas(his_name, his_name, 600, 400);

	double min_angle = -180, max_angle = 180; // 360 deg plot range
	if (periodic) { min_angle = -540; max_angle = 540; } // 1080 deg plot range
	TH1* his = new TH1F(his_name, his_name, nbins, min_angle, max_angle);

	// loop through all events and compute phi_ew
	float lightyield, anglevaluex, anglevaluey, phi_ew = 0;
	for (int i = 0; i < nevents; i++) {
		if (!skip_event[i]) {
			for (int j = 0; j < sipmnum; j++) { //loop through all SiPM-channels
				TH1F* hisly = ((TH1F*)rundata->At(i * nchannels + ch_index[j]));
				lightyield = GetPeakIntegral(hisly, windowmin, windowmax, maxfrom, maxto, 0); //lightyield as the integral around maximum
				anglevaluex += cos(phi_chx[j] * TMath::Pi() / 180) * lightyield / ly_corr[j]; //x part of vectorial addition
				anglevaluey += sin(phi_chx[j] * TMath::Pi() / 180) * lightyield / ly_corr[j]; //y part of vectorial addition
			}

			phi_ew = atan2(anglevaluey, anglevaluex) * 180 / TMath::Pi(); //vectorial addition for phi_ew --> fill in histo
			his->Fill(phi_ew);
			if (periodic) {
				his->Fill(phi_ew + 360);
				his->Fill(phi_ew - 360);
			}
			anglevaluex = 0, anglevaluey = 0; //reset the x and y parts
		}
	}

	//make histogram fancy + printing
	his->GetXaxis()->SetTitle("#phi_ew (deg.)"); his->GetYaxis()->SetTitle("#Entries"); //titling of axes
	his->Draw();

	if (periodic) {
		Fitf_periodic_gauss fitf;
		int n_par = 4;
		TF1* f = new TF1("fitf", fitf, min_angle, max_angle, n_par); f->SetLineColor(2); f->SetNpx(1000);

		double max = his->GetMaximum();
		double min = his->GetMinimum();
		his->GetXaxis()->SetRangeUser(-180, 180); //only look at one period for starting values
		double phi_est = his->GetXaxis()->GetBinCenter(his->GetMaximumBin());
		his->GetXaxis()->SetRangeUser(0, 0); //reset scale

		f->SetParName(0, "A");				f->SetParameter(0, max - min);		f->SetParLimits(0, 1, 1e9);
		f->SetParName(1, "#Phi_{ew}");		f->SetParameter(1, phi_est);		f->SetParLimits(1, -180, 180);
		f->SetParName(2, "#sigma");			f->SetParameter(2, 40);				f->SetParLimits(2, 5, 360);
		f->SetParName(3, "offset");			f->SetParameter(3, min);			f->SetParLimits(3, TMath::Min(min - 1, 0.1), 1e9);

		TFitResultPtr fresults = his->Fit(f, "LRS");
		fit_results.push_back(fresults);
	}

	hisc->Update();
	root_out->WriteObject(hisc, "Phi_ew");
	root_out->WriteObject(his, "Phi_ew_his");

	string pdf_name = "phi_ew_spectrum";
	if (corr) pdf_name += "_corr.pdf"; else pdf_name += "_uncorr.pdf";
	hisc->SaveAs(pdf_name.c_str()); //write the histogram to a .pdf-file
}