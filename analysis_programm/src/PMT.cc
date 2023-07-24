#include "PMT.h"

/// @brief "Charge" spectrum optimized for PMT signals
/// 
/// Just for plotting. To analyze the data use PrintChargeSpectrum() with Fitf_PMT_pedestal() for low number of photons 
/// and Fitf_langaus() for >10-15 photons. \n 
/// See PrintChargeSpectrum() for parameters.
/// 
void PMT::PrintChargeSpectrumPMT(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {

	PrintChargeSpectrumPMT_cnt++;

	string ctitle("charge spectra PMT" + to_string(PrintChargeSpectrumPMT_cnt));
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
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
			two_gauss->SetParName(2, "#sigma_{pedestal}");	two_gauss->SetParLimits(2, 1e-9, 1e3);
			two_gauss->SetParName(3, "A_{SPE}");
			two_gauss->SetParName(4, "#mu_{SPE}");
			two_gauss->SetParName(5, "#sigma_{SPE}");		two_gauss->SetParLimits(5, 1e-9, 1e3);

			if (!PrintChargeSpectrumPMT_pars.empty()) {
				for (int j = 0; j < static_cast<int>(PrintChargeSpectrumPMT_pars.size()); j++) two_gauss->SetParameter(j, PrintChargeSpectrumPMT_pars[j]);
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
void PMT::PrintChargeSpectrumPMTthreshold(float windowlow, float windowhi, float rangestart, float rangeend, int nbins, double threshold, bool calculate_SiPM_DCR) {

	PrintChargeSpectrumPMTthreshold_cnt++;

	gStyle->SetOptStat(0); // 11 is title + entries

	// show fraction of events above 0.5 pe charge = pedestal + gain/2
	// dark count rate for SiPMs (currently only automated for fit function Fitf)
	// need to call the SiPM fit function before this one for this functionality
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
	TCanvas* chargec = new TCanvas(ctitle.c_str(), ctitle.c_str(), 600, 400);
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