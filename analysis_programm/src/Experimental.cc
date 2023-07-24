#include "Experimental.h"

/// @brief Rebin the data to test bandwidth effects. Will combine an integer number of bins into a new, 
/// wider bin and divide the new bin content by the integer number to preserve the shape and integral. 
/// See [ROOT::TH1::Rebin()](https://root.cern.ch/doc/master/classTH1.html).
/// 
/// CAUTION: Not tested with all functions. Make sure to adjust bin numbers e. g. in 
/// CorrectBaselineMinSlopeRMS() and bin size e. g. in SmoothArray(). 
/// 
/// Please note that PlotChannelSums() is calculated while parsing the data - **before** rebinning. 
/// Use PlotChannelAverages() instead.
/// 
/// @param ngroup Integer number of bins to combine.
/// @param noise_level Add gaussian noise with ```sigma = noise_level``` to rebinned data.
/// @param seed Seed for the random number generator.
void Experimental::RebinAll(int ngroup, float noise_level, unsigned long seed) {

	SP *= static_cast<float>(ngroup);
	binNumber /= ngroup;
	float norm = 1. / static_cast<float>(ngroup);
	cout	<< "\nRebinning the data to a new sampling rate of " << 1. / SP 
			<< " GS/s which corresponds to a bin size of " << SP << " ns and the data now has " << binNumber << " bins\n";

	for (int j = 0; j < nwf; j++) {
		TH1F* his = Getwf(j);
		his->Rebin(ngroup);
		
		his->Scale(norm);

		if (noise_level != 0.) {
			auto noise = new TRandom3();
			noise->SetSeed(seed);
			for (int i = 1; i <= his->GetNbinsX(); i++) his->SetBinContent(i, his->GetBinContent(i) + noise->Gaus(0, noise_level));
		}
	}
}
/// @example timing_example_rebin.cc

/// @brief derivative of all waveforms (for measurements w/o pole-zero cancellation)
///
/// Experimental!
void Experimental::DerivativeAll() {
	// just for testing
	cout << "\nderivative of wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = Getwf(j);
		double* yvals = gety(his);
		for (int i = 1; i <= his->GetNbinsX() - 1; i++) his->SetBinContent(i, yvals[i + 1] - yvals[i]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}