#include "FFT_WF.h"

/// @brief Print Fourier transform of waveforms
///  
/// For analysis of noise/signal frequencies, potential filters, or input data for machine learning.
/// 
/// @param eventnr Event number
/// @param xmin Range min
/// @param xmax Range max
/// @param multiplier Multiplier for resolution of plot
void FFT_WF::PrintFFTWF(int eventnr, float xmin, float xmax, int multiplier) {
	// plot waveforms of all channels for a given event number eventnr and add the determined integration windows to the plot
	TString name(Form("fft_waveforms_event__%04d", eventnr));
	TCanvas* fftc = new TCanvas(name.Data(), name.Data(), 600, 400);
	SplitCanvas(fftc);

	TString imname(Form("fft_im_waveforms_event__%04d", eventnr));
	TCanvas* imfftc = new TCanvas(imname.Data(), imname.Data(), 600, 400);
	SplitCanvas(imfftc);

	int size = binNumber * multiplier;
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
			if (j < binNumber) yvals[j] = his->GetBinContent(j + 1);
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