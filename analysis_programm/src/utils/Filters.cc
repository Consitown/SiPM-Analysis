#include "Filters.h"

/// @brief Helper to perform convolution of two 1D arrays
/// 
/// WARNING: All arrays must be of the same size. \n \n
/// Used for smoothing etc.
/// 
/// @param[in,out] result Array containing convolution result
/// @param first First array for convolution
/// @param second Second array for convolution 
/// @param size Size of arrays
void Filters::Convolute(double*& result, double* first, double* second, int size) {

	// FFT real -> im
	auto fft_re_im = [&size](double* orig, double*& re, double*& im) {
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

	fft_re_im(first, refirst, imfirst);
	fft_re_im(second, resecond, imsecond);

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

/// @brief Helper to perform convolution of two 1D arrays
/// 
/// WARNING: All arrays must be of the same size. \n \n
/// Used for smoothing etc.
/// 
/// @param[in,out] first Array to be convoluted
/// @param second Second array for convolution 
/// @param size Size of arrays
void Filters::Convolute(double*& first, double* second, int size) {
	Convolute(first, first, second, size);
}

/// @brief Helper to perform partial deconvolution of two 1D arrays
///
/// WARNING: All arrays must be of the same size. \n \n
/// 
/// @param[in,out] result Array containing deconvolution result
/// @param first First array for deconvolution
/// @param second Second array for deconvolution
/// @param size Size of arrays
/// @param sigma Gaussian sigma in ns
/// @param x0 Gaussian center in ns
/// @param bin_size Bin size in ns
void Filters::Deconvolute(double*& result, double* first, double* second, int size, double sigma, double x0, double bin_size) {

	// FFT real -> im
	auto fft_re_im = [&size](double* orig, double*& re, double*& im) {
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
	double* regauss = new double[size];
	double* imgauss = new double[size];
	double* reres = new double[size];
	double* imres = new double[size];
	double* gauss = new double[size];
	
	fft_re_im(first, refirst, imfirst);
	fft_re_im(second, resecond, imsecond);
	if (sigma > 0.) {
		for (int i = 0; i < size; i++) {
			gauss[i] = TMath::Gaus(static_cast<double>(i) * bin_size, x0, sigma);
		}
		fft_re_im(gauss, regauss, imgauss);
	}

	TComplex cofirst;
	TComplex cosecond;
	TComplex cogauss;
	TComplex cores;

	for (int i = 0; i < size; i++) {
		cofirst(refirst[i], imfirst[i]);
		cosecond(resecond[i], imsecond[i]);
		if (sigma > 0.) cogauss(regauss[i], imgauss[i]);
		else cogauss(1., 0.);

		cores = cofirst * cogauss / cosecond / static_cast<double>(size);

		reres[i] = cores.Re();
		imres[i] = cores.Im();
	}

	TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &size, "C2R ES");
	fft_back->SetPointsComplex(reres, imres);
	fft_back->Transform();
	fft_back->GetPoints(result);
	delete fft_back;

	delete[] imres; delete[] reres; delete[] refirst; delete[] imfirst; 
	delete[] resecond; delete[] imsecond; delete[] regauss; delete[] imgauss; delete[] gauss;
}

/// @brief Apply smoothing array of double with length nbins
/// 
/// Use with care. Method="Gaus" is preferred. \n \n
///
/// Please note that if you want to use gaussian smoothing for data with a binning different from 0.3125 ns/bin 
/// you need to set the variable bin_size to the new bin size.
/// 
/// \image html use_functions_wo_measurement.png "Gaussian smoothing of a simple array with 15 entries. Code in example." width=50%
/// 
/// @param[in,out] ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Number of bins before and after central bin for "Box" and "Median". \n
/// Gaussian sigma in ns for "Gaus" and "GausFFT". \n
/// Sigma in x for "Bilateral" and "Bilateral2".
/// @param method "Box": Use running average (box kernel smoothing). Simple, very fast. \n 
/// "GausFFT": Use 5 sigma FFT gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// "Gaus": Use 3 sigma gaussian kernel smoothing. Preferred method, fast. 
/// "Median": Use median filter. Fast. \n
/// "Bilateral": Use bilateral filter. Slow. \n
/// "Bilateral2": Use bilateral filter with different kernel. Slow. \n
/// @param bin_size Bin width of the array to smooth for gauss sigma. Default is .3125 for wavecatcher sampling rate. 
/// Set to 1 to change sigma unit to number of bins.
/// @param sigma2 Sigma in y for "Bilateral" and "Bilateral2".
void Filters::SmoothArray(double*& ar, int nbins, double sigma, string method, double bin_size, double sigma2) {

	if (method == "Box" || method == "Mean" || method == "Average")	BoxFilter(ar, nbins, static_cast<int>(sigma));
	else if	(method == "GausFFT")									GausFFTFilter(ar, nbins, sigma, bin_size);
	else if	(method == "Gaus")										GausFilter(ar, nbins, sigma, bin_size);
	else if (method == "Bilateral")									BilateralFilter(ar, nbins, sigma, sigma2, bin_size);
	else if (method == "Bilateral2")								Bilateral2Filter(ar, nbins, sigma, sigma2, .05, bin_size);
	else if (method == "Median")									MedianFilter(ar, nbins, static_cast<int>(sigma));
	else {
		cout << "Smooth method " << method << " not known. Using default Gaus method." << endl;
		GausFilter(ar, nbins, sigma, bin_size);
	}
}
/// @example use_functions_wo_measurement.cc

/// @brief Apply smoothing array of double with length nbins
/// @param method If 0: Use running average (box kernel smoothing). Simple, very fast. \n 
/// If 1: Use 5 sigma FFT gaussian smoothing. This method is not central and will shift peaks. Very slow. \n
/// Else: Use 3 sigma gaussian kernel smoothing. Preferred method, fast. 
void Filters::SmoothArray(double*& ar, int nbins, double sigma, int method, double bin_size, double sigma2) {

	switch (method) {
	case 0:
		SmoothArray(ar, nbins, sigma, "Box");
		break;
	case 1:
		SmoothArray(ar, nbins, sigma, "GausFFT", bin_size);
		break;
	case 3: 		
		SmoothArray(ar, nbins, sigma, "Bilateral", bin_size, sigma2);
		break;
	case 4:
		SmoothArray(ar, nbins, sigma, "Bilateral2", bin_size, sigma2);
		break;
	case 5: 
		SmoothArray(ar, nbins, sigma, "Median");
		break;
	default:
		SmoothArray(ar, nbins, sigma, "Gaus", bin_size);
		break;
	}
}

/// @brief Simple running average (box) filter
/// @param ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Number of bins before and after central bin.
void Filters::BoxFilter(double*& ar, int nbins, int sigma) {
	// calculate running average from -sigma until +sigma (sigma = number of bins)
	for (int i = 0; i < nbins; i++) {
		double mean = 0.;
		int nmean = 0;
		int start = max(0, i - sigma);
		int end = min(i + sigma, nbins - 1);

		for (int k = start; k <= end; k++) {
			mean += ar[k];
			nmean++;
		}
		if (nmean != 0) {
			ar[i] = mean / static_cast<double>(nmean);
		}
	}
}

/// @brief Gaussian smoothing with simple 3 sigma kernel
/// @param ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma Gauss sigma in ns.
void Filters::GausFilter(double*& ar, int nbins, double sigma, double bin_size) {
	// gauss kernel 3*sigma
	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	// calculate running average from -sigma until +sigma (sigma = number of bins)
	int nbins_6sigma = static_cast<int>(ceil(6. * sigma / bin_size));
	if (nbins_6sigma % 2 == 0) nbins_6sigma++;
	if (nbins_6sigma > 1) {
		double gauss_offset = floor(static_cast<double>(nbins_6sigma) / 2.) * bin_size;
		double gauss[nbins_6sigma];
		for (int i = 0; i < nbins_6sigma; i++) {
			gauss[i] = TMath::Gaus(static_cast<double>(i) * bin_size, gauss_offset, sigma);
		}

		double res = 0, norm = 0;
		for (int i = 0; i < nbins; i++) {
			res = 0., norm = 0.;
			int start = max(0, nbins_6sigma / 2 - i);
			int end = min(nbins - i + nbins_6sigma / 2, nbins_6sigma);

			for (int j = start; j < end; j++) {
				res += gauss[j] * artmp[i + j - nbins_6sigma / 2];
				norm += gauss[j];
			}
			if (norm != 0.) ar[i] = res / norm;
		}
	}
	delete[] artmp;
}

/// @brief Gaussian smoothing with FFT convolution
/// @param ar Array to be smoothed.
///	@param nbins Number of bins of input.
/// @param sigma Gauss sigma in ns.
void Filters::GausFFTFilter(double*& ar, int nbins, double sigma, double bin_size) {
	// convolution with gauss clipped at +-5 sigma (very inefficient and slow)
	double* gauss = new double[nbins];

	double sum = 0.;
	double five_sigma = 5 * sigma;
	double denom1 = -2. * sigma * sigma;
	double denom2 = sigma * 2.506628;

	for (int i = 0; i < nbins; i++) {
		double position = static_cast<double>(i) * bin_size;
		if (position < five_sigma) gauss[i] = exp(pow((position - five_sigma), 2) / denom1) / denom2;
		else gauss[i] = 0.;
		sum += gauss[i];
	}

	for (int i = 0; i < nbins; i++) {
		gauss[i] /= sum;
	}

	Convolute(ar, gauss, nbins);
	delete[] gauss;
}

/// @brief Bilateral filter (non-linear) with 3 sigma kernel
/// 
/// For noise reduction while preserving edges, i .e. for baseline correction. Slow.
/// 
/// @param ar Array to be smoothed.
///	@param nbins Number of bins of input.
/// @param sigma_x Gauss sigma along x in ns.
/// @param sigma_y Gauss sigma along y in mV.
/// @param bin_size Bin size in ns.
void Filters::BilateralFilter(double*& ar, int nbins, double sigma_x, double sigma_y, double bin_size) {
	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	int nbins_3sigma = static_cast<int>(ceil(3 * sigma_x / bin_size));

	for (int i = 0; i < nbins; i++) {
		double weighted_sum = 0.0;
		double sum_of_weights = 0.0;

		for (int j = max(0, i - nbins_3sigma); j < min(nbins - 1, i + nbins_3sigma); j++) {
			// Calculate spatial and intensity distances
			double spatial_dist = abs(i - j) * bin_size;
			double intensity_dist = abs(artmp[i] - artmp[j]);

			// weights
			double w_x = exp(-0.5 * (spatial_dist * spatial_dist) / (sigma_x * sigma_x));
			double w_y = exp(-0.5 * (intensity_dist * intensity_dist) / (sigma_y * sigma_y));
			double weight = w_x * w_y;

			weighted_sum += artmp[j] * weight;
			sum_of_weights += weight;
		}

		if (sum_of_weights != 0.) ar[i] = weighted_sum / sum_of_weights;
	}
	delete[] artmp;
}

/// @brief Nonlinear filter based on bilateral filter
/// 
/// Uses data and slope of data to weight the smoothing. x is not used. Very slow.
/// 
/// @param ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param sigma_x Window along x in bins.
/// @param sigma_y Sigma along y in mV.
/// @param sigma_slope Sigma of slope in mV.
/// @param bin_size Bin size in ns.
void Filters::Bilateral2Filter(double*& ar, int nbins, int sigma_x, double sigma_y, double sigma_slope, double bin_size) {
	// Asymmetric smoothing using a gauss of the change of the data relative to the central bin in a window for weighting
	double* artmp = new double[nbins]; 
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	double denom = -2. * sigma_y * sigma_y;
	int x_shift = static_cast<int>(sigma_x / bin_size);

	double weightedSum = 0.;
	double sumOfWeights = 0.;

	double slope[nbins - 1];
	double slope_weight[nbins - 1];
	for (int i = 0; i < nbins - 1; i++) {
		slope[i] = ar[i + 1] - ar[i];
		slope_weight[i] = 1. / (1. + pow(.5 * slope[i] / sigma_slope, 2.));
	}

	for (int i = 0; i < nbins; i++) {
		weightedSum = 0.;
		sumOfWeights = 0.;
		int start =	max(0, i - x_shift);
		int end =	min(i + x_shift, nbins - 1);

		for (int j = start; j <= end; j++) {
			double diff = artmp[i] - artmp[j];
			double weight = exp(pow(diff, 2) / denom);
			if (j < nbins - 1) weight *= slope_weight[j];
			weightedSum += weight * artmp[j];
			sumOfWeights += weight;
		}
		if (sumOfWeights != 0) ar[i] = weightedSum / sumOfWeights;
	}
	delete[] artmp;
}

/// @brief Median filter
/// @param ar Array to be smoothed.
/// @param nbins Number of bins of input.
/// @param window_size Window size for median in bins.
void Filters::MedianFilter(double*& ar, int nbins, int window_size) {
	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	int half_window = window_size / 2;

	for (int i = 0; i < nbins; i++) {
		double* window = new double[window_size];
		int window_index = 0;

		for (int j = i - half_window; j <= i + half_window; j++) {
			if (j >= 0 && j < nbins) {
				window[window_index] = artmp[j];
				window_index++;
			}
		}

		// calculate the median
		sort(window, window + window_index);
		double median = (window_index % 2 == 0) ? 
			(window[window_index / 2 - 1] + window[window_index / 2]) / 2. :
			window[window_index / 2];

		ar[i] = median;

		delete[] window;
	}
	delete[] artmp;
}

/// @brief Custom filter emulating primitive response function
/// 
/// Can be used to highlight peaks and suppress long tails (see pole-zero cancellation). 
/// Emulates response function of an amplifier setup. \n
/// Uses constant padding to avoid edge effects. \n
/// Use Filter_test.ipynb to test parameters.
/// 
/// @param[in,out] ar Array to be filtered.
/// @param nbins Number of bins of input.
/// @param sigma1 First.
/// @param sigma2 Second.
/// @param factor Factor for negative part (<=1).
/// @param bin_size Bin width. Default is .3125. Set to 1 to get sigma in units of bins.
void Filters::ResponseFilter(double*& ar, int nbins, double sigma1, double sigma2, double factor, double bin_size) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	// shifted difference of two gauss functions (~smoothed derivative)
	int nbins_23sigma = static_cast<int>(ceil((2. * sigma1 + 3. * sigma2) / bin_size));
	int nbins_3sigma = static_cast<int>(ceil(3. * sigma2 / bin_size));
	double sdog[nbins_23sigma];

	double denom1 = 2. * sigma1 * sigma1;
	double denom2 = 2. * sigma2 * sigma2;
	double norm = 0.;
	for (int i = 0; i < nbins_23sigma; i++) {
		sdog[i] = exp(-1. * pow(static_cast<double>(i) * bin_size - 3. * sigma2, 2.) / denom1) -
			factor * exp(-1. * pow(static_cast<double>(i) * bin_size - 2. * sigma2, 2.) / denom2);
		if (sdog[i] > 0.) norm += sdog[i];
	}
	if (norm == 0.) norm = 1.;

	for (int i = 0; i < nbins; i++) {
		double res = 0.;
		for (int j = -1 * nbins_3sigma; j < nbins_23sigma - nbins_3sigma; j++) {
			if (i + j >= 0 && i + j < nbins) res += sdog[j + nbins_3sigma] * artmp[i + j];
			else if (i + j < 0) res += sdog[j + nbins_3sigma] * artmp[0];
			else res += sdog[j + nbins_3sigma] * artmp[nbins - 1];
		}
		ar[i] = res / norm;
	}
	delete[] artmp;
}

/// @brief Shifted second order underdamped filter
/// 
/// Emulating 2nd order underdamped response (simple damped harmonic oscillator). \n
/// Will only consider 3 periods of the response function. \n
/// In development, for testing only.
/// 
/// @param[in,out] ar Array to be filtered.
/// @param nbins Number of bins of input.
/// @param period Period of response function in ns.
/// @param damping Damping time constant of the response function in ns. Needs to be > 0.
/// @param shift Shift of response function in ns.
/// @param bin_size Bin width of the input.
/// @param do_plot Plot response function and save to file.
void Filters::SecondOrderUnderdampedFilter(double*& ar, int nbins, double period, double damping, double shift, double bin_size, bool do_plot) {

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	int shift_n = static_cast<int>(ceil(shift / bin_size));
	int nbins_response = shift_n + static_cast<int>(3. * period / bin_size);
	
	if (nbins_response > nbins) {
		cout << "ERROR: nbins_response > nbins. Reduce shift or period." << endl;
		nbins_response = nbins - 1;
	}
	if (damping <= 0) {
		cout << "ERROR: damping <= 0. Setting to 1 ns" << endl;
		damping = 1;
	}

	if (nbins_response % 2 == 0) nbins_response++;
	double souf[nbins_response];

	double exp_factor = -1. / damping;
	double omega = 2. * M_PI / period;
	double sum_norm = 0.;
	for (int i = 0; i < nbins_response; i++) {
		double position = static_cast<double>(i) * bin_size;
		
		if (i < shift_n) souf[i] = 0.;
		else {
			souf[i] = exp(exp_factor * (position - shift)) * sin(omega * (position - shift));
			sum_norm += souf[i];
		}
	}
	if (sum_norm == 0.) sum_norm = 1.;

	for (int i = 0; i < nbins; i++) {
		double res = 0.;
		for (int j = max(0, nbins_response / 2 - i); j < min(nbins_response, nbins + nbins_response / 2 - i); j++) {
			res += souf[nbins_response - j - 1] * artmp[i + j - nbins_response / 2];
		}
		ar[i] = res / sum_norm;
	}

	// plot the results
	if (do_plot) {
		double x[nbins_response];
		double x_full[nbins];
		for (int i = 0; i < nbins; i++) {
			x_full[i] = static_cast<double>(i) * bin_size;
			if (i < nbins_response) x[i] = static_cast<double>(i) * bin_size;
		}

		TCanvas* c_response = new TCanvas("response", "response", 800, 600);
		TGraph* g_original = new TGraph(nbins, x_full, artmp);
		g_original->SetTitle("original");
		g_original->SetMarkerStyle(34);
		g_original->Draw("AP");

		TGraph* g_response = new TGraph(nbins_response, x, souf);
		g_response->SetTitle("response");
		g_response->SetLineColor(kBlue);
		g_response->SetLineWidth(5); 
		g_response->SetLineStyle(7);
		g_response->Draw("L same");

		TGraph* g_result = new TGraph(nbins, x_full, ar);
		g_result->SetTitle("result");
		g_result->SetLineColor(kRed);
		g_result->SetLineWidth(3);
		g_result->Draw("L same");

		g_response->Scale(TMath::MaxElement(nbins, ar) / TMath::MaxElement(nbins_response, souf));
		c_response->Update(); c_response->Modified();
		gPad->BuildLegend(.5, .7, .9, .9);
		c_response->SaveAs("response.png");
	}
	delete[] artmp;
}

/// @brief Shifted second order underdamped filter
/// See above for details.
/// @param[in,out] his Histogram to be filtered.
void Filters::SecondOrderUnderdampedFilter(TH1F*& his, int nbins, double period, double damping, double shift, double bin_size, bool do_plot) {
	double* yvals = Helpers::gety(his);
	SecondOrderUnderdampedFilter(yvals, nbins, period, damping, shift, bin_size, do_plot);
	for (int i = 1; i <= his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i - 1]);
	delete[] yvals;
}