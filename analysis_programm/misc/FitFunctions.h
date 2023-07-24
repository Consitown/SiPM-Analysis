/// @brief Default fit function for SiPMs missing after-pulses and dark counts
/// 
/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
/// See https://root.cern/manual/fitting/ for ROOT fitting.
/// 
class Fitf {
public:
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// @return Function value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		double mu = p[1];
		double lambda = p[2];

		double sigma0 = p[3];
		double sigma1 = p[4];

		double G = p[5];
		double B = p[6];

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2.));
		}
		return sum;
	};
};

/// @brief Default fit function for SiPMs with after-pulses missing dark counts
/// 
/// Still missing dark counts in integration window (3.3 in paper). \n \n
/// 
/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
/// See https://root.cern/manual/fitting/ for ROOT fitting.
/// 
class Fitf_full {
public:
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// 7 - alpha: after-pulsing probability \n 
	/// 8 - beta: the inverse of the exponential slope of the after-pulse PH distribution
	/// @return Function value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;
		
		double mu = p[1];
		double lambda = p[2];

		double sigma0 = p[3];
		double sigma1 = p[4];

		double G = p[5];
		double B = p[6];

		double alpha = p[7];
		double beta = p[8];

		for (int kint = 0; kint <= kmax; kint++) {
			//numbe of fired cells
			double k = static_cast<double>(kint);
			//pulse width
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			double gausnormsigmak = 1. / (sqrt(2. * TMath::Pi()) * sigmaK);
			//gauss peak
			double gauss = gausnormsigmak * TMath::Exp(-1. * TMath::Power(x[0] - (k * G + B), 2.) / (sigmaK * sigmaK * 2.));

			//generalized poisson envelope
			double gp = p[0] * mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += (gp * gauss);
			}
			else {
				//after-pulse + delayed cross talk
				double bk0 = TMath::Power(1. - alpha, k);
				double bk1 = TMath::Factorial(kint) / TMath::Factorial(kint - 1) * alpha * TMath::Power(1. - alpha, k - 1);
				double pk1 = TMath::Exp(-1. * (x[0] - (k * G + B)) / beta) * gausnormsigmak / beta * sigmaK * sqrt(TMath::Pi() / 2) * (TMath::Erf((x[0] - (k * G + B)) / (sqrt(2) * sigmaK)) + 1.);


				if (kint == 1) sum += gp * (bk0 * gauss + bk1 * pk1);
				else {
					double api2k = 0.;
					for (int ii = 2; ii < kint; ii++) {
						double iid = static_cast<double>(ii);
						double bkialpha = TMath::Factorial(kint) / (TMath::Factorial(ii) * TMath::Factorial(kint - ii)) * TMath::Power(alpha, iid) * TMath::Power((1. - alpha), (k - iid));
						double dpkidph = 0;
						if (x[0] > (k * G + B)) dpkidph = TMath::Power((x[0] - (k * G + B)), (iid - 1.)) / (TMath::Factorial(ii - 1) * TMath::Power(beta, iid)) * TMath::Exp(-1. * (x[0] - (k * G + B)) / beta);
						api2k += bkialpha * dpkidph;
					}
					sum += gp * (bk0 * gauss + bk1 * pk1 + api2k);
				}
			}
		}
		return sum;
	};
};

/// @brief Fit function for SiPMs missing after-pulses and dark counts but with biased pedestal
/// 
/// See https://arxiv.org/abs/1609.01181 for explanation of fit function. \n 
/// See https://root.cern/manual/fitting/ for ROOT fitting.
/// 
class Fitf_biased {
public:
	/// @param x 
	/// @param p 
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Virtual pedestal shift of pe peaks \n 
	/// 7 - Pedestal scaling for biased pedestal \n 
	/// 8 - Position of biased pedestal
	/// @return Function value 
	/// @return Func value
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		double mu = p[1];
		double lambda = p[2];

		double sigma0 = p[3];
		double sigma1 = p[4];

		double G = p[5];
		double B = p[6];

		double a_ped = p[7];
		double x_ped = p[8];

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += a_ped * p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - x_ped) / sqrt(2) / sigmaK), 2));
			}
			else {
				sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
			}
		}
		return sum;
	};
};

/// @brief Gauss-Poisson distribution for fit of PMT charge spectra
/// 
/// See https://doi.org/10.1016/0168-9002(94)90183-X 
///
class Fitf_PMT {
public: 
	/// @param x 
	/// @param p 
	/// 0 - A:		normalization to number of events in fit region \n 
	/// 1 - w:		probability for type II BG \n 
	/// 2 - alpha:	coefficient of exponential decrease of typ II BG \n 
	/// 3 - sigma0:	sigma of pedestal \n 
	/// 4 - Q0:		position of pedestal \n 
	/// 5 - mu:		mean number of PE \n 
	/// 6 - sigma1:	width of 1 PE peak \n 
	/// 7 - Q1:		position of 1 PE peak
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;
		int kmax = static_cast<int>(ceil(p[5])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else if (x[0] != p[4]) gn_term *= 0.; // delta function

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

/// @brief Gauss-Poisson distribution for fit of PMT charge spectra
/// 
/// With biased pedestal peak. \n
/// See https://doi.org/10.1016/0168-9002(94)90183-X 
/// 
class Fitf_PMT_pedestal {
public:
	/// @param x 
	/// @param p
	/// 0 - A:		normalization to number of events in fit region \n 
	/// 1 - w:		probability for type II BG \n 
	/// 2 - alpha:	coefficient of exponential decrease of typ II BG \n 
	/// 3 - sigma0:	sigma of pedestal \n 
	/// 4 - Q0:		position of pedestal \n 
	/// 5 - mu:		mean number of PE \n 
	/// 6 - sigma1:	width of 1 PE peak \n 
	/// 7 - Q1:		position of 1 PE peak \n 
	/// 8 - norm0:	norm of 0 PE peak 
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;
		int kmax = static_cast<int>(ceil(p[5])) * 10;

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[5], k) * TMath::Exp(-p[5]) / TMath::Factorial(kint);

			double normgn_term = 1.;
			if (kint != 0) normgn_term /= (p[6] * TMath::Sqrt(2. * TMath::Pi() * k));
			else normgn_term *= p[8] / (p[3] * TMath::Sqrt(2. * TMath::Pi()));
			double gn_term = (1. - p[1]) * normgn_term;
			if (kint != 0) gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - k * p[7] - p[4], 2.) / (2. * k * p[6] * p[6]));
			else  gn_term *= TMath::Exp(-1. * TMath::Power(x[0] - p[4], 2.) / (2. * p[3] * p[3]));

			double Qn = p[4] + k * p[7];
			double sigman = TMath::Sqrt(p[3] * p[3] + k * p[6] * p[6]);
			double ignxe_term_exp = p[2] / 2. * TMath::Exp(-1. * p[2] * (x[0] - Qn - p[2] * sigman * sigman));

			double sgn_arg = x[0] - Qn - sigman * sigman * p[2];
			double arg_sgn = 1.;
			if (sgn_arg < 0.) arg_sgn = -1.;
			double ignxe_term_erf = TMath::Erf(fabs(p[4] - Qn - sigman * sigman * p[2]) / (sigman * TMath::Sqrt(2.))) + arg_sgn * TMath::Erf(fabs(sgn_arg) / (sigman * TMath::Sqrt(2.)));

			pmt_charge_spectrum += p[0] * poiss * (gn_term + p[1] * ignxe_term_exp * ignxe_term_erf);
		}
		if (pmt_charge_spectrum < 0.) pmt_charge_spectrum = 0.;
		return pmt_charge_spectrum;
	};
};

/// @brief Ideal PMT charge spectrum
/// 
/// Gives very good fit but will not describe pedestal well. \n 
/// See https://doi.org/10.1016/0168-9002(94)90183-X 
///
class Fitf_PMT_ideal {
public: 
	/// @param x 
	/// @param p 
	/// 0 - A_s:		norm. of PE spectrum \n 
	/// 1 - mu:		mean number of PE \n 
	/// 2 - sigma1:	width of 1st PE peak \n 
	/// 3 - Q1:		position (gain*e) of 1st peak
	/// @return Func value
	double operator() (double* x, double* p) {
		double pmt_charge_spectrum = 0.;

		int kmax = static_cast<int>(ceil(p[1])) * 10;

		for (int kint = 1; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);

			double poiss = TMath::Power(p[1], k) * TMath::Exp(-1. * p[1]) / TMath::Factorial(kint);

			double norm = 1. / (p[2] * TMath::Sqrt(2. * TMath::Pi() * k));

			double gauss = TMath::Exp(-1. * TMath::Power(x[0] - k * p[3], 2) / (2. * k * p[2] * p[2]));

			pmt_charge_spectrum += p[0] * poiss * norm * gauss;
		}
		return pmt_charge_spectrum;
	};
};

/// @brief Landau-Gauss-convolution
/// 
/// Used for large photon yields. \n 
/// From https://root.cern.ch/doc/master/langaus_8C.html => \n
/// In the Landau distribution (represented by the CERNLIB approximation), \n 
/// the maximum is located at x=-0.22278298 with the location parameter=0. \n 
/// This shift is corrected within this function, so that the actual \n 
/// maximum is identical to the MP parameter. 
/// 
class Fitf_langaus {
public:
	/// @param x 
	/// @param par
	/// par[0]=Width (scale) parameter of Landau density \n 
	/// par[1]=Most Probable (MP, location) parameter of Landau density \n 
	/// par[2]=Total area (integral -inf to inf, normalization constant) \n 
	/// par[3]=Width (sigma) of convoluted Gaussian function
	/// @return Func value
	double operator() (double* x, double* par) {
		// Numeric constants
		double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
		double mpshift = -0.22278298;       // Landau maximum location

		// Control constants
		double np = 100.0;      // number of convolution steps
		double sc = 5.0;      // convolution extends to +-sc Gaussian sigmas

		// Variables
		double xx;
		double mpc;
		double fland;
		double sum = 0.0;
		double xlow, xupp;
		double step;
		double i;


		// MP shift correction
		mpc = par[1] - mpshift * par[0];

		// Range of convolution integral
		xlow = x[0] - sc * par[3];
		xupp = x[0] + sc * par[3];

		step = (xupp - xlow) / np;

		// Convolution integral of Landau and Gaussian by sum
		for (i = 1.0; i <= np / 2; i++) {
			xx = xlow + (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);

			xx = xupp - (i - .5) * step;
			fland = TMath::Landau(xx, mpc, par[0]) / par[0];
			sum += fland * TMath::Gaus(x[0], xx, par[3]);
		}

		return (par[2] * step * sum * invsq2pi / par[3]);
	};
};

/// @brief Fit function for SiPMs missing after-pulses and dark counts but including additional dark count spectrum
/// 
/// Sum of two spectra for event spectrum + dark count (background trigger) spectrum. 
/// To be used if the trigger and filter selections are not perfect and there are still empty events.
/// 
class Fitf_plus_DC {
public:
	/// @param x 
	/// @param p
	/// 0 - N0: Normalization (~Number of events) \n 
	/// 1 - mu: for generalized poisson distribution \n 
	/// 2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda) \n 
	/// 3,4 -sigma0, sigma1 \n 
	/// 5 - G: gain \n 
	/// 6 - B: Pedestal \n 
	/// 7 - mu_dc: mu for dark count rate in dark events (dark events = noise/background triggers without photons in SiPMs) \n 
	/// 8 - N0_dc: fraction (dark events)/(total number of events) 
	/// @return 
	double operator() (double* x, double* p) {
		double sum = 0;
		int kmax = static_cast<int>(ceil(p[1])) * 10;

		double mu = p[1];
		double lambda = p[2];

		double sigma0 = p[3];
		double sigma1 = p[4];

		double G = p[5];
		double B = p[6];

		double mu_dc = p[7];
		double N0_dc = p[8];

		for (int kint = 0; kint <= kmax; kint++) {
			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = ((1. - N0_dc) * (mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda))) + N0_dc * (mu_dc * TMath::Power((mu_dc + k * lambda), k - 1) * TMath::Exp(-(mu_dc + k * lambda)))) / TMath::Factorial(kint);

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};

/// @brief Periodic gauss fit function for angular reconstruction of SiPM array signal (SHiP SBT WOM tube)
class Fitf_periodic_gauss {
public:
	/// @param x 
	/// @param p
	/// 0 - A: Normalization \n 
	/// 1 - mu: Gauss pos. \n 
	/// 2 - sigma: Gauss sigma \n 
	/// 3 - B: Constant offset \n 
	/// @return 
	double operator() (double* x, double* p) {
		double k = 0;
		int k_pm = 2; // how many Gaussians to consider -> k_pm=2 means 2+1+2=5 Gaussians
		double A = p[0];
		double mu = p[1];
		double sigma = p[2];
		double B = p[3];

		double gauss = 0.;

		for (int kint = -k_pm; kint <= k_pm; kint++) {
			k = static_cast<double>(kint);
			// periodic gauss
			gauss += A*TMath::Exp(-TMath::Power(x[0] - mu + k * 360., 2.)/(2.*sigma*sigma));
		}
		return gauss + B;
	};
};

/// @brief Integral of Landau-Gauss-Convolution ("S-curve")
class Fitf_langaus_int {
public:
	/// @param x
	/// @param par
	/// par[0]=Width (scale) parameter of Landau density \n 
	/// par[1]=Most Probable (MP, location) parameter of Landau density \n 
	/// par[2]=Total area (integral -inf to inf, normalization constant) \n 
	/// par[3]=Width (sigma) of convoluted Gaussian function
	/// par[4]=Maximum value of integral (usually 1)
	/// @return Func value
	double operator() (double* x, double* par) {
		// Control constants
		Double_t sc = 5.;        // integral extends to MPV-sc*(Gaussian sigma + Landau width)

		// Variables
		Double_t start = par[1] - sc * (par[0] + par[3]);

		Fitf_langaus func;
		TF1* fun = new TF1("fun", func, start - 1, x[0] + 1, 4);
		fun->SetParameter(0, par[0]);
		fun->SetParameter(1, par[1]);
		fun->SetParameter(2, par[2]);
		fun->SetParameter(3, par[3]);

		return (par[4] * (1. - fun->Integral(start, x[0])));
	};
};