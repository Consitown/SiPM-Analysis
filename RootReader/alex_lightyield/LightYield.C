//root
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TFitResult.h>
#include <TLine.h>
#include <TMatrixD.h>

//C, C++
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

using namespace std;






void strip_ext(char *fname)
{
	char *end = fname + strlen(fname);

	while (end > fname && *end != '.')
	{
		--end;
	}

	if (end > fname)
	{
		*end = '\0';
	}
}

vector<string> split(const string &str, const string &delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos)
			pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty())
			tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}
void extractIntegerWords(string str) 
{ 
    stringstream ss;     
  
    /* Storing the whole string into string stream */
    ss << str; 
  
    /* Running loop till the end of the stream */
    string temp; 
    int found; 
    while (!ss.eof()) { 
  
        /* extracting word by word from stream */
        ss >> temp; 
  
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found) 
            cout << found << " "; 
  
        /* To save from space at the end of string */
        temp = ""; 
    } 
} 

/*******************
__ FUNCTIONS ______
*******************/
Double_t poisson(Double_t *x, Double_t *par) //Poisson called with par[0] = normalization, par[1] = lamda
{
	Double_t fitval = par[0] * TMath::Power(par[1], *x) * TMath::Exp(-par[1]) / TMath::Gamma(*x + 1);

	return fitval;
}


Double_t generalized_poisson(Double_t *x, Double_t *par) //Poisson called with par[0] = normalization, par[1] = theta, par[2] = lambda
{
	Double_t fitval = par[0] * par[1] * TMath::Power((par[1] + *x * par[2]), (*x - 1)) * TMath::Exp(-par[1] - *x * par[2]) / TMath::Gamma(*x + 1);

	return fitval;
}

Double_t langaufun(Double_t *x, Double_t *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) 
	  {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
 
   Int_t i;
   Char_t FunName[100];
 
   sprintf(FunName,"Fitfcn_%s",his->GetName());
 
   //TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   //if (ffitold) delete ffitold;
 
   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
 
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }
 
   his->Fit(FunName,"RB");   // fit within specified range, use ParLimits
 
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf
 
   return (ffit);              // return fit function
 
}


 
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &maxx_err, Double_t &FWHM) {
 
   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.
 
   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;
 
 
   // Search for maximum
 
   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-1);
 
   maxx = x;
   maxx_err = step; // the error of this algorithm is the step width
 
   fy = l/2;
 
 
   // Search for right x location of fy
 
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-2);
 
   fxr = x;
 
 
   // Search for left x location of fy
 
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-3);
 
 
   fxl = x;
 
   FWHM = fxr - fxl;
   return (0);
}
 
 

/*******************
__ MAIN_____________
*******************/

int main(int argc, char *argv[]){

	cout<< "VALUE: "<<argv[1]<<endl;

	int fit_function = 5; //1 for Poisson, 2 for generalized Poisson, 3 for Gaussian, 4 for Landau, 5 for Landau*Gaussian

	string fileLocation=string(argv[1]);
	vector<string> splitted=split(string(argv[1]),"//");
	string run_name=splitted[1];
	vector<string> runInformations=split(run_name.c_str(),"_"); 



	int Xmin = -10, Xmax = 400; // histogram range
	double range_lo = 3, range_hi = 140; // landau fit range
	float frac_l=0.5,frac_u=0.5; // gauss fit range, threshold fraction
	int n_bins = 1000;
	double bin_error = (Xmax - Xmin) / n_bins;
	
	
	int n_ch = 8 ; // number of SiPM channels
	TH1F * h_vec[n_ch]; // histograms
	TLine * ln_vec[n_ch]; // vertical lines
	TF1 * fit_vec_g[n_ch], * fit_vec_l[n_ch], * fit_vec_v[n_ch]; // gauss fit functions
	TF1 * fit_vec_p[n_ch]; // poisson fit functions
	double amp_max[n_ch], peak[n_ch], peak_err[n_ch], rchi2_g[n_ch]; // gauss fit results
	double amp_max_p[n_ch], expect_p[n_ch], expect_err_p[n_ch], rchi2_p[n_ch], max_p[n_ch], max_err_p[n_ch]; // poisson fit results
	int h_entries[n_ch];

	TF1 *fitsnr; //fits for Landau*Gauss
	Double_t chisqr_lg;
	Int_t    ndf_lg; //Landau*Gauss fit results
	Double_t SNRPeak, SNRPeak_err, SNRFWHM, MPV_langau, MPV_langau_err;

	/***** 
	__ Lines for distribution and fit things _______
	*****/
	TLine * poisson_max;
	TLine * langaus_max;
	TLine * langaus_mpv;
	TLine * poisson_exp;

	 
	/***** 
	__ for the sum of all histograms _______
	*****/
	TH1F * h_sum;
	h_sum = new TH1F("Sum","Sum",n_bins,Xmin-0.5,Xmax-0.5);

	TCanvas *C1;
	C1 = new TCanvas("ly_dist","",1200,750);
	C1->Divide(4,3);

	/***** 
	__ File export _______
	*****/

	// for channel export
	ofstream channel_save_file_poisson ("Channel_values_poisson.txt", ios::app);
	ofstream channel_save_file_langau ("Channel_values_langau.txt", ios::app);
	ofstream channel_save_file_mean ("Channel_values_mean.txt", ios::app);
	ofstream channel_save_file_trunc_mean ("Channel_values_trunc_mean.txt", ios::app);



	//26_pos7_angle0_e52_ch32.root
	//1_calib_vb58_sipm5d


	// parse input arguments
	//int run_nr = atoi(runInformations[1]);
	/*
	if (s1.find("pos") != std::string::npos) {
    std::cout << "found!" << '\n';
	}
	int pos = atoi(extractIntegerWords[runInformations[2]]);
	string WC_v = (string)argv[5];
	//string particle = (string)argv[5];
	double energy = (double)atoi(extractIntegerWords[runInformations[2]])/10;
	//string target_wom = (string)argv[7];
*/

	/***** 
	__ INITIALIZE ___________________________
	*****/

	// ROOT GUI
	bool is_interactive = 0;
	// amp (0), charge (1)
	bool is_charge = 1;	
	// store results
	bool store_result = 1;
	bool legend=true;

	bool store_channel = true; //to store results from the single channels


	TApplication * ROOTapp;
	if (is_interactive){ROOTapp = new TApplication("ROOT interactive", &argc, argv);}
	

	// style settings
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1); // export high resolution .pdf
	float alpha = 0.1;
	int lineW = 2;

	

	// open tree
	// string root_dir = "..//runs/"+run_name+"/"+run_name+".root";
	string out_dir = "./ly_histograms/charge/"; // dir. to export plots and txt file
	// export file names
	string out_pdf_name = Form("light_yield_%s.pdf",run_name.c_str()) ;// histogram pdf filename
	//string out_values_list_name = Form("mpv_%s.txt",target_wom.c_str()) ; // MPV list filename
	//string out_values_err_list_name = Form("mpv_err_%s.txt",target_wom.c_str()) ; // MPV err list


	TFile* file = new TFile(fileLocation.c_str());
	TTree* tree;
	file->GetObject("T",tree);


	// default hist x-axis maximum
	int upperA = 20000, upperB = 100, upperC = 450, upperD = 450;



	//save the run information
	if (store_channel && fit_function == 1) channel_save_file_poisson << runInformations[0] << "\t";
	if (store_channel && fit_function == 5) channel_save_file_langau << runInformations[0] << "\t";
	if (store_channel) channel_save_file_mean << runInformations[0] << "\t";
	if (store_channel) channel_save_file_trunc_mean << runInformations[0] << "\t";

	//___ LOOP OVER CHANNELS ___
	for (int i = 0; i < n_ch; ++i)
	{
		cout << "\n\nIch bin jetzt in Loop " << i << endl;
		
		TString h_name1, h_title1, draw_cmnd1, cut_cmnd1;

		// set current wom id
		string wom_id;
		if (i<8)
		{wom_id = "WOM-D";}
		else if (i>=8 && i<16)
		{wom_id = "WOM-C";}
		else if (i>=16 && i<24)
		{wom_id = "WOM-A";}
		else if (i>=24 && i<31)
		{wom_id = "WOM-B";}


		/***** 
		__ DRAW HISTOGRAM ______________________________
		*****/


		C1->cd(i+1);
		h_name1.Form("h%d",i);
		h_title1.Form("%s, ch%d",run_name.c_str(),i);

		draw_cmnd1.Form("Integral[%d]>>h%d",i,i);
		cut_cmnd1.Form("");

		gPad->SetRightMargin(0.00);
		gPad->SetLeftMargin(.12);
		

		h_vec[i] = new TH1F(h_name1, h_title1, n_bins, Xmin-0.5, Xmax-0.5);

		h_vec[i]->GetXaxis()->SetTitle("Integrated charge [mV #times ns]");
		// h_vec[i]->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}]");
		//h_vec[i]->GetXaxis()->SetTitleOffset(1.3);
		h_vec[i]->GetYaxis()->SetTitle("Entries");
		h_vec[i]->GetXaxis()->SetLabelSize(0.05);
		h_vec[i]->GetYaxis()->SetLabelSize(0.05);
		h_vec[i]->GetXaxis()->SetTitleSize(0.05);
		h_vec[i]->GetYaxis()->SetTitleSize(0.05);
		h_vec[i]->SetLineColorAlpha(kBlack,0.7);
		h_vec[i]->SetFillColorAlpha(kBlack,alpha);
		h_vec[i]->SetMarkerStyle(8);
		h_vec[i]->SetMarkerSize(0.2);
		h_vec[i]->SetMarkerColorAlpha(kBlack,0.6);
		
		cout << "\n\n" << endl;
		/*cout << "draw_cmnd1 ist " << draw_cmnd1 << endl;
		cout << "cut_cmnd1 ist " << cut_cmnd1 << endl;*/

		tree->Draw(draw_cmnd1,cut_cmnd1,"HISTE");


		h_sum->Add(h_vec[i]);


		
	
		/***** 
		__ MEAN_________________
		*****/
		double h_max_y = h_vec[i]->GetMaximum();
		double h_mean =  h_vec[i]->GetMean();
		double h_mean_err = h_vec[i]->GetStdDev();
		h_mean_err /= sqrt(h_vec[i]->GetEntries());
		

		/***** 
		__ TRUNCATED MEAN_________________
		*****/
		h_vec[i]->GetXaxis()->SetRangeUser(0, 150); // set the axis range to 0 to 150
		double h_trunc_mean =  h_vec[i]->GetMean();
		double h_trunc_mean_err = h_vec[i]->GetStdDev();
		h_trunc_mean_err /= sqrt(h_vec[i]->GetEntries()-100);
		h_vec[i]->GetXaxis()->SetRange(0,0); //unsetting the axis range
    	
		/***** 
		__ MAXIMUM_________________
		*****/
    	// get histogram maximum value & position in range (buggy...)
		// h_vec[i]->GetXaxis()->SetRangeUser(0,100);
		
		amp_max[i] = h_vec[i]->GetXaxis()->GetBinCenter(h_vec[i]->GetMaximumBin());
		float max = h_vec[i]->GetMaximum();
		int lower_bin = h_vec[i]->GetMaximumBin();
		int upper_bin = h_vec[i]->GetMaximumBin();
		// h_vec[i]->GetXaxis()->UnZoom();
		// printf("ch%d max amp: %f\n",i,amp_max[i]);		

		ln_vec[i] = new TLine(amp_max[i],0.1,amp_max[i],h_max_y);
		ln_vec[i]->SetLineColor(8);
		ln_vec[i]->SetLineWidth(1);
		ln_vec[i]->SetLineStyle(5);
		ln_vec[i]->Draw("same");


		/***** 
		__ MEDIAN ______________________________
		*****/	
		Double_t median_x, q;
		q = 0.5;
		h_vec[i]->GetQuantiles(1,&median_x, &q);
		cout << "Der Median wurde zu " << median_x << " berechnet." << endl;
		TLine * ln_median = new TLine(median_x, 0, median_x ,h_max_y);
		ln_median->SetLineColor(9);
		ln_median->SetLineStyle(9);
		ln_median->Draw("same");


		/***** 
		__ FIT PEAK: GAUSS ______________________________
		*****/
		/*printf("max %1.1f\n",max);
		printf("h_max_y %1.1f\n",h_max_y );

		// fit range ± 0.5x maximum
		// while (h_vec[i]->GetBinContent(lower_bin) > max*frac_l && (h_vec[i]->GetBinLowEdge(lower_bin))>0.5 ) {lower_bin--;}
		while (h_vec[i]->GetBinContent(lower_bin) > max*frac_l) {lower_bin--;}
		while (h_vec[i]->GetBinContent(upper_bin) > max*frac_u) {upper_bin++;}

		float r1,r2;
		r1 = h_vec[i]->GetBinCenter(lower_bin);
		r2 = h_vec[i]->GetBinCenter(upper_bin);
		printf("%1.1f %1.1f\n",r1,r2 );

		// fit_vec_g[i] = new TF1("g_fit","gaus",0.5,30);
		fit_vec_g[i] = new TF1("g_fit","gaus",r1,r2);
		fit_vec_g[i]->SetParameter(1,10);
		fit_vec_g[i]->SetLineColor(4);
		fit_vec_g[i]->SetLineStyle(1);
		// fit_vec_g[i]->SetLineStyle(5);
		fit_vec_g[i]->SetNpx(1000); // draw function with high resolution 
		h_vec[i]->Fit("g_fit","RQM");
		h_vec[i]->Draw("same");

		//Mean of the Gaussian
		Double_t gauss_mean = fit_vec_g[i]->Mean(r1,r2);

		Double_t par_g[3]; // to store fit results
		fit_vec_g[i]->GetParameters(&par_g[0]);

		// results
		peak[i] = fit_vec_g[i]->GetParameter(1);
		peak_err[i] = fit_vec_g[i]->GetParError(1);

		double sigma_g = fit_vec_g[i]->GetParameter(2);
		double sigma_g_err = fit_vec_g[i]->GetParError(2);
		double rchi2_g = fit_vec_g[i]->GetChisquare()/fit_vec_g[i]->GetNDF();

		// draw line for Gaussian Maximum
		TLine * ln_gauss_mpv = new TLine(peak[i],0,peak[i],h_max_y);
		ln_gauss_mpv->SetLineColor(4);
		ln_gauss_mpv->SetLineStyle(3);
		ln_gauss_mpv->Draw("same");

		// draw line for Gaussian Mean
		/*TLine * ln_gauss_mean = new TLine(gauss_mean,0,gauss_mean,h_max_y);
		ln_gauss_mean->SetLineColor(2);
		ln_gauss_mean->SetLineStyle(4);
		ln_gauss_mean->Draw("same");*/



		/***** 
		__ FIT PEAK: LANDAU ______________________________
		*****/

		/*fit_vec_l[i] = new TF1("fit_l","landau",0.5,lan_range_hi);
		fit_vec_l[i]->SetLineColor(2);
		fit_vec_l[i]->SetNpx(1000); // draw function with high resolution 
		// fit_vec_l[i]->SetParameters(&par_g[0]);

		h_vec[i]->Fit("fit_l","RQMN");
		// h_vec[i]->Draw("sameFUNC");

		 // store fit results
		double mpshift  = -0.22278298;
		double lan_sigma = fit_vec_l[i]->GetParameter(2);
		double lan_sigma_err = fit_vec_l[i]->GetParError(2);
		float lan_mpv = fit_vec_l[i]->GetParameter(1);
		float lan_mpv_err = fit_vec_l[i]->GetParError(1);
		float lan_rchi2 = fit_vec_l[i]->GetChisquare()/fit_vec_l[i]->GetNDF();
		
		double lan_mpv_c = lan_mpv + mpshift*lan_sigma ;
		double lan_mpv_c_err = sqrt(lan_mpv_err*lan_mpv_err + (mpshift*lan_sigma_err)*(mpshift*lan_sigma_err)) ;
		
		

		// // draw line for MPV
		// TLine * ln_lan_mpv = new TLine(lan_mpv_c,0,lan_mpv_c,h_max_y);
		// ln_lan_mpv->SetLineColor(2);
		// ln_lan_mpv->SetLineStyle(2);
		// ln_lan_mpv->Draw("same");

		// // custom histogram legend with fit results
		// TLegend *h_lan_leg = new TLegend(0.53,0.30,1.0,0.9);
		// h_lan_leg->SetTextSize(0.04);
		// h_lan_leg->AddEntry(h_vec[i],Form("#bf{data}"),"elpf");
		// h_lan_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i]->GetEntries()),"");
		// h_lan_leg->AddEntry(ln_vec[i],Form("dist. max. = %1.2f ",amp_max[i]),"l");
		// h_lan_leg->AddEntry(fit_vec_l[i],Form("Landau fit"),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}_{L}/ndf = %1.1f",lan_rchi2),"");
		// h_lan_leg->AddEntry(ln_lan_mpv,Form("MPV_{L} = %1.2f #pm %1.2f ",lan_mpv_c,lan_mpv_c_err),"l");
		// // h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{L} = %1.2f #pm %1.2f ",lan_sigma,lan_sigma_err),"");
		// h_lan_leg->AddEntry(fit_vec_g[i],Form("Gauss fit"),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}_{G}/ndf = %1.1f",rchi2_g),"");
		// h_lan_leg->AddEntry(ln_gauss_mpv,Form("MPV_{G} = %1.2f #pm %1.2f ",peak[i],peak_err[i]),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{G} = %1.2f #pm %1.2f ",sigma_g,sigma_g_err),"");
		// h_lan_leg->Draw();



		/***** 
		__ FIT: LANDAU GAUSS ______________________________
		*****/
		if (fit_function == 5)
		{
			// set all the parameters for the fit
			Double_t fr[2] {range_lo, range_hi}; //fit range
			Double_t sv[4], pllo[4], plhi[4];
			Double_t fp[4], fpe[4]; // the final fit parameters and their errors are saved here

			pllo[0]=0.5; pllo[1]=2.0; pllo[2]=0.5; pllo[3]=0.4; // lower parameter limits
			plhi[0]=50.; plhi[1]=90.0; plhi[2]=1000000.0; plhi[3]=50.0; // upper parameter limits
			//   par[0]=Width (scale) parameter of Landau density
			//   par[1]=Most Probable (MP, location) parameter of Landau density
			//   par[2]=Total area (integral -inf to inf, normalization constant)
			//   par[3]=Width (sigma) of convoluted Gaussian function
			sv[0]=40; sv[1]=40.0; sv[2]=1000.0; sv[3]=25.; // starting values for the fit

			fitsnr = langaufit(h_vec[i], fr, sv, pllo, plhi, fp, fpe, &chisqr_lg, &ndf_lg);


			// to get the x- and y-values of the maximum
			langaupro(fp, SNRPeak, SNRPeak_err, SNRFWHM);
			if (SNRPeak_err < 0) SNRPeak_err *= -1;
			
			printf("Fitting done\nPlotting results...\n");

			cout << "Langaus Peak = " << SNRPeak << "#pm" << SNRPeak_err << endl;

			fitsnr->Draw("same");

			//draw line for maximum
			langaus_max = new TLine(SNRPeak,0,SNRPeak,h_max_y);
			langaus_max->SetLineColor(2);
			langaus_max->SetLineStyle(4);
			langaus_max->Draw("same");
		}

		/***** 
		__ FIT: POISSON ______________________________
		*****/

		if (fit_function == 1)
		{
			fit_vec_p[i] = new TF1("p_fit", poisson, range_lo, range_hi, 2); // 0 just normalization, 1 = #lambda = expected value and variance
			//fit_vec_p[i] = new TF1("p_fit","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", poisson_range_lo, poisson_range_hi);
			//fit_vec_p[i] = new TF1("p_fit", "[0] * TMath::Power([1],[2]) / TMath::Factorial([2]) * TMath::Exp(-[1])", 13, 130); // 0 just normalization, 1 = #lambda = expected value and variance
			//fit_vec_p[i] = new TF1("p_fit", "[0] * TMath::Poisson(x, [1])", poisson_range_lo, poisson_range_hi); // 0 just normalization, 1 = #lambda = expected value and variance
			cout << "h_mean = " << h_mean << endl;
			fit_vec_p[i]->SetParameters(1., h_mean);
			fit_vec_p[i]->SetLineColor(6);
			fit_vec_p[i]->SetLineStyle(1);
			fit_vec_p[i]->SetNpx(1000); // draw function with high resolution 
			h_vec[i]->Fit("p_fit","RQ");
			h_vec[i]->Draw("same");

			Double_t par_p[3]; // to store fit results
			fit_vec_p[i]->GetParameters(&par_p[0]);

			// results
			expect_p[i] = fit_vec_p[i]->GetParameter(1);
			expect_err_p[i] = fit_vec_p[i]->GetParError(1);

			double rchi2_p = fit_vec_p[i]->GetChisquare()/fit_vec_p[i]->GetNDF();

			// draw line for Poisson Maximum
			poisson_max = new TLine(expect_p[i],0,expect_p[i],h_max_y);
			poisson_max->SetLineColor(6);
			poisson_max->SetLineStyle(4);
			poisson_max->Draw("sameFUNC");
		}

		/***** 
		__ FIT: GENERALIZED POISSON ______________________________
		*****/
		if (fit_function == 2)
		{
			fit_vec_p[i] = new TF1("p_fit", generalized_poisson, range_lo, range_hi, 3); // 0 just normalization, 1 = theta, 2 = lambda
			fit_vec_p[i]->SetParameters(1000, h_mean, 1);
			fit_vec_p[i]->SetParNames("Norm.", "Theta", "Lambda");
			fit_vec_p[i]->SetParLimits(2, 0., 1.);
			fit_vec_p[i]->SetLineColor(6);
			fit_vec_p[i]->SetLineStyle(1);
			fit_vec_p[i]->SetNpx(1000); // draw function with high resolution 
			h_vec[i]->Fit("p_fit","RM+");
			h_vec[i]->Draw("same");

			/*if (i!=0)
			{
				TFitResultPtr r = h_vec[i]->Fit("p_fit", "S");
				TMatrixD cov = r->GetCovarianceMatrix();
				cov.Print();
			}*/
			Double_t par_p[3]; // to store fit results
			Double_t par_err_p[3];
			fit_vec_p[i]->GetParameters(&par_p[0]);
			for (int j = 0; j<3; j++) par_err_p[j] = fit_vec_p[i]->GetParError(j);


			/*for (int j=0; j<3; j++) {
				cout << "par_p[" << j << "] = " << par_p[j] << endl; // to print fit results
				cout << "par_error_p[" << j << "] = " << par_error_p[j] << endl; // to print fit errors
			}*/

			// results
			expect_p[i] = par_p[1] / (1 - par_p[2]); //formula for the expectation value of the generalized Poisson
			expect_err_p[i] = sqrt(1./(pow((1 - par_p[2]), 2)) * pow(par_err_p[1], 2) + par_p[1] / pow((1 - par_p[2]), 4) * pow(par_err_p[2], 2)); //formula from Gaussian error propagation

			cout << "expect_p[" << i << "] = " << expect_p[i] << endl;
			cout << "expect_err_p[" << i << "] = " << expect_err_p[i] << endl;
			cout << "Das Maximum des Poissonfits ist " << fit_vec_p[i]->GetMaximumX(range_lo + 1, range_hi - 1) << endl;

			double rchi2_p = fit_vec_p[i]->GetChisquare()/fit_vec_p[i]->GetNDF();

			// draw line for Poisson expectation value
			poisson_exp = new TLine(expect_p[i], 0, expect_p[i], h_max_y);
			poisson_exp->SetLineColor(6);
			poisson_exp->SetLineStyle(6);
			if (i!=0) poisson_exp->Draw("same");

			// maximum
			max_p[i] = fit_vec_p[i]->GetMaximumX(range_lo + 1, range_hi - 1);
			max_err_p[i] = (range_hi - range_lo) / 1000; // error of the max is one bin of the function

			// draw line for Poisson maximum value
			poisson_max = new TLine(max_p[i], 0, max_p[i], h_max_y);
			poisson_max->SetLineColor(6);
			poisson_max->SetLineStyle(6);
			poisson_max->Draw("same");
		}


		/***** 
		__  LEGEND	_______________________________________
		*****/
		TLegend *h_lan_leg = new TLegend(0.45,0.42,1.0,0.90);
		h_lan_leg->SetTextSize(0.045);
		h_lan_leg->AddEntry(h_vec[i],Form("#bf{data}"),"elpf");
		h_lan_leg->AddEntry((TObject*)0,Form("entries = %1.f", h_vec[i]->GetEntries()),"");		
		h_lan_leg->AddEntry(ln_vec[i],Form("Dist. max. = %1.1f #pm %1.1f",amp_max[i], bin_error),"l");
		h_lan_leg->AddEntry((TObject*)0,Form("Dist. mean = %1.2f #pm %1.2f", h_mean, h_mean_err),"");
		h_lan_leg->AddEntry((TObject*)0,Form("Dist. trunc. mean = %1.2f #pm %1.2f", h_trunc_mean, h_trunc_mean_err),"");
		h_lan_leg->AddEntry(ln_median,Form("Dist. median: %1.1f #pm %1.1f",median_x, bin_error),"l");
		if (fit_function == 1){ //for Poisson
			h_lan_leg->AddEntry(fit_vec_p[i],Form("Poisson fit"),"l");
			h_lan_leg->AddEntry((TObject*)0,Form("Poisson #chi^{2}/ndf = %1.1f",rchi2_p),"");
			h_lan_leg->AddEntry(poisson_max,Form("Poisson exp. value = %1.2f #pm %1.2f", expect_p[i], expect_err_p[i]),"l");
		}

		if (fit_function == 5){ //for Landau*Gauss
			h_lan_leg->AddEntry(fitsnr,Form("Landau*Gauss Fit"),"l");
			h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.2f", chisqr_lg/ndf_lg),"");
			h_lan_leg->AddEntry(langaus_max,Form("Fit Maximum = %1.2f #pm %1.2f", SNRPeak, SNRPeak_err), "l");
		}
		// h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{L} = %1.2f #pm %1.2f ",lan_sigma,lan_sigma_err),"");
		// h_lan_leg->AddEntry(fit_vec_g[i],Form("Gaussian fit"),"l");
		//h_lan_leg->AddEntry((TObject*)0,Form("Gauss #chi^{2}/ndf = %1.1f",rchi2_g),"");
		//h_lan_leg->AddEntry(ln_gauss_mpv,Form("Gaussian mean = %1.2f #pm %1.2f ",peak[i], peak_err[i]), "l");
		//h_lan_leg->AddEntry((TObject*)0,Form("#sigma = %1.2f #pm %1.2f",sigma_g,sigma_g_err),"");
		
		//h_lan_leg->AddEntry(poisson_max,Form("Poisson max = %1.2f #pm %1.2f", max_p[i], max_err_p[i]),"l");
		
		if(legend) h_lan_leg->Draw();


		// print results
		if (i<n_ch)
		{printf("%s, ch%d\n",run_name.c_str(),1 );}
		//printf("Gauss Fit:   max: %1.2f ± %1.2f | red. chi2 = %1.2f \n",peak[i],peak_err[i],rchi2_g);
		//printf("Landau Fit:   MPV: %1.2f ± %1.2f | red. chi2 = %1.2f \n",lan_mpv_c,lan_mpv_err,lan_rchi2);
    	// printf("Vavilov Fit: MPV: %1.2f         | red. chi2 = %1.2f \n",vav_max,vav_rchi2);

		/***** 
		__  STORE RESULTS IN TXT FILE	_______________________________________
		*****/
		//saves channel number, poisson exp value, poisson exp error
		if (store_channel)
		{
			//saves channel number, mean, mean error
			channel_save_file_mean  << std::fixed << std::setprecision(2) << i << "\t";
			channel_save_file_mean << h_mean << "\t";
			channel_save_file_mean << h_mean_err << "\t";

			//saves channel number, truncated mean, mean error
			channel_save_file_trunc_mean  << std::fixed << std::setprecision(2) << i << "\t";
			channel_save_file_trunc_mean << h_trunc_mean << "\t";
			channel_save_file_trunc_mean << h_trunc_mean_err << "\t";
		}

		if (store_channel && fit_function == 1) //for Poisson
		{
			channel_save_file_poisson  << std::fixed << std::setprecision(2) << i << "\t"; //channel number
			channel_save_file_poisson << expect_p[i] << "\t";
			channel_save_file_poisson << expect_err_p[i] << "\t";
		}

		if (store_channel && fit_function == 5) //for Landau*Gauss
		{
			channel_save_file_langau << std::fixed << std::setprecision(2) << i << "\t"; //channel number
			channel_save_file_langau << SNRPeak << "\t"; //MPV/peak of the fit
			channel_save_file_langau << SNRPeak_err << "\t"; //peak err
		}
	} // ---> end loop over channels + sum

	if (store_channel && fit_function == 1) channel_save_file_poisson << "\n";
	if (store_channel && fit_function == 5) channel_save_file_langau << "\n";
	if (store_channel) channel_save_file_mean << "\n";
	if (store_channel) channel_save_file_trunc_mean << "\n";
	channel_save_file_poisson.close(); //close the file for the channel results
	channel_save_file_langau.close(); //close the file for the channel results
	channel_save_file_mean.close(); //close the file for the channel results
	channel_save_file_trunc_mean.close(); //close the file for the channel results


	// draw the sum
	double h_mean =  h_sum->GetMean();
    double h_max_y = h_sum->GetMaximum();
	TLine * ln_sum; // vertical lines
	TF1 * fit_vec_g_sum, * fit_vec_l_sum, * fit_sum_p; // fit functions
	double amp_max_sum, peak_sum, peak_err_sum, rchi2_g_sum; // gauss fit results
	double expect_p_sum, expect_err_p_sum, rchi2_p_sum, max_p_sum, max_err_p_sum; // Poisson fit results

	C1->cd(n_ch+1);
	h_sum->GetXaxis()->SetTitle("Integrated charge [mV #times ns]");
	h_sum->GetXaxis()->SetTitleOffset(1.3);
	h_sum->GetYaxis()->SetTitle("Entries");
	h_sum->SetLineColorAlpha(kBlack,0.7);
	h_sum->SetFillColorAlpha(kBlack,alpha);
	h_sum->SetMarkerStyle(8);
	h_sum->SetMarkerSize(0.2);
	h_sum->SetMarkerColorAlpha(kBlack,0.6);
	h_sum->Draw("HISTE");


	//Maximum in sum


    // get histogram maximum value & position in range (buggy...)
	amp_max_sum = h_sum->GetXaxis()->GetBinCenter(h_sum->GetMaximumBin());
	float max = h_sum->GetMaximum();
	int lower_bin = h_sum->GetMaximumBin();
	int upper_bin = h_sum->GetMaximumBin();	

	ln_sum = new TLine(amp_max_sum,0.1,amp_max_sum,h_max_y);
	ln_sum->SetLineColor(8);
	ln_sum->SetLineWidth(1);
	ln_sum->SetLineStyle(5);
	ln_sum->Draw("same");

	/***** 
	__ FIT SUM PEAK: GAUSS ______________________________
	*****/
	/*printf("max %1.1f\n",max);
	printf("h_max_y %1.1f\n",h_max_y );

	// fit range ± 0.5x maximum
	while (h_sum->GetBinContent(lower_bin) > max*frac_l) {lower_bin--;}
	while (h_sum->GetBinContent(upper_bin) > max*frac_u) {upper_bin++;}

	float r1,r2;
	r1 = h_sum->GetBinCenter(lower_bin);
	r2 = h_sum->GetBinCenter(upper_bin);
	printf("%1.1f %1.1f\n",r1,r2 );


	fit_vec_g_sum = new TF1("g_fit","gaus",r1,r2);
	fit_vec_g_sum->SetParameter(1,10);
	fit_vec_g_sum->SetLineColor(4);
	fit_vec_g_sum->SetLineStyle(1);
	fit_vec_g_sum->SetNpx(1000); // draw function with high resolution 
	h_sum->Fit("g_fit","RQM");
	h_sum->Draw("sameFUNC");

	//Mean of the Gaussian
	Double_t gauss_mean = fit_vec_g_sum->Mean(r1,r2);

	Double_t par_g[3]; // to store fit results
	fit_vec_g_sum->GetParameters(&par_g[0]);
	// results
	peak_sum = fit_vec_g_sum->GetParameter(1);
	peak_err_sum = fit_vec_g_sum->GetParError(1);

	double sigma_g = fit_vec_g_sum->GetParameter(2);
	double sigma_g_err = fit_vec_g_sum->GetParError(2);
	rchi2_g_sum = fit_vec_g_sum->GetChisquare()/fit_vec_g_sum->GetNDF();


	// draw line for Gaussian Maximum
	TLine * ln_gauss_mpv = new TLine(peak_sum,0,peak_sum,h_max_y);
	ln_gauss_mpv->SetLineColor(4);
	ln_gauss_mpv->SetLineStyle(3);
	ln_gauss_mpv->Draw("same");

	// draw line for Gaussian Mean
	/*TLine * ln_gauss_mean = new TLine(gauss_mean,0,gauss_mean,h_max_y);
	ln_gauss_mean->SetLineColor(2);
	ln_gauss_mean->SetLineStyle(4);
	ln_gauss_mean->Draw("same");*/


	/***** 
	__ FIT SUM: LANDAU GAUSS ______________________________
	*****/
	Double_t SNRPeak_sum, SNRPeak_err_sum, SNRFWHM_sum;
	if (fit_function == 5)
	{
		cout << "\n\n" << endl;
		// set all the parameters for the fit
		Double_t fr[2] {6, 141}; //fit range
		Double_t sv[4], pllo[4], plhi[4];
		Double_t fp[4], fpe[4]; // the final fit parameters and their errors are saved here

		pllo[0]=0.5; pllo[1]=2.0; pllo[2]=0.5; pllo[3]=0.4; // lower parameter limits
		plhi[0]=50.; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=75.0; // upper parameter limits
		//   par[0]=Width (scale) parameter of Landau density
		//   par[1]=Most Probable (MP, location) parameter of Landau density
		//   par[2]=Total area (integral -inf to inf, normalization constant)
		//   par[3]=Width (sigma) of convoluted Gaussian function
		sv[0]=40; sv[1]=20.0; sv[2]=63000.0; sv[3]=25.; // starting values for the fit

		Double_t chisqr_lg;
		Int_t    ndf_lg;
		TF1 *fitsnr = langaufit(h_sum, fr, sv, pllo, plhi, fp, fpe, &chisqr_lg, &ndf_lg);


		// to get the x- and y-values of the maximum
		
		langaupro(fp, SNRPeak_sum, SNRPeak_err_sum, SNRFWHM_sum);
		if (SNRPeak_err_sum < 0) SNRPeak_err_sum *= -1;

		printf("Fitting of the sum histogram done\nPlotting results...\n");

		cout << "Langaus Peak of the sum = " << SNRPeak_sum << "#pm" << SNRPeak_err_sum << endl;

		fitsnr->Draw("lsame");

		//draw line for maximum
		TLine * langaus_max = new TLine(SNRPeak_sum,0,SNRPeak_sum,h_max_y);
		langaus_max->SetLineColor(2);
		langaus_max->SetLineStyle(4);
		langaus_max->Draw("same");
	}


	/***** 
	__ FIT SUM: GENERALIZED POISSON ______________________________
	*****/
	if (fit_function == 2)
	{
		double poisson_range_lo_sum = 9, poisson_range_hi_sum = 130;
		
		fit_sum_p = new TF1("p_fit", generalized_poisson, poisson_range_lo_sum, poisson_range_hi_sum, 3); // 0 just normalization, 1 = theta, 2 = lambda
		fit_sum_p->SetParameters(1000, h_mean, 1);
		fit_sum_p->SetParNames("Norm.", "Theta", "Lambda");
		fit_sum_p->SetParLimits(2, 0., 1.);
		fit_sum_p->SetLineColor(6);
		fit_sum_p->SetLineStyle(1);
		fit_sum_p->SetNpx(1000); // draw function with high resolution 
		h_sum->Fit("p_fit","RM+");
		h_sum->Draw("same");

		Double_t par_p[3]; // to store fit results
		Double_t par_err_p[3];
		fit_sum_p->GetParameters(&par_p[0]);
		for (int j = 0; j<3; j++) par_err_p[j] = fit_sum_p->GetParError(j);


		/*for (int j=0; j<3; j++) {
			cout << "par_p[" << j << "] = " << par_p[j] << endl; // to print fit results
		cout << "par_error_p[" << j << "] = " << par_error_p[j] << endl; // to print fit errors
		}*/

		// expectation value results
		expect_p_sum = par_p[1] / (1 - par_p[2]); //formula for the expectation value of the generalized Poisson
		expect_err_p_sum = sqrt(1./(pow((1 - par_p[2]), 2)) * pow(par_err_p[1], 2) + par_p[1] / pow((1 - par_p[2]), 4) * pow(par_err_p[2], 2)); //formula from Gaussian error propagation

		cout << "expect_p_sum = " << expect_p_sum << endl;
		cout << "expect_err_p_sum = " << expect_err_p_sum << endl;

		// draw line for Poisson Expectation value
		TLine * poisson_exp = new TLine(expect_p_sum, 0, expect_p_sum, 1100);
		poisson_exp->SetLineColor(6);
		poisson_exp->SetLineStyle(4);
		poisson_exp->Draw("same");

		rchi2_p_sum = fit_sum_p->GetChisquare()/fit_sum_p->GetNDF();

		// maximum
		max_p_sum = fit_sum_p->GetMaximumX(poisson_range_lo_sum + 1, poisson_range_hi_sum - 1);
		max_err_p_sum = (poisson_range_hi_sum - poisson_range_lo_sum) / 1000; // error of the max is one bin of the function

		// draw line for Poisson maximum value
		/*TLine * poisson_max = new TLine(max_p_sum, 0, max_p_sum, h_max_y);
		poisson_max->SetLineColor(6);
		poisson_max->SetLineStyle(5);
		poisson_max->Draw("same");


		



		cout << "Das Maximum des Poissonfits ist " << max_p_sum << "#pm"  << max_err_p_sum << endl;*/
	}

	/***** 
	__ MEDIAN OF THE SUM ______________________________
	*****/
	Double_t median_x_sum, q;
	q = 0.5;
	h_sum->GetQuantiles(1,&median_x_sum, &q);
	//cout << "Der Median der Summe wurde zu " << median_x_sum << " berechnet." << endl;
	TLine * ln_median = new TLine(median_x_sum, 0, median_x_sum ,h_max_y);
	ln_median->SetLineColor(9);
	ln_median->SetLineStyle(9);
	ln_median->Draw("same");


	/***** 
	__  LEGEND	_______________________________________
	*****/
	TLegend *h_lan_leg = new TLegend(0.33,0.50,0.9,0.90);
	h_lan_leg->SetTextSize(0.035);
	h_lan_leg->AddEntry(h_sum,Form("#bf{data}"),"elpf");
	h_lan_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_sum->GetEntries()),"");
	h_lan_leg->AddEntry(ln_sum,Form("Dist. max.: %1.1f",amp_max_sum),"l");
	h_lan_leg->AddEntry((TObject*)0,Form("Dist. mean = %1.2f #pm %1.2f",h_sum->GetMean(),h_sum->GetStdDev()/sqrt(h_sum->GetEntries())),"");
	h_lan_leg->AddEntry(ln_median,Form("Dist. median: %1.2f",median_x_sum),"l");
	if (fit_function == 5) //for Landau*Gauss
	{
		h_lan_leg->AddEntry(fitsnr,Form("Landau*Gauss Fit"),"l");
		h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.2f", chisqr_lg/ndf_lg),"");
		h_lan_leg->AddEntry(langaus_max,Form("Langau fit Maximum = %1.2f #pm %1.2f", SNRPeak_sum, SNRPeak_err_sum), "l");
	}
	//h_lan_leg->AddEntry(fit_vec_g_sum,Form("Gaussian fit"),"l");
	//h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.1f",rchi2_g_sum),"");
	//h_lan_leg->AddEntry(ln_gauss_mpv,Form("Gaussian mean = %1.2f #pm %1.2f ",peak_sum,peak_err_sum),"l");
	//h_lan_leg->AddEntry((TObject*)0,Form("#sigma = %1.2f #pm %1.2f",sigma_g,sigma_g_err),"");
	if (fit_function == 2) //for Generalized Poisson
	{
		h_lan_leg->AddEntry(fit_sum_p,Form("Generalized Poisson fit"),"l");
		h_lan_leg->AddEntry((TObject*)0,Form("Generalized Poisson #chi^{2}/ndf = %1.1f",rchi2_p_sum),"");
		h_lan_leg->AddEntry(poisson_exp,Form("Gen. Poisson exp. value = %1.2f #pm %1.2f", expect_p_sum, expect_err_p_sum),"l");
	}
	//h_lan_leg->AddEntry(poisson_max,Form("Poisson max = %1.2f #pm %1.2f", max_p_sum, max_err_p_sum),"l");
	//h_lan_leg->AddEntry(ln_gauss_mean,Form("Gaussian mean = %1.2f ",gauss_mean),"l");
	
	if(legend) h_lan_leg->Draw();


	/***** 
	__ EXPORT ______________________________
	*****/

//	printf("%s %d %d %s \n",particle.c_str(),run_nr,pos,WC_v.c_str() );

	if (store_result)
	{
		/*FILE *mpv_list, *mpv_err_list;
		mpv_list = fopen(Form("%s%s",out_dir.c_str(),out_values_list_name.c_str()),"a");
		mpv_err_list =fopen(Form("%s%s",out_dir.c_str(),out_values_err_list_name.c_str()),"a");

		//   // print date to file
		// time_t now;
		// time(&now); 
		// fprintf(mpv_list,"\nlight yield analysis - %s\n",ctime(&now));
		// fprintf(mpv_err_list,"\nlight yield analysis - %s\n",ctime(&now));

		for (int i = 0; i < n_ch; ++i){fprintf(mpv_list, "%f ",peak[i] );}
		for (int i = 0; i < n_ch; ++i){fprintf(mpv_err_list, "%f ",peak_err[i] );}

		//fprintf(mpv_list, "%s %d %d %s\n",particle.c_str(),run_nr,pos,WC_v.c_str());
		//fprintf(mpv_err_list, "%s %d %d %s\n",particle.c_str(),run_nr,pos,WC_v.c_str());

*/
		gErrorIgnoreLevel = kError; // suppress root terminal output 
		C1->SaveAs(Form("%s%s",out_dir.c_str(),out_pdf_name.c_str()));
		gErrorIgnoreLevel = kUnset;

	}


	// saving to file (in that order) of the sum histogram: run number, distribution max, distribution mean, distribution median, Gaussian max, Gaussian max error,
	// Poisson expectation value, Poisson expectation value error, Poisson max value, Poisson max error
	ofstream save_file ("sum_values.txt", ios::app);
	if (save_file.is_open())
	{
		//run number, i.e. first three letters of run_name
		save_file << run_name[0];
		save_file << run_name[1];
		save_file << run_name[2] << "\t";

		//distribution max
		save_file << std::fixed << std::setprecision(2) << amp_max_sum << "\t";

		//Distribution mean
		save_file << h_sum->GetMean() << "\t";

		//Distribution median
		save_file << median_x_sum << "\t";

		//Gaussian max
		/*save_file << peak_sum << "\t";

		//Gaussian max error 
		save_file << peak_err_sum << "\t";

		//Poisson expectation value
		save_file << expect_p_sum << "\t";

		//Gaussian expectation value error 
		save_file << expect_err_p_sum << "\t";

		//Poisson max value
		save_file << max_p_sum << "\t";

		//Poisson max value error 
		save_file << max_err_p_sum << "\n";*/

		if (fit_function == 5)
		{
			//Langau MPV
			save_file << SNRPeak_sum << "\t";

			//Langau MPV error
			save_file << SNRPeak_err_sum << "\n";
		}

		save_file.close(); 
	}
	else cout << "Unable to open file" << endl;

	if (is_interactive){ROOTapp->Run();}
	

  	return 0;
}