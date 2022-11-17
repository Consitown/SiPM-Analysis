//root
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
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

//C, C++
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

//Wie viele Peaks man anschauen will.
const int n_peaks = 6;

//Values-> Blaue Punkte schwankungen sind indikator für Fehler des Roten Wertes (Abschätzung, hier Standardabweichung), mehr Peaks-> höhrere Werte-> systematisches PRoblem auf Elektronik zurückführbar.

//Continousfit-> was wir im Endeffekt nehmen, wir nehmen fitf-> trotz 
//Alt Fit-> addition der 6 Gauß Peaks-> abstände in Valus.pdf Blau->Steigen-> Schwierig, man könnte auch Mittelwert bilden, aber alles korreliert-> aufwendig
//Rot-> FitF wert

// define a fit function, with single parameter corresponding to ptp distance
//Eine Art Vergleichsfunktion, noch mit Gaußfits.
//Nimmt Abstand der 0-1 Peaks und nimmt an, dass der überall gleich ist.-> Darauf stützen wir unsere Kalibration (Gain=Const, sonst für jede Intensität neue Kalib).
//NEW FIT, SCHLECHTER
//PArameter in PDF-> p hier 
Double_t fitf(Double_t *x,Double_t *p)
{
	
	Double_t pos[n_peaks];
	Double_t gaus_single[n_peaks];

	Double_t dP = p[1]; // distance of second minus first peak

	for (int i = 0; i < n_peaks; ++i)
	{
		pos[i]=p[0]+dP*i; 
		gaus_single[i] = p[2+2*i]*TMath::Exp(-(pos[i]-x[0])*(pos[i]-x[0])/2/p[3+2*i]/p[3+2*i]);
	}
	double gaus_comb = 0;
	for (int i = 0; i < n_peaks; ++i)
	  {
	  	gaus_comb = gaus_comb + gaus_single[i];
	  }  
	 return gaus_comb;
}

// fit function, sum of n_peaks Gaussians 
//ALT FIT
//Addition von Gaußpeaks, jeder Parameter muss vorher Range haben, sonst wird das nix. -> Zerschneide Histogramm, fitte Peaks einzeln, speichere Daten und nehme das als Range
Double_t alt_f(Double_t *x,Double_t *p)
{
	Double_t gaus_single[n_peaks];
	for (int i = 0; i < n_peaks; ++i)
	{
		gaus_single[i] = p[0+3*i]*TMath::Exp(-TMath::Power( p[1+3*i]-x[0],2) / (2*TMath::Power(p[2+3*i],2)) );
	}

	double gaus_comb = 0;
	for (int i = 0; i < n_peaks; ++i)
	  {
	  	gaus_comb = gaus_comb + gaus_single[i];
	  }  
	return gaus_comb;
}

/********************
__ FIT ROUTINE______
********************/

int main(int argc, char *argv[]){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetStatW(0.2); gStyle->SetStatH(0.1); // stats box size
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);


	/***** 
	__ INITIALIZE ___________________________
	*****/

	string path = "./"+(string)argv[1]+".root";
	printf("PATH NAME ________________ %d", path);
	int channel = atoi(argv[2]);
	TFile* file = new TFile(path.c_str());
	TTree* tree;
	file->GetObject("T",tree);
	TCanvas* c = new TCanvas("c","",800,600);
 
	int nBins = 135;
	int xmin = 0;
	int xmax = 45;

	TString cut(Form("BL_Chi2_lower[%d]<1.4",channel));
	TH1F* h = new TH1F("h","PE spectrum - amplitude window 20 ns",nBins,xmin,xmax);
	// amp_inRange: max. amplitude in 20 ns window from 145 ns tp 165 ns
	tree->Draw(Form("amp_inRange[%d]>>h",channel));

	TCanvas *C2 = new TCanvas("C2","alt fit",800,600);
	TH1F* h2 = new TH1F("h2",";amp_inRange [mV];",nBins,xmin,xmax);
	h2 = (TH1F*)h->Clone();
	// tree->Draw(Form("amp_inRange[%d]>>h2",channel),"","goff");

	/***** 
	__ FIT RANGE ___________________________
	*****/
	
	// read off spectrum, ranges for individual peaks
	

	//Range für die Peaks, ablesen an Histogramm.
	std::vector<int> ranges = {5,10,15,21,27,33,39}; // run 65
	// std::vector<int> ranges = {3,7,12,16,21,26,31}; // run 77,76

	// upper and lower range 
	//Für ContinousFit, für 
	float l_range = 1.5;
  	float u_range = 1.5;

	/***** 
	__ PRE FIT ___________________________
	*****/
	
	// get starting value for gauss mean parameter
	float pos_peak[n_peaks];
	for (int i = 0; i < n_peaks; ++i)
	{
		h->GetXaxis()->SetRange( h->GetXaxis()->FindBin(ranges[i]) ,h->GetXaxis()->FindBin(ranges[i+1]));
		pos_peak[i] = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	}	
	h->GetXaxis()->SetRange(0,nBins);

	TF1 * peak_single[n_peaks];
	Double_t par_single[3*n_peaks];
	for (int i = 0; i < n_peaks; ++i)
	{	
		float range = 1.5; // symmetric fit range for individual peaks
		peak_single[i] = new TF1("peak","gaus",pos_peak[i]-range,pos_peak[i]+range);
		h->Fit("peak","RNQ+");
		peak_single[i]->GetParameters(&par_single[0+3*i]);
	}

	/***** 
	__ CONTINUOUS FIT ___________________________
	*****/

	TF1 *comb = new TF1("comb",fitf, pos_peak[0]-l_range, pos_peak[n_peaks-1]+u_range, 2+2*n_peaks);
	comb->SetParameter(0,par_single[1]); // peak_1, param_1
	comb->SetParameter(1,par_single[4]-par_single[1]); // distance: second-first peak

	for (int i = 0; i < n_peaks; ++i)
	{	
		comb->SetParameter(2+2*i,par_single[0+3*i]); // peak_i, param_0
		comb->SetParameter(3+2*i,par_single[2+3*i]); // peak_i, param_2
	}
	
	comb->SetLineColor(3);
	h->Fit("comb","RQM");
	h->Draw();

	/***** 
	__ STORE RESULTS ___________________________
	*****/

	// save parameters form contiuous fit
	Double_t par_comb[2+2*n_peaks];
	Double_t par_comb_err[2+2*n_peaks];
	comb->GetParameters(&par_comb[0]);
	for (int i = 0; i < 2+2*n_peaks; ++i)
	{
		par_comb_err[i] = comb->GetParError(i);
	}

	Double_t calib_factor = comb->GetParameter(1);
	Double_t calib_factor_err = comb->GetParError(1);

	double comb_chi2_ndof = comb->GetChisquare()/comb->GetNDF();

	printf("ch%d : 1pe dist: %1.2f ± %1.2f | red. chi2 = %1.2f\n",channel,calib_factor,calib_factor_err,comb_chi2_ndof );
	
	// print indiviual peak fits
	for (int i = 0; i < n_peaks; ++i)
	{
		peak_single[i]->Draw("same");
	}

	// draw individual peaks from contiuous fit
	for (int i = 0; i < n_peaks; ++i)
	{
		TF1* s_peak = new TF1("s_peak","gaus",par_comb[0]+par_comb[1]*i-5,par_comb[0]+par_comb[1]*i+5);
		s_peak->SetParameter(0,par_comb[2+2*i]);
		s_peak->SetParameter(1,par_comb[0]+par_comb[1]*i);
		s_peak->SetParameter(2,par_comb[3+2*i]);

		s_peak->SetLineColor(4);
		s_peak->Draw("same");
	}

	comb->Draw("same");	

	/***** 
	ALT FIT
	__ ALTERNATIVE FIT ___________________________
	*****/

	TF1 *alt = new TF1("alt",alt_f, pos_peak[0]-l_range, pos_peak[n_peaks-1]+u_range, 3*n_peaks);
	alt->SetParameters(&par_single[0]);
	// for (int i = 0; i < 3*n_peaks+3; ++i)
	// {
	// 	alt->SetParameter(i,par_single[i]);
	// }
	h2->Fit("alt","RQM");

	alt->SetLineColor(3);
	
	h2->Draw();

	Double_t par_alt[3+n_peaks];
	Double_t par_alt_err[3+n_peaks];
	// alt->GetParameters(&par_alt[0]);
	for (int i = 0; i < 3*n_peaks; ++i)
	{
		par_alt[i] = alt->GetParameter(i);
		par_alt_err[i] = alt->GetParError(i);
	}

	double pos_alt[n_peaks], u_pos_alt[n_peaks];
	for (int i = 0; i < n_peaks; ++i)
	{
		pos_alt[i] = alt->GetParameter(1+3*i);
		u_pos_alt[i] = alt->GetParError(1+3*i);
	}
	
	// printf("ch%d: ",channel );
	// for (int i = 0; i < n_peaks; ++i)
	// {
	// 	printf("%1.2f±%1.2f ",pos_alt[i],u_pos_alt[i] );
	// }
	// printf("\n");

	double diff[n_peaks-1], u_diff[n_peaks-1], wght[n_peaks-1];
	double w_mean, u_w_mean, sum1, sum2;
	for (int i = 0; i < n_peaks-1; ++i)
	{
		diff[i] = pos_alt[i+1]-pos_alt[i];
		u_diff[i] = sqrt( u_pos_alt[i+1]*u_pos_alt[i+1] + u_pos_alt[i]*u_pos_alt[i] );
		wght[i] = 1./( u_diff[i] * u_diff[i] );
		sum1 += diff[i]*wght[i];
		sum2 += wght[i];
	}

	w_mean = sum1/sum2;
	u_w_mean = sqrt (1./sum2);

	/***** 
	__ COMPARE ALT <-> COMB ___________________________
	*****/

	double x_values[n_peaks-1];
  	double x_values_err[n_peaks-1];
  	for (int i = 0; i < n_peaks-1; ++i)
	{
	x_values[i] = i+1;
	x_values_err[i] = 0;
	}


	TCanvas *C1 = new TCanvas("C1","Compare Fit Results",700,500);
	// C1->SetGrid();
	gPad->SetGridx(); gPad->SetGridy();
	TGraphErrors * gr_alt = new TGraphErrors(n_peaks-1,x_values,diff,x_values_err,u_diff);
	gr_alt->SetName("gr_alt");
	gr_alt->SetTitle("PE Spectrum: Peak-To-Peak Distance");
   	gr_alt->SetMarkerColor(4);
   	gr_alt->SetMarkerStyle(21);
   	gr_alt->GetXaxis()->SetTitle("peak number");
   	gr_alt->GetYaxis()->SetTitle("peak-to-peak distance in mV");
   	gr_alt->Draw("AP");


   	double x_pint= 2.5;
   	TGraphErrors * gr_comb = new TGraphErrors(1,&x_pint,&calib_factor,&x_values_err[1],&calib_factor_err);
   	gr_comb->SetName("gr_comb");
   	gr_comb->SetMarkerColor(2);
   	gr_comb->SetMarkerStyle(11);
   	gr_comb->Draw("P");
	
   	TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
	// leg->SetHeader("The Legend Title");
	// leg->AddEntry(h1,"Histogram filled with random numbers","f");
	// leg->AddEntry("gr_alt","Function abs(#frac{sin(x)}{x})","l");
	leg->AddEntry("gr_comb","#splitline{average distance}{single parameter gauss fit}","lep");
	leg->AddEntry("gr_alt","#splitline{individual peak distances}{reg. cont. gauss fit}","lep");
	leg->Draw();



	// ____ Print Info _____

	// bool print_more = false; 
	// if (print_more)
	// {
	// 	cout << Form("ch%d:Peak Maxima (max. Bin): %f %f %f %f %f",channel, pos_peak1, pos_peak2, pos_peak3, pos_peak4, pos_peak5) << endl;
	// 	cout << Form("ch%d:Peak Maxima (Gauss Mean): %f %f %f %f %f",channel, pos_peak1, pos_peak2, pos_peak3, pos_peak4,pos_peak5) << endl;
	// 	cout << Form("ch%d:Peak Maxima Uncert. (Gauss Mean): %f %f %f %f %f",channel, pos_peak1_err, pos_peak2_err, pos_peak3_err, pos_peak4_err, pos_peak5_err) << endl;
	// }
	
	/***** 
	__ Export Result ___________________________
	*****/

	string target_path = "./calib_histograms/amplitude/";
	string pdf_filename1 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_new_fit.pdf";
	string pdf_filename2 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_alt_fit.pdf";
	string pdf_filename3 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_values.pdf";
	string list_filename = target_path+"calib_factor.txt";

	FILE * factor_list;
	factor_list = fopen(list_filename.c_str(),"a");
	fprintf(factor_list, "ch%d : 1pe dist: %f ± %f | red. chi2 = %f\n",channel,calib_factor,calib_factor_err,comb_chi2_ndof );
	fclose(factor_list);

	gErrorIgnoreLevel = kError; // suppress root terminal output 
  	c->Print(pdf_filename1.c_str());
  	C1->Print(pdf_filename3.c_str());
  	C2->Print(pdf_filename2.c_str());

  	return 0;
}