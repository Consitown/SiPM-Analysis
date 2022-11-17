//root
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>		  // root verbosity level
#include <TApplication.h> // open root gui

//C, C++
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

	/***** 
	INITIALIZE 
	*****/

	// ROOT GUI
	bool is_interactive = 0;
	// amp (0), charge (1)
	bool is_charge = 1;
	// store results
	bool store_result = 1;

	TApplication *ROOTapp;
	if (is_interactive)
	{
		ROOTapp = new TApplication("ROOT interactive", &argc, argv);
	}

	// style settings
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleW(1.0);
	// gStyle->SetTitleH(1.0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	// gStyle->SetImageScaling(30.);
	gStyle->SetLineScalePS(1);

	// parse arguments
	string run_name = (string)argv[1];
	int run_nr = atoi(argv[2]);
	double x_lim = atof(argv[3]);
	string fileLocation = string(argv[4]);

	string in_dir = fileLocation;				 // path to current root file
	string out_dir = "./BL_histograms/";		 // to export plots and values
	string out_list = "BL_values_run_" + run_nr; // filename values list

	printf("\nANALYZE CHARGE: run%d\n", run_nr);

	// open root file
	TFile *file = new TFile(in_dir.c_str());
	TTree *tree;
	file->GetObject("T", tree);

	int Xmin = 0;
	int Xmax = 10;
	int binX = 320;
	float bin_w = (float)Xmax / binX;

	float alpha = 0.1;
	int lineW = 1;

	/***** 
	AVERAGE BASELINE 
	*****/

	TCanvas *C_vec[16];
	TH1F *h_vec[16];
	TLine *ln_vec[8];
	TF1 *fit_BL_vec[8];
	float chi2_min[8], mean_BL[8], mean_BL_err[8], sigma_BL[8], sigma_BL_err[8];

	TCanvas *C1 = new TCanvas("BL_Chi2_Dist", "", 1000, 1000);
	TCanvas *C2 = new TCanvas("BL_Value_Dist", "", 1000, 1000);

	C1->Divide(3, 3);
	C2->Divide(3, 3);

	for (int i = 0; i < 8; ++i)
	{

		TString c_name1, c_name2, h_name1, h_name2, draw_cmnd1, draw_cmnd2, cut_cmnd2;
		TString h_title1, h_title2;
		TString h_xtile1, h_xtile2, h_ytitle1, h_ytitle2;
		TString g_fit;
		TString xTitle2;

		int channel = i;

		/***** 
		__ SHOW BL_CHI2 ______________________________
		*****/

		C1->cd(i + 1);
		gPad->SetLogy();
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetRightMargin(0.00);
		gPad->SetLeftMargin(.12);
		gPad->SetBottomMargin(.12);

		c_name1.Form("c%d", i);
		h_name1.Form("h%d", i);
		h_title1.Form("ch%d", i + 1);
		draw_cmnd1.Form("BL_Chi2_used[%d]>>h%d", i, i);

		h_vec[i] = new TH1F(h_name1, h_title1, binX, Xmin, Xmax);
		h_vec[i]->SetLineColor(kBlack);
		h_vec[i]->SetLineWidth(lineW);
		h_vec[i]->SetFillColorAlpha(kBlack, alpha);
		h_vec[i]->SetMarkerStyle(8);
		h_vec[i]->SetMarkerSize(0.1);
		h_vec[i]->SetMarkerColorAlpha(kBlack, 0.6);
		h_vec[i]->SetTitle(Form("baseline fit, red. #chi^{2} , ch%d", channel));
		h_vec[i]->GetXaxis()->SetTitle("#chi_{red}^{2}");
		h_vec[i]->GetYaxis()->SetTitle("counts");
		h_vec[i]->GetYaxis()->SetTitleSize(0.05);
		h_vec[i]->GetXaxis()->SetTitleSize(0.05);

		tree->Draw(draw_cmnd1, "", "HISTE");

		/***** 
		__ MINIMUM BL_CHI2 in range_________________
		*****/

		h_vec[i]->GetXaxis()->SetRange(1. / bin_w, x_lim / bin_w);
		chi2_min[i] = h_vec[i]->GetXaxis()->GetBinCenter(h_vec[i]->GetMinimumBin());
		// printf("chi2_min_1: %f\n",chi2_min[i]);
		h_vec[i]->GetXaxis()->UnZoom();
		// h_vec[i]->GetXaxis()->SetRange(0,5/bin_w);

		ln_vec[i] = new TLine(chi2_min[i], 10, chi2_min[i], 2000);
		ln_vec[i]->SetLineColor(2);
		ln_vec[i]->SetLineWidth(lineW);
		ln_vec[i]->Draw("same");

		TLegend *h_chi2_leg = new TLegend(0.53, 0.70, 1.0, 0.9);
		h_chi2_leg->SetTextSize(0.04);
		h_chi2_leg->AddEntry(h_vec[i], Form("#bf{#chi_{red}^{2}, baseline fit}"), "elpf");
		//h_chi2_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i]->GetEntries()),"");
		h_chi2_leg->AddEntry(ln_vec[i], Form("#chi_{cut}^{2}= %1.2f", chi2_min[i]), "l");
		h_chi2_leg->Draw("");

		/***** 
		__ SHOW BL ___________________________________
		*****/

		c_name2.Form("c%d", i + 8);
		h_name2.Form("h%d", i + 8);
		h_title2.Form("ch%d", i + 1);
		draw_cmnd2.Form("BL_used[%d]>>h%d", i, i + 8);
		cut_cmnd2.Form("BL_Chi2_used[%d]<%f", i, chi2_min[i]);
		xTitle2.Form("baseline (cut: #chi_{red}^{2}<%1.2f) [mV]", chi2_min[i]);
		g_fit.Form("fit ch%d", i);

		C2->cd(i + 1);
		gPad->SetRightMargin(0.00);
		gPad->SetLeftMargin(.12);
		// gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
		gPad->SetBottomMargin(.12);

		h_vec[i + 8] = new TH1F(h_name2, h_title2, binX, -5, 8);
		h_vec[i + 8]->SetLineColor(kBlue);
		h_vec[i + 8]->SetLineWidth(lineW);
		h_vec[i + 8]->SetFillColorAlpha(kBlue, alpha);
		h_vec[i + 8]->SetLineColor(kBlack);
		h_vec[i + 8]->SetLineWidth(lineW);
		h_vec[i + 8]->SetFillColorAlpha(kBlack, alpha);
		h_vec[i + 8]->SetMarkerStyle(8);
		h_vec[i + 8]->SetMarkerSize(0.5);
		h_vec[i + 8]->SetMarkerColorAlpha(kBlack, 0.6);
		h_vec[i + 8]->SetTitle(Form("baseline fit, ch%d", channel));

		h_vec[i + 8]->GetYaxis()->SetTitle("counts");
		h_vec[i + 8]->GetXaxis()->SetTitle(xTitle2);
		h_vec[i+8]->GetYaxis()->SetTitleSize(0.05);
		h_vec[i+8]->GetXaxis()->SetTitleSize(0.05);
		tree->Draw(draw_cmnd2, cut_cmnd2, "goff");

		/***** 
		__ FIT AVERAGE BL ___________________________________
		*****/

		// gauss fit around histogram mean
		float h_mean = h_vec[i + 8]->GetMean();
		float fit_range = 1.0;

		fit_BL_vec[i] = new TF1(g_fit, "gaus", h_mean - fit_range, h_mean + fit_range);
		h_vec[i + 8]->Fit(g_fit, "RQ");

		mean_BL[i] = fit_BL_vec[i]->GetParameter(1);
		mean_BL_err[i] = fit_BL_vec[i]->GetParError(1);
		sigma_BL[i] = fit_BL_vec[i]->GetParameter(2);
		sigma_BL_err[i] = fit_BL_vec[i]->GetParError(2);


		//Alternative
		/*mean_BL[i] = h_vec[i + 8]->GetMean();
		mean_BL_err[i] = h_vec[i + 8]->GetMeanError();
		sigma_BL[i] = 0;*/
		sigma_BL_err[i] = 0;

    	
		h_vec[i + 8]->Draw();


		TLine *ln4 = new TLine(mean_BL[i] , h_vec[i + 8]->GetMinimum(), mean_BL[i] ,h_vec[i + 8]->GetMaximum());
		ln4->SetLineColor(2);
        ln4->SetLineWidth(3);
		//ln4->Draw();

		TLegend *h_BL_leg = new TLegend(0.55, 0.65, 1.0, 0.9);
		h_BL_leg->SetTextSize(0.035);
		h_BL_leg->AddEntry(h_vec[i + 8], Form("#bf{B_{dyn}(#chi_{red}^{2} #leq #chi_{cut}^{2}) }"), "elpf");
		//h_BL_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i+8]->GetEntries()),"");
		h_BL_leg->AddEntry(fit_BL_vec[i], Form("Gaussian fit:"), "l");
		h_BL_leg->AddEntry((TObject *)0, Form("B_{const}= %1.3f #pm %1.3f", mean_BL[i], mean_BL_err[i]), "");
	//	h_BL_leg->AddEntry((TObject *)0, Form("#sigma_{BL} = %1.3f #pm %1.3f", sigma_BL[i], sigma_BL_err[i]), "");
		h_BL_leg->Draw("");

	} // end loop over channels

	/***** 
	__ PRINT RESULTS ___________________________________
	*****/

	// for (int i = 0; i < 8; ++i)
	// {
	// 	printf("Minimum Chi2 in range(0.7,1.5) ch%d = %f\n",i,chi2_min[i] );
	// }
	// for (int i = 0; i < 8; ++i)
	// {
	// 	printf("Mean Baseline ch%d = %f\n",i+1,mean_BL[i] );
	// }

	printf("Mean Baseline ch1-8 = %f,%f,%f,%f,%f,%f,%f,%f\n", mean_BL[0], mean_BL[1], mean_BL[2], mean_BL[3], mean_BL[4], mean_BL[5], mean_BL[6], mean_BL[7]);

	FILE *pFile = fopen("Baselines.txt", "a");
	fprintf(pFile, "BL_%s={%f,%f,%f,%f,%f,%f,%f,%f}\n", string(run_name).c_str(), mean_BL[0], mean_BL[1], mean_BL[2], mean_BL[3], mean_BL[4], mean_BL[5], mean_BL[6], mean_BL[7]);
	fclose(pFile);

	string pdf_filename1, pdf_filename2;
	pdf_filename1 = out_dir + run_name + Form("_fit.pdf");
	pdf_filename2 = out_dir + run_name + Form("_values.pdf");

	gErrorIgnoreLevel = kError; // suppress root terminal output
	C1->Print(pdf_filename1.c_str());
	C2->Print(pdf_filename2.c_str());
	gErrorIgnoreLevel = kUnset;
	delete C1;
	delete C2;

	if (is_interactive)
	{
		ROOTapp->Run();
	}
	return 0;
}