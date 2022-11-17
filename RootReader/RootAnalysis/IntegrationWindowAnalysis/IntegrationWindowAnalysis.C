//root
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TRandom.h>
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
#include <experimental/filesystem>
#include <TMultiGraph.h>
#include <TF1Convolution.h>
#include <TAttMarker.h>

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
#include <sys/types.h>
#include <sys/stat.h>
#include <TGraph.h>
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TApplication.h> // open root gui

/**
 *  Tutorial
 *  1.constant Baseline -> read in binary data with dynamic BL -> everything else is arbitrary
 * 	2.create a constant Baseline for this run -> Baseline Script-> Move Baselines.txt to RunHelper/src
 * 	3.now, every run is done with the constant BL setting on (RunHelper)
 * 	4.create correction factor -> If working with integrals
 * 		4.1 ->use any generated root files (only using the sumHistograms inside the file) and run this script,
 * 			  it will create multiple files: IntegrationWindows.txt (contains the windows) and IW.pdf -> shows if the windows are calculated correctly
 * 		4.2 ->move the IntegrationWindows.txt to RunHelper/src, the new runs will now use the correct integration windows bases on the sumHistograms
 * 		4.3 ->run the RunHelper one more time, this time to calculate the Percentage Distribution, used to figure out the Correction factor uncertainty
 * 		4.4 ->Inside of your root file there is now a branch called "Percentages" holding a distribution of SignalPeak/SignalAll
 * 		4.5 ->Run this script, it will use this branch to calculate the uncertainty of the Percentagevalue and the correctionvalue + uncertainty
 * 	5.move CorrectionValues.txt to RunHelper/src, it will automatically use the correctionValues and multiply them with the integral
 * 
 * 
 * 
 * 
 * REMINDER:
 * There is a calibration Runnumber hardcoded below. You should make sure, that this calibration run is always present in the txt files.
 * 
 * 
 * */

namespace fs = std::experimental::filesystem;
using namespace std;
struct stat info;
float SP = 0.3125; // ns per bin
float IntegralHist(TH1F *hWave, float t1, float t2, float BL)
{
	float BW = hWave->GetXaxis()->GetBinWidth(1);
	int bin1 = hWave->FindBin(t1);
	int bin2 = hWave->FindBin(t2);
	float c1 = hWave->GetBinContent(bin1);
	float c2 = hWave->GetBinContent(bin2);
	return hWave->Integral(bin1, bin2, "width") - BL * (t2 - t1) - c1 * (t1 - hWave->GetXaxis()->GetBinLowEdge(bin1)) - c2 * (hWave->GetXaxis()->GetBinUpEdge(bin2) - t2);
}

float *BL_fit(TH1F *hWave, float *BL_chi2, float t1, float t2)
{

	// TF1 *f_const = new TF1("f_const","[0]",t1,t2);
	// hWave->Fit("f_const","RN");
	TF1 *f_const = new TF1("f_const", "pol0", t1, t2);
	hWave->Fit("f_const", "RNWQ");

	BL_chi2[0] = f_const->GetParameter(0);
	BL_chi2[1] = f_const->GetParError(0);
	BL_chi2[2] = f_const->GetChisquare() / f_const->GetNDF();
	BL_chi2[3] = f_const->GetProb();

	return BL_chi2;
}

float t_max_inRange(TH1F *hWave, float t1, float t2)
{

	hWave->GetXaxis()->SetRange(t1 / SP, t2 / SP);
	float t_max = hWave->GetXaxis()->GetBinCenter(hWave->GetMaximumBin());
	hWave->GetXaxis()->SetRange(1, 1024);
	return t_max;
}
float t_min_inRange(TH1F *hWave, float t1, float t2)
{

	hWave->GetXaxis()->SetRange(t1 / SP, t2 / SP);
	float t_max = hWave->GetXaxis()->GetBinCenter(hWave->GetMinimumBin());
	hWave->GetXaxis()->SetRange(1, 1024);
	return t_max;
}
string extractValues(string str)
{
	unsigned first = str.find("{");
	unsigned last = str.find("}");
	string strNew = str.substr(first + 1, last - first - 1);
	return strNew;
}
float stringToFloat(string text)
{
	try
	{
		if (text.length() == 0)
			return 0.0;
		return stof(text);
	}
	catch (const std::exception &e)
	{
		return 0.0;
	}
}
pair<vector<float>, vector<float>> readPairs(string line, double initValue)
{
	vector<float> peakSignals;
	vector<float> allSignals;

	peakSignals.clear();
	allSignals.clear();

	string s = extractValues(line);
	string columnDelimiter = ",";
	string valueDelimiter = "/";

	size_t pos = 0;
	std::string token;
	string calibValue;
	while ((pos = s.find(columnDelimiter)) != std::string::npos)
	{
		token = s.substr(0, pos); //LIKE: 33/178

		int posOfDelimiter = token.find("/");

		string peak = token.substr(0, posOfDelimiter);
		string all = token.substr(posOfDelimiter + 1, token.length() - 1);

		peakSignals.push_back(stringToFloat(peak));
		allSignals.push_back(stringToFloat(all));

		s.erase(0, pos + columnDelimiter.length());
	}
	int posOfDelimiter = s.find("/");
	string peak = s.substr(0, posOfDelimiter);
	string all = s.substr(posOfDelimiter + 1, s.length() - 1);

	peakSignals.push_back(stringToFloat(peak));
	allSignals.push_back(stringToFloat(all));

	for (std::size_t i = peakSignals.size(); i < 32; i++)
	{
		peakSignals.push_back(initValue);
		allSignals.push_back(initValue);
	}

	return pair<vector<float>, vector<float>>(peakSignals, allSignals);
}

float min_inRange(TH1F *hWave, float t1, float t2)
{
	TF1 *f1 = new TF1("f1", "pol0", t1, t2);
	double r1 = 0;
	double r2 = 0;

	hWave->GetXaxis()->SetRange(t1 / SP, t2 / SP); //window
	r1 = hWave->GetMinimumBin() * SP - 0.5;
	r2 = r1 + 1;

	hWave->Fit("f1", "QN", "", r1, r2);
	float max = f1->GetParameter(0);
	hWave->GetXaxis()->SetRange(1, 1024);

	return max;
}

/**
 *  The goal of this class is to determine the perfect integration window for each run and each channel.
 *  For this the integration interval is calculated by the sum Histograms -> Needs SumHistogram attached to the rootfile -> "Full Runmode"
 * 	The uncertainty is measured from the difference between the mean in the Percentages taken from sumhistogram and each waveform
 **/

string outputFolder = "./IntegrationWindowPlots/";
int smoothDefault = 1500;

int peakLeftThreshold = 100;
int peakRightThreshold = 280;
bool singlePrinter = true;
bool hardMode =true; //There are runs with extremely low light intensity (45,46,47 e.g) without a sharp first APO peak -> get some integration window for them
float labelTextSize = 0.02;
int lineSize = 2;

int main(int argc, char *argv[])
{

	cout << "HALLO COMPUTER 2" << endl;
	bool printPercentage = true;
	gStyle->SetOptStat(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
	gStyle->SetLineStyleString(11, "5 5");

	gErrorIgnoreLevel = kError;

	if (singlePrinter)
	{
		labelTextSize = 0.03;
		lineSize = 4;
	}

	string runName = (string)argv[1];
	int runNr = atoi(argv[2]);
	string inDir = string(argv[3]);
	printf("Runname IW Analysis %s \n", runName.c_str());
	printf("runNr IW Analysis %d \n", runNr);
	printf("In File IW Analysis %s \n", inDir.c_str());

	//Output Plot Folder
	cout << 1 << endl;
	string outDir = Form((outputFolder + "%s/").c_str(), runName.c_str());
	cout << outDir << endl;
	cout << outputFolder << endl;
	if (!fs::is_directory(outputFolder) || !fs::exists(outputFolder))
	{
		fs::create_directory(outputFolder);
		cout << 2 << endl;
	}
	if (!fs::is_directory(outDir) || !fs::exists(outDir))
	{
		fs::create_directory(outDir);
		cout << 3 << endl;
	}
	cout << "Doing: " << inDir << endl;
	//Open Root File
	TFile *file = new TFile(inDir.c_str());
	if (file->IsZombie())
	{
		cout << "PROBLEM with the initialization of the output ROOT ntuple "
			 << inDir << ": check that the path is correct!!!"
			 << endl;

		exit(-1);
	}
	TTree *tree;
	file->GetObject("T", tree);

	//Create 2 Canvas -> 1 for sumHistograms, 1 for Percentage Plots
	TCanvas *sumCanvas = new TCanvas("sumCanvas", "Sum Histogram", 1000, 1000);
	TCanvas *percentageCanvas = new TCanvas("percentageCanvas", "Percentages", 1000, 1000);

	int channelNumber = 8;
	file->GetObject("T", tree);
	//tree->SetBranchAddress("nCh", &channelNumber);
	//tree->GetEntry(1);

	int plotGrid = ceil(sqrt(channelNumber));

	if (!singlePrinter)
	{
		sumCanvas->Divide(plotGrid, plotGrid);
		percentageCanvas->Divide(plotGrid, plotGrid);
	}
	else
	{
		sumCanvas->Divide(1, 1);
		percentageCanvas->Divide(1, 1);
	}

	sumCanvas->SetLeftMargin(0.15);
	FILE *file_list;
	string list_filename = outputFolder + "/IntegrationWindows.txt";
	file_list = fopen(list_filename.c_str(), "a");

	FILE *file_listP;
	string list_filenameP = outputFolder + "/PercentageValues.txt";
	file_listP = fopen(list_filenameP.c_str(), "a");

	int integrationLeftOffset = 20;
	if (runName.find("calib") != string::npos)
	{
		integrationLeftOffset = 10;
	}

	fprintf(file_listP, "PV_%s={", string(runName).c_str());
	fprintf(file_list, "IW_%s={", string(runName).c_str());

	for (int i = 0; i < channelNumber; i++)
	{
		//OscillationAnalysis ---------------------------------------------------------
		cout << " " << endl;

		if (!singlePrinter)
		{
			sumCanvas->cd(i + 1);
		}
		else
		{
			sumCanvas->cd(1);
		}
		sumCanvas->SetGrid();
		//	gPad->SetGridx();
		//	gPad->SetGridy();

		TH1F *sumHist;
		TH1F *sumHistCalculation; //Without smoothing

		TLegend *sumHistLegend = new TLegend(0.60, 0.60, 0.99, 1);
		sumHistLegend->SetTextFont(62);

		sumHistLegend->SetHeader(Form("SumHistogram Channel: %d", i), "c");
		sumHistLegend->SetTextSize(0.02);

		string name = Form("hChSum_%d", i);
		file->GetObject(name.c_str(), sumHist);
		file->GetObject(name.c_str(), sumHistCalculation);

		cout << "Entries: " << sumHist->GetEntries() / 1024.0 << endl;
		int smooth = smoothDefault;
		switch (i)
		{
		case 3:
			//	smooth = 3000;
			break;

		default:
			smooth = 0;
			break;
		}
		smooth = smoothDefault;

		sumHist->Smooth(smooth);

		sumHist->Draw("hist");
		sumHist->SetLineColor(1);

		sumHistLegend->AddEntry(sumHist, Form("%s", runName.c_str()), "l");
		sumHistLegend->SetTextFont(42);

		TAxis *yaxis = sumHist->GetYaxis();
		TAxis *xaxis = sumHist->GetXaxis();
		yaxis->SetLabelSize(labelTextSize);
		xaxis->SetLabelSize(labelTextSize);

		sumCanvas->Update();
		float t_amp = t_max_inRange(sumHist, 0.0, 320.0);
		float minY = gPad->GetUymin();
		float maxY = gPad->GetUymax();

		sumHist->SetAxisRange(30, 290); //IMPORTANT -> EXCLUDE OSZILLATIONS AT END OR PF WONT FIND MANY PEAKS

		Int_t npeaks = 10; //number of peaks to fit
		cout << "npeaks: " << npeaks << endl;

		TSpectrum *s = new TSpectrum(npeaks, 15); // TSpectrum class finds candidates peaks
		Int_t nfound = s->Search(sumHist, 0, "goff", 0.0001);
		sumHist->SetAxisRange(0, 320);

		Double_t par[100]; //Array of fitting parameters
		npeaks = 0;
		int p = 0;

		//Find Peaks in sumHistogram
		Double_t *xpeaks = s->GetPositionX(); //array with X-positions of the centroids found by TSpectrum

		//GetParameters of found peaks
		for (p = 0; p < nfound; p++)
		{
			Double_t xp = xpeaks[p];
			Int_t bin = sumHist->GetXaxis()->FindBin(xp);
			Double_t yp = sumHist->GetBinContent(bin);

			if (xp < (t_amp - 100)) //Too far left
			{
				continue;
			}
			if (xp > (peakRightThreshold)) //Too far right, skip weird
			{
				continue;
			}
			if (xp < (peakLeftThreshold)) //Too far right, skip weird
			{
				continue;
			}
			if (p > 0)
			{
				if (abs(xp - xpeaks[p - 1]) < 20)
				{
					//Double fit of a single peak-> Skip
					continue;
				}
			}
			cout << "PEAKS: " << i << "  " << xp << "  " << yp << endl;

			par[3 * npeaks + 2] = yp; //height
			par[3 * npeaks + 3] = xp; //centroid
			par[3 * npeaks + 4] = 50; //unused
			npeaks++;
		}

		printf("IntegrationWindow: : %s, Channel: %d, %d useful peaks\n", runName.c_str(), i, npeaks);

		//Gauß Fits on the peaks
		TF1 *peak_single[npeaks];
		Double_t par_single[3 * npeaks];
		Double_t chi2ndf[npeaks];
		Double_t meanError[npeaks]; //Parametric Error

		for (int i = 0; i < npeaks; ++i)
		{
			double pos_peak = par[3 + 3 * i];
			float range = 5; //Width of gauß Fit
			peak_single[i] = new TF1("peak", "gaus", pos_peak - range, pos_peak + range);
			peak_single[i]->SetLineStyle(11);
			peak_single[i]->SetLineColor(kRed);
			if (singlePrinter)
				peak_single[i]->SetLineWidth(0);
			sumHist->Fit("peak", "RQ+");
			if (!singlePrinter)
				peak_single[i]->Draw("same");
			chi2ndf[i] = peak_single[i]->GetChisquare() / peak_single[i]->GetNDF();
			meanError[i] = peak_single[i]->GetParError(1);
			peak_single[i]->GetParameters(&par_single[3 * i]);
		}

		int j = 0;
		float periods[npeaks - 1];
		float periodErrors[npeaks - 1]; //From the parametric mean error -> unused
		float meanPeriod = 0;
		float meanPeriodError = 0; //Unused
		float meanPeriodStdD = 0;  //Unused
		float means[npeaks];	   //Gauß Fit Means

		for (j = 1; j <= npeaks; j++)
		{
			means[j - 1] = par_single[1 + (j - 1) * 3];

			if (j < npeaks)
			{
				float period = par_single[1 + (j)*3] - par_single[1 + (j - 1) * 3];
				cout << "PERIOD FOR: " << i << " " << j << " PERIODE: " << period << "  " << npeaks << " In1: " << par_single[1 + (j)*3] << "  In2: " << par_single[1 + (j - 1) * 3] << endl;

				if (hardMode)
				{
					if (period > 70)
					{
						//One peak is missing, so 1 is skipped -> Wrong period
						period = period / 2.0;
					}
				}

				float periodError = sqrt(pow(meanError[j], 2) + pow(meanError[j - 1], 2));
				periods[j - 1] = period;
				periodErrors[j - 1] = periodError;
			}
		}
		sumHistLegend->AddEntry(peak_single[0], Form("Amplitude %1.2f +- %1.2f", means[0], meanError[0]), "l");

		//Determine order and integration Windows
		float minSearchLeft = 1000;
		float minSearchRight = 1000;
		float meanPeriodCounter = 0;
		for (j = 0; j < npeaks; j++)
		{
			if (periods[j] < 200 && periods[j] > 20)
			{ //Prevent bugs

				meanPeriod = meanPeriod + periods[j];
				meanPeriodError = meanPeriodError + pow(periodErrors[j], 2) + pow(periodErrors[j - 1], 2);
				meanPeriodCounter++;
			}
			if (means[j] < minSearchLeft)
			{
				minSearchLeft = means[j];
			}
			else if (means[j] < minSearchRight)
			{
				minSearchRight = means[j];
			}
			else
			{
			}
		}

		meanPeriod = meanPeriod / (meanPeriodCounter);
		cout << "CHANNEL: " << i << " MEAN PERIOD: " << meanPeriod << "  Counter: " << meanPeriodCounter << endl;
		//meanPeriod=periods[0];
		for (j = 0; j < npeaks; j++)
		{
			meanPeriodStdD = meanPeriodStdD + pow(meanPeriod - periods[j], 2);
		}
		meanPeriodStdD = sqrt(meanPeriodStdD / (npeaks - 1));

		float statisticalMeanPeriodError = meanPeriodStdD / sqrt((npeaks - 1));		  //Makes no sense, cause only 3 values
		float systematicMeanPeriodError = sqrt(meanPeriodError) / sqrt((npeaks - 1)); //Only from parametric erros -> not the real error
		//meanPeriodError = sqrt(pow(statisticalMeanPeriodError, 2) + pow(systematicMeanPeriodError, 2));
		meanPeriodError = sqrt(pow(systematicMeanPeriodError, 2));

		//	sumHistLegend->AddEntry((TObject *)0, Form("Mean Period: %1.2f +- %1.2f", meanPeriod, meanPeriodError), "");

		//Find first Minimum -> Search Between 0,1 Maxima
		float sumPercentage = 0;
		if (npeaks > 1)
		{

			float fitLeft = par_single[1 + 3] - 30;
			float fitRight = par_single[1 + 3] - 5;
			float searchLeft = fitLeft;
			float searchRight = fitRight;

			if (hardMode)
			{
				switch (runNr)
				{
				case 45:
					fitRight = 120;
					searchLeft = 150;
					searchRight = 170;
					break;
				case 46:
					fitRight = 180;
					searchLeft = 160;
					searchRight = 170;
					break;
				case 47:
					fitRight = 180;
					searchLeft = 160;
					searchRight = 170;
					break;
				default:
					break;
				}
			}

			TF1 *putc = new TF1("inversepeak", "pol4", fitLeft, fitRight);
			putc->SetLineColor(8);
			putc->SetLineStyle(11);
			putc->SetParameter(0, 10);
			putc->SetParameter(1, 155);
			putc->SetParameter(2, 5);
			sumHist->Fit("inversepeak", "RQ+");
			if (!singlePrinter)
				putc->Draw("same");

			float minInRange = sumHist->GetFunction("inversepeak")->GetMinimumX(searchLeft, searchRight);
			cout << "DEBUG: " << minInRange << "  XX " << fitLeft << "  y " << fitRight << endl;

			float minInRangeError = 320.0 / (sumHist->GetNbinsX());
			TLine *minLine = new TLine(minInRange, minY, minInRange, maxY);
			minLine->SetLineColor(3);
			minLine->SetLineWidth(lineSize);

			float entireSignalRight = 3 * meanPeriod + minInRange;

			TLine *entireSignalRightLine = new TLine(entireSignalRight, minY, entireSignalRight, maxY);
			entireSignalRightLine->SetLineColor(9);
			entireSignalRightLine->SetLineWidth(lineSize);
			entireSignalRightLine->Draw();

			TLine *leftLine = new TLine(means[0] - integrationLeftOffset, minY, means[0] - integrationLeftOffset, maxY);
			leftLine->SetLineColor(2);
			leftLine->SetLineWidth(lineSize);

			float allSignalError = sqrt(pow(3 * meanPeriodError, 2) + pow(meanError[0], 2));

			//	sumHistLegend->AddEntry((TObject *)0, Form("Interval [%1.2f,%1.2f]", -20.0, minInRange - means[0]), "");

			//CALCULATIONS
			float BL_output[4];
			BL_fit(sumHistCalculation, BL_output, 30.0, 70.0);
			float BL_shift = BL_output[0];

			TLine *baselineUsed = new TLine(30, BL_shift, 70, BL_shift);
			baselineUsed->SetLineColor(9);
			baselineUsed->SetLineWidth(lineSize);

			float integralPeak = IntegralHist(sumHistCalculation, means[0] - integrationLeftOffset, minInRange, BL_shift);
			float integralPeakAll = IntegralHist(sumHistCalculation, means[0] - integrationLeftOffset, entireSignalRight, BL_shift);
			sumPercentage = integralPeak / integralPeakAll;

			//	sumHistLegend->AddEntry((TObject *)0, Form("Interval w Errors [%1.2f,%1.2f]", -(20 + meanError[0]), (minInRange - means[0]) + minInRangeError), "");
			sumHistLegend->AddEntry(peak_single[0], Form("Mean Period: %1.2f", meanPeriod), "l");
			sumHistLegend->AddEntry(baselineUsed, Form("Baseline: %1.2f", BL_shift), "l");
			sumHistLegend->AddEntry(leftLine, Form("Left %1.2f +- %1.2f", means[0] - integrationLeftOffset, meanError[0]), "l");
			sumHistLegend->AddEntry(minLine, Form("Signal %1.2f +- %1.2f", minInRange - means[0], minInRangeError), "l");
			sumHistLegend->AddEntry(entireSignalRightLine, Form("All %1.2f +- %1.2f", entireSignalRight - means[0], allSignalError), "l");

			sumHist->SetMarkerSize(0.4);
			sumHist->SetMarkerStyle(21);

			TH1 *dummyforMarker = new TH1I("h1", "h1 title", 0, 0.0, 0.0);
			dummyforMarker->SetMarkerSize(0.4);
			dummyforMarker->SetMarkerStyle(21);
			dummyforMarker->SetMarkerColorAlpha(30, 0.9);

			sumHist->GetXaxis()->SetRangeUser(means[0] - integrationLeftOffset, minInRange);

			sumHist->SetFillColorAlpha(kTeal + 2, 0.9);

			sumHist->SetFillStyle(1001);
			sumHist->DrawClone("same hist");
			sumHistLegend->AddEntry(dummyforMarker, Form("Integral: Peak: %1.1e", integralPeak), "P");
			sumHistLegend->AddEntry(dummyforMarker, Form("IWindow Size: %1.2f", minInRange - (means[0] - integrationLeftOffset)), "l");

			sumHist->GetXaxis()->SetRangeUser(means[0] - integrationLeftOffset, entireSignalRight);

			//sumHist->SetFillColorAlpha(kAzure-8, 0.35);
			sumHist->SetMarkerColorAlpha(44, 0.35);

			sumHist->SetFillStyle(3342);
			sumHist->DrawClone("same hist");
			sumHistLegend->AddEntry(sumHist, Form("Integral: All: %1.1e", integralPeakAll), "P");

			sumHistLegend->AddEntry((TObject *)0, Form("Percentage (Signal/All): %1.4f", sumPercentage), "");

			//Reset -> Needs to be done
			sumHist->GetXaxis()->SetRangeUser(0, 320);
			sumHist->GetXaxis()->SetRangeUser(0, 320);

			sumHist->SetFillColorAlpha(0, 0.0);

			for (int i = 0; i < npeaks; ++i)
			{
				peak_single[i]->Draw("same");
				//Redraw over Integral
			}
			sumHist->SetFillStyle(3001);
			minLine->Draw();
			leftLine->Draw();
			if (!singlePrinter)
				baselineUsed->Draw();
			if (!singlePrinter)
				putc->Draw("same");

			fprintf(file_list, "%f/%f,", minInRange - means[0], entireSignalRight - means[0]);
		}
		else
		{
			fprintf(file_list, "%f/%f,", 0.0, 0.0);
		}

		//if(!singlePrinter)
		sumHistLegend->Draw();

		/***
 *      _    _                     _        _       _         
 *     | |  | |                   | |      (_)     | |        
 *     | |  | |_ __   ___ ___ _ __| |_ __ _ _ _ __ | |_ _   _ 
 *     | |  | | '_ \ / __/ _ \ '__| __/ _` | | '_ \| __| | | |
 *     | |__| | | | | (_|  __/ |  | || (_| | | | | | |_| |_| |
 *      \____/|_| |_|\___\___|_|   \__\__,_|_|_| |_|\__|\__, |
 *                                                       __/ |
 *                                                      |___/ 
 */

		printf("Uncertainty: %s, Channel: %d\n", runName.c_str(), i);
		if (!singlePrinter)
		{
			percentageCanvas->cd(i + 1);
		}
		else
		{
			percentageCanvas->cd(1);
		}

		float xmin = -2.0; //EXTREMELY IMPORTANT -> include all values in the histogram, otherwise the CF is biased
		float xmax = +2;
		int nBins = 200;
		//gPad->SetGridx();
		//gPad->SetGridy();

		TH1F *h = new TH1F("h", "Percentage", nBins, xmin, xmax);
		h->SetTitle("");
		h->SetLineColorAlpha(kBlack, 1);
		h->SetLineWidth(3);

		h->SetMarkerStyle(7);
		h->SetMarkerColorAlpha(kBlack, 1);

		TString cut("");
		tree->Draw(Form("IntegralDifference[%d]>>h", i), cut);
		float entries = h->GetMean();
		if (entries == 1.000f)
			printPercentage = false;
		percentageCanvas->Update();
		minY = gPad->GetUymin();
		maxY = gPad->GetUymax();
		TLegend *h_leg = new TLegend(0.54, 0.75, 0.90, 0.9);
		h_leg->SetTextFont(62);
		if (!singlePrinter)
			h_leg->SetHeader(Form("Percentages Channel: %d", i), "c");
		h_leg->SetTextSize(labelTextSize);
		if (!singlePrinter)
			h_leg->AddEntry(h, Form("%s", runName.c_str()), "l");
		h_leg->SetTextFont(42);
		//h_leg->AddEntry(h, Form("entries: %1.0lf  std: %1.2f", h->GetEffectiveEntries(), h->GetStdDev()), "l");

		TAxis *yaxisP = h->GetYaxis();
		TAxis *xaxisP = h->GetXaxis();
		yaxisP->SetLabelSize(labelTextSize);
		yaxisP->SetTitle("counts");
		xaxisP->SetLabelSize(labelTextSize);

		xaxisP->SetTitle(Form("f_{W}"));
		xaxisP->SetTitleSize(labelTextSize);
		yaxisP->SetTitleSize(labelTextSize);

		double meanP = h->GetMean();
		double meanErrorP = h->GetMeanError();

		TLine *meanLine = new TLine(meanP, minY, meanP, maxY);
		meanLine->SetLineColor(2);
		meanLine->SetLineWidth(lineSize);

		meanLine->Draw();

		TLine *sumPLine = new TLine(sumPercentage, minY, sumPercentage, maxY);
		sumPLine->SetLineColor(4);
		sumPLine->SetLineWidth(lineSize);
		sumPLine->Draw();

		double pError = sqrt(pow(abs(meanP - sumPercentage), 2) + pow(meanErrorP, 2));

		TLine *meanLineError = new TLine(sumPercentage, minY, sumPercentage, maxY);
		meanLineError->SetLineWidth(pError * 1125);
		meanLineError->SetLineColorAlpha(4, 0.14);
		//	meanLineError->Draw();

		//double mean=gauss->GetParameter(1);
		h_leg->AddEntry(meanLine, Form("#bar{f}_{W}: %1.3lf +- %1.3lf", meanP, meanErrorP), "l");
		h_leg->AddEntry(sumPLine, Form("f_{W}^{CSW}: %1.3lf", sumPercentage), "l");

		TLine *dummyP = new TLine(0, 0, 0, 0);
		dummyP->SetLineColorAlpha(4, 0.14);
	//	h_leg->AddEntry(dummyP, Form("Err D.  : %1.3lf  ", pError), "l");
		
		h->GetXaxis()->SetRangeUser(0.2,1.6);
		h_leg->Draw();




		if (printPercentage)
			fprintf(file_listP, "%f/%f,", meanP, pError);

		if (singlePrinter)
		{
			sumCanvas->Print((outDir + runName + "_" + to_string(i) + "_IntegrationWindow.pdf").c_str());
			if (printPercentage)
				percentageCanvas->Print((outDir + runName + "_" + to_string(i) + "_PercentageDistribution.pdf").c_str());
		}
	}

	fprintf(file_list, "}\n");
	fclose(file_list);
	if (printPercentage)
		fprintf(file_listP, "}\n");
	fclose(file_listP);

	sumCanvas->Print((outDir + runName + "_IntegrationWindow.pdf").c_str());
	if (printPercentage)
		percentageCanvas->Print((outDir + runName + "_PercentageDistribution.pdf").c_str());

	/***
 *       _____                         _   _              __           _             
 *      / ____|                       | | (_)            / _|         | |            
 *     | |     ___  _ __ _ __ ___  ___| |_ _  ___  _ __ | |_ __ _  ___| |_ ___  _ __ 
 *     | |    / _ \| '__| '__/ _ \/ __| __| |/ _ \| '_ \|  _/ _` |/ __| __/ _ \| '__|
 *     | |___| (_) | |  | | |  __/ (__| |_| | (_) | | | | || (_| | (__| || (_) | |   
 *      \_____\___/|_|  |_|  \___|\___|\__|_|\___/|_| |_|_| \__,_|\___|\__\___/|_|   
 *                                                                                   
 *                                                                                   
 */

	string runNameOfCalibration = "dummy"; //"7_calib_vb58_tune8700_pcbd";
	//Find Content of calibration Run inside text file
	ifstream data_store(list_filenameP);
	string lineWithDataCalib;
	string lineWithDataMeasurement;

	for (string line; getline(data_store, line);)
	{
		if (line.find(runNameOfCalibration) != string::npos)
		{
			lineWithDataCalib = line; //Use the last line on the txt matching the condition
		}
		if (line.find(runName) != string::npos)
		{
			lineWithDataMeasurement = line; //Use the last line on the txt matching the condition
		}
	}
	cout << "Using Calibration for CF: " << lineWithDataCalib << endl;
	cout << "Using Measurement for CF: " << lineWithDataMeasurement << endl;

	cout << "Parsing and Calculating..." << endl;

	pair<vector<float>, vector<float>> calibCombined = readPairs(lineWithDataCalib, 0.0);
	pair<vector<float>, vector<float>> measurementCombined = readPairs(lineWithDataMeasurement, 0.0);

	vector<float> calibPValues = calibCombined.first;
	vector<float> calibPValuesError = calibCombined.second;
	vector<float> measurementPValues = measurementCombined.first;
	vector<float> measurementPValuesError = measurementCombined.second;

	FILE *file_listCF;
	FILE *file_listCFF;

	FILE *file_listfW;

	string list_filenameCF = outputFolder + "/CorrectionValues.txt";
	string list_filenameFW = outputFolder + "/FWValues.txt";
	string list_filenameCFF = outputFolder + "/CorrectionValuesFormated.txt";

	file_listCFF = fopen(list_filenameCFF.c_str(), "a");
	file_listCF = fopen(list_filenameCF.c_str(), "a");
	file_listfW = fopen(list_filenameFW.c_str(), "a");

	if (printPercentage)
		fprintf(file_listCF, "CF_%s={", string(runName).c_str());
	fprintf(file_listfW, "%d,", runNr);
	if (printPercentage)
		fprintf(file_listCFF, "%d,", runNr);

	for (int i = 0; i < channelNumber; i++)
	{
		float calibPVal = calibPValues[i];
		float calibPValError = calibPValuesError[i];
		float measurementPVal = measurementPValues[i];
		float measurementPValError = measurementPValuesError[i];

		float correctionFactor = calibPVal / measurementPVal;
		float correctionFactorError = sqrt(pow((calibPValError / measurementPVal), 2) + pow(measurementPValError * calibPVal / (pow(measurementPVal, 2)), 2));

		cout << "Correction Factor: Channel: " << i << " :: " << correctionFactor << " +- " << correctionFactorError << endl;
		if (printPercentage)
			fprintf(file_listCF, "%f/%f,", correctionFactor, correctionFactorError);
		if (printPercentage)
			fprintf(file_listCFF, "%f/%f,", correctionFactor, correctionFactorError);
		fprintf(file_listfW, "%f/%f,", measurementPVal, measurementPValError);
	}
	if (printPercentage)
		fprintf(file_listCF, "}\n");
	if (printPercentage)
		fprintf(file_listCFF, "}\n");
	fprintf(file_listfW, "\n");

	if (printPercentage)
		fclose(file_listCFF);
	if (printPercentage)
		fclose(file_listCF);
	fclose(file_listfW);

	return 0;
}