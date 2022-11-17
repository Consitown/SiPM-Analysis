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
#include "TEfficiency.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include <algorithm>

namespace fs = std::experimental::filesystem;
using namespace std;
struct stat info;
float SP = 0.3125; // ns per bin

string outputFolder = "./EfficiencyPlots";
string thresholdsFile = "./Thresholds_Integral.txt";
string thresholdsFileErrorP = "./Thresholds_IntegralErrorP.txt";
string thresholdsFileErrorM = "./Thresholds_IntegralErrorM.txt";

string runNameDC = "331_dc_vb58_pcbd";



string extractValues(string str)
{
	unsigned first = str.find("{");
	unsigned last = str.find("}");
	string strNew = str.substr(first + 1, last - first - 1);
	return strNew;
}
float stringToFloat(string text)
{
	if (text.length() == 0)
		return -10.0;
	return stof(text);
}
pair<vector<float>, vector<float>> readPairs(string line, double initValue)
{
	cout << "hello world";
	
	vector<float> p1;
	vector<float> p2;

	p1.clear();
	p2.clear();

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
		float p = stringToFloat(peak);
		if (p >= 0)
		{
			p1.push_back(p);
			p2.push_back(stringToFloat(all));
		}

		s.erase(0, pos + columnDelimiter.length());
	}
	int posOfDelimiter = s.find("/");
	string peak = s.substr(0, posOfDelimiter);
	string all = s.substr(posOfDelimiter + 1, s.length() - 1);

	float p = stringToFloat(peak);
	if (p >= 0)
	{
		p1.push_back(stringToFloat(peak));
		p2.push_back(stringToFloat(all));
	}

	for (std::size_t i = p1.size(); i < 9; i++)
	{
		p1.push_back(initValue);
		p2.push_back(initValue);
	}

	return pair<vector<float>, vector<float>>(p1, p2);
}
string vectorToString(vector<float> vec)
{
	std::ostringstream vts;

	if (!vec.empty())
	{
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec.begin(), vec.end() - 1,
				  std::ostream_iterator<float>(vts, ", "));

		// Now add the last element with no delimiter
		vts << vec.back();
	}
	return vts.str();
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

int main(int argc, char *argv[])
{
	gStyle->SetOptStat(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
	gStyle->SetLineStyleString(11, "5 5");
	gErrorIgnoreLevel = kError;

	string runName = (string)argv[1];
	int runNr = atoi(argv[2]);
	string inDir = string(argv[3]);

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

	gStyle->SetOptTitle(0);
	int channelNumber = 16;
	int entryNumber = 16;

	file->GetObject("T", tree);
	tree->SetBranchAddress("nCh", &channelNumber);
	tree->GetEntry(1);

	int plotGrid = ceil(sqrt(channelNumber));
	bool zoomMode=false;
	int numberToScaleUp=10010;
	float dcBinSize = 0.1;
	int maxX = 4500;
	int min = -100;
	int numberOfNPE=maxX-min;
	int bins = numberOfNPE*(1.0/dcBinSize); //0.1 Photoelectron
	bins=450;

	if(zoomMode){
		bins=50;
		maxX=60;
		min=-10;
	}
	float threshold = -1;
	float thresholdErrorP = -1;
	float thresholdErrorM = -1;


	switch (runNr)
	{
	case 86:
		numberToScaleUp=30030;
		break;
	case 87:
		numberToScaleUp=11011;
		break;
	
	default:
		numberToScaleUp=10010;
		break;
	}







	std::ifstream inFile(thresholdsFile);
	int numberOfLines = std::count(std::istreambuf_iterator<char>(inFile),
								   std::istreambuf_iterator<char>(), '\n');

	//find threshold
	ifstream data_store(thresholdsFile);
	ifstream data_storeErrorP(thresholdsFileErrorP);
	ifstream data_storeErrorM(thresholdsFileErrorM);

	TCanvas *effCanvas = new TCanvas("effCanvas", "Sum Histogram", 1000, 650);
	effCanvas->SetGrid();

	bool foundThreshold=false;
	string lineWithDataMeasurement;
	for (string lineWithDataMeasurement_; getline(data_store, lineWithDataMeasurement_);)
	{
		
		string d=split(lineWithDataMeasurement_, "_")[2];

		string r=split(d,"=")[0].erase(0,2);

	/*	if (d.find("dc") != std::string::npos) {
				continue;
			}*/
			int runNummerOfThreshold = stoi(r);
			if (runNr != runNummerOfThreshold)
			{
				continue;
			}
		foundThreshold=true;
		lineWithDataMeasurement=lineWithDataMeasurement_;
	}
		
		cout<<""<<endl;
		cout<<"\nCURRENT RUN: "<<runName<<endl;
		cout<<"FOUND THRESHOLD: "<<foundThreshold<<endl;
		cout<<"LINE: "<<lineWithDataMeasurement<<endl;

		if(!foundThreshold)return 0;

			//Output Plot Folder
		string outDir = Form((outputFolder + "/%s/").c_str(), runName.c_str());
		if (!fs::is_directory(outputFolder) || !fs::exists(outputFolder))
		{
			fs::create_directory(outputFolder);
		}
		if (!fs::is_directory(outDir) || !fs::exists(outDir))
		{
			fs::create_directory(outDir);
		}
	
		string thresholdName = lineWithDataMeasurement.substr(0, lineWithDataMeasurement.find("="));
		string lineWithDataMeasurementErrorP;
		string lineWithDataMeasurementErrorM;

		for (string lineWithDataMeasurementErrorMTemp; getline(data_storeErrorM, lineWithDataMeasurementErrorMTemp);)
		{

			if (lineWithDataMeasurementErrorMTemp.substr(0, lineWithDataMeasurementErrorMTemp.find("=")) == thresholdName)
			{
				lineWithDataMeasurementErrorM = lineWithDataMeasurementErrorMTemp;
			}
		}
		for (string lineWithDataMeasurementErrorPTemp; getline(data_storeErrorP, lineWithDataMeasurementErrorPTemp);)
		{

			if (lineWithDataMeasurementErrorPTemp.substr(0, lineWithDataMeasurementErrorPTemp.find("=")) == thresholdName)
			{
				lineWithDataMeasurementErrorP = lineWithDataMeasurementErrorPTemp;
			}
		}

		pair<vector<float>, vector<float>> combined = readPairs(lineWithDataMeasurement, 0.0);
		vector<float> npeCuts = combined.first;
		vector<float> percentages = combined.second;

		pair<vector<float>, vector<float>> combinedErrorP = readPairs(lineWithDataMeasurementErrorP, 0.0);
		vector<float> npeCutsErrorP = combinedErrorP.first;
		vector<float> percentagesErrorP = combinedErrorP.second;

		pair<vector<float>, vector<float>> combinedErrorM = readPairs(lineWithDataMeasurementErrorM, 0.0);
		vector<float> npeCutsErrorM = combinedErrorM.first;
		vector<float> percentagesErrorM = combinedErrorM.second;

		float sumCut = npeCuts.back();
		float sumCutErrorP = npeCutsErrorP.back();
		float sumCutErrorM = npeCutsErrorM.back();

		float sumPercentage = percentages.back() * 100;
		float sumPercentageErrorP = percentagesErrorP.back() * 100;
		float sumPercentageErrorM = percentagesErrorM.back() * 100;

	
		sumCut = sumCut - dcBinSize / 2;
		threshold = sumCut;
		sumCutErrorP = sumCutErrorP - dcBinSize / 2;
		thresholdErrorP = sumCutErrorP;
		sumCutErrorM = sumCutErrorM - dcBinSize / 2;
		thresholdErrorM = sumCutErrorM;

		printf("Sum Cut: %1.2f (%1.2f,%1.2f) ::: Percentages: %1.2f (%1.2f,%1.2f) \n", sumCut, sumCutErrorP, sumCutErrorM, sumPercentage, sumPercentageErrorP, sumPercentageErrorM);

		TH1D *allHist = new TH1D("allHist", "", bins, min, maxX);
		TH1D *allHistErrorP = new TH1D("allHistErrorP", "", bins, min, maxX);
		TH1D *allHistErrorM = new TH1D("allHistErrorM", "", bins, min, maxX);

		allHist->SetFillColorAlpha(kBlue, 0.4);
		allHist->SetLineColorAlpha(1, 0);

		TH1D *cuttedHist = new TH1D("cuttedHist", "", bins, min, maxX);
		TH1D *cuttedHistErrorP = new TH1D("cuttedHistErrorP", "", bins, min, maxX);
		TH1D *cuttedHistErrorM = new TH1D("cuttedHistErrorM", "", bins, min, maxX);

		//Set threshold according to DC measurement results
		TString cutNum1(Form("(chargeChannelSumWOM[3])>=%f", threshold));
		TString cutNum2(Form("(chargeChannelSumWOMErrorP[3])>=%f", thresholdErrorP));
		TString cutNum3(Form("(chargeChannelSumWOMErrorM[3])>=%f", thresholdErrorM));

		//Fill
		tree->Draw("chargeChannelSumWOM[3]>>allHist", "");
		tree->Draw("chargeChannelSumWOMErrorP[3]>>allHistErrorP", "");
		tree->Draw("chargeChannelSumWOMErrorM[3]>>allHistErrorM", "");

		tree->Draw("chargeChannelSumWOM[3]>>cuttedHist", cutNum1);
		tree->Draw("chargeChannelSumWOMErrorP[3]>>cuttedHistErrorP", cutNum2);
		tree->Draw("chargeChannelSumWOMErrorM[3]>>cuttedHistErrorM", cutNum3);


		int entry_temp=allHist->GetEntries();
		int entryP_temp=allHistErrorP->GetEntries();
		int entryM_temp=allHistErrorM->GetEntries();

		int entryC_temp=cuttedHist->GetEntries();
		int entryCP_temp=cuttedHistErrorP->GetEntries();
		int entryCM_temp=cuttedHistErrorM->GetEntries();
		//There are maybe entries skipped due to exceeding events, since they need to be counted as detected, the numbers need to be adjusted
		
		cout<<"Raw Entries All: "<<entry_temp<<"/"<<entryC_temp<<endl;
		cout<<"Raw Entries P: "<<entryP_temp<<"/"<<entryCP_temp<<endl;
		cout<<"Raw Entries M: "<<entryM_temp<<"/"<<entryCM_temp<<endl;




		int skippedEntries=numberToScaleUp-entry_temp;
		int skippedEntriesP=numberToScaleUp-entryP_temp;
		int skippedEntriesM=numberToScaleUp-entryM_temp;

		if(skippedEntries<0){
			cout<<"Scale UP ERROR, histogram has: "<<entry_temp<<" entries, but numberToScaleUp is: "<<numberToScaleUp<<endl;
			assert(0);
		}
		cout<<"Skipped: All: "<<skippedEntries<<"  P: "<<skippedEntriesP<<" M: "<< skippedEntriesM<<endl;

		if(!((skippedEntries == skippedEntriesM) && (skippedEntries== skippedEntriesP))){
			cout<<"The entries in the Error distributions do not match. Something went wrong on analysis."<<endl;
			
			assert(0);
		}

		//Adjusting the numbers -> Adding the skipped ones
		float entry=entry_temp+skippedEntries;
		float entryP=entryP_temp+skippedEntries;
		float entryM=entryM_temp+skippedEntries;

		float entryC=entryC_temp+skippedEntries;
		float entryCP=entryCP_temp+skippedEntries;
		float entryCM=entryCM_temp+skippedEntries;


		float efficiency =  (entryC/entry)* 100;
		float efficiencyP = (entryCP/entryP)* 100;
		float efficiencyM = (entryCM/entryM)* 100;




		float deviationP = abs(efficiencyP - efficiency);
		float deviationM = abs(efficiencyM - efficiency);
		//float maxDeviation = max(deviationM, deviationP);

		double upperStatErr = abs(efficiency - TEfficiency().ClopperPearson(entry, entryC, 0.682689492137, true) * 100);
		double lowerStatErr = abs(efficiency - TEfficiency().ClopperPearson(entry, entryC, 0.682689492137, false) * 100);

		float combinedUpperError =sqrt(pow(upperStatErr, 2) + pow(deviationP, 2));
		
		 //std::min(sqrt(pow(upperStatErr, 2) + pow(deviationP, 2)),(double)(100-efficiency)); //max 100%, not needed


		float combinedLowerError = sqrt(pow(lowerStatErr, 2) + pow(deviationM, 2));
		cout << "Run: " << runName << "  Unscaled All: " << entry_temp << "  Cut: " << entryC_temp << endl;
		cout << "Run: " << runName << "  EffErrP: " << efficiencyP << "  Unscaled All: " << entryP_temp << "  Cut: " << entryCP_temp << endl;
		cout << "Run: " << runName << "  EffErrM: " << efficiencyM << "  Unscaled All: " << entryM_temp << "  Cut: " << entryCM_temp << endl;
		//cout << "Run: " << runName <<"  Efficiency: "<<efficiency <<"+/- "<<deviationP<<"/"<<deviationM<<"  All: " <<allHist->GetEntries()<<"  Cut: " <<cuttedHist->GetEntries()<< endl;

		printf("Run: %s, Efficiency: %1.4f [%1.4f,%1.4f] ", runName.c_str(), efficiency, deviationP, deviationM);

		TAxis *yaxisP = allHist->GetYaxis();
		TAxis *xaxisP = allHist->GetXaxis();
		yaxisP->SetLabelSize(0.04);
		yaxisP->SetTitle("counts");
		yaxisP->SetTitleSize(0.04);
		xaxisP->SetLabelSize(0.04);
		xaxisP->SetTitle("integral");
		xaxisP->SetTitleSize(0.04);

		allHist->Draw("hist");

		TLatex latex;
		latex.SetNDC(true);
		latex.SetTextSize(0.04);
		latex.SetTextAlign(13); //align at top
		//latex.DrawLatex(0.1, 0.93, Form("Efficiency: %s", runName.c_str()));

		//hs->Draw("nostack");
		effCanvas->Update();
		float minY = gPad->GetUymin();
		float maxY = gPad->GetUymax();
		TLine *limit = new TLine(threshold, minY, threshold, maxY);
		limit->SetLineColor(kRed);
		limit->SetLineWidth(3);

		int binOfThreshold = allHist->FindBin(threshold);

		TLegend *h_leg = new TLegend(0.47, 0.62, 0.90, 0.90);
		h_leg->SetTextSize(0.03);
		h_leg->AddEntry(allHist, Form("Entries: %1.1f, Cutted: %1.1f ", entry,entryC), "f");
		//h_leg->AddEntry(limit, Form("Threshold from: %s", thresholdName.c_str()), "l");
		h_leg->AddEntry(limit, Form("#Lambda_{thr}: %1.2f (+%1.2f -%1.2f)", threshold, sumCutErrorP,sumCutErrorM), "l");
		h_leg->AddEntry(limit, Form("P_{thr}: %1.2f (+%1.2f -%1.2f%s)", sumPercentage, sumPercentageErrorP,sumPercentageErrorM, "%"), "l");

		h_leg->AddEntry((TObject *)0, Form("Err Up :#Delta #varepsilon_{sys}:%1.4f, #Delta #varepsilon_{stat}: %1.4f", deviationP, upperStatErr), "");
		h_leg->AddEntry((TObject *)0, Form("Err Low :#Delta #varepsilon_{sys}:%1.4f, #Delta #varepsilon_{stat}: %1.4f", deviationM, lowerStatErr), "");
		h_leg->AddEntry((TObject *)0, Form("#bf{Efficiency: %1.4f%s (+%1.4lf,-%1.4lf)}",efficiency,"%", combinedUpperError, combinedLowerError), "");

		h_leg->SetTextFont(42);

		//if(!zoomMode)h_leg->Draw();
		//limit->Draw();

		//TCanvas *statCanvas = new TCanvas("statCanvas", "Stat Histogram", 1000, 1000);
		//statCanvas->SetGrid();
		/*TEfficiency *eff1;
		eff1 = new TEfficiency(*cuttedHist, *allHist);
		eff1->SetStatisticOption(TEfficiency::kFCP);
		eff1->Draw("AP");*/
		//float statEff=eff1->GetEfficiency(binOfThreshold);

		FILE *file_listP;
		string list_filenameP = outputFolder + "/Efficiencies.txt";
		file_listP = fopen(list_filenameP.c_str(), "a");

		fprintf(file_listP, "%s=%1.4f/%1.4lf/%1.4lf) \n", string(runName).c_str(), efficiency, combinedUpperError, combinedLowerError);
		fclose(file_listP);

	
		effCanvas->Print((outDir + runName + "_Efficiency.pdf").c_str());
	
	



	return 0;
}