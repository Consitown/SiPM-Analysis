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
#include <TGraphAsymmErrors.h>
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

string outputFolder = "./EfficiencyPlots";
string effFile = outputFolder + "/Efficiencies.txt";
pair<int, int> pos0 = make_pair(0, 0);
pair<int, int> pos1 = make_pair(160, 160);
pair<int, int> pos2 = make_pair(320, 320);
pair<int, int> pos3 = make_pair(-320, 320);
pair<int, int> pos4 = make_pair(-160, 160);
pair<int, int> pos5 = make_pair(-160, -160);
pair<int, int> pos6 = make_pair(-320, -320);
pair<int, int> pos7 = make_pair(320, -320);
pair<int, int> pos8 = make_pair(160, -160);
pair<int, int> pos9 = make_pair(160, -510);
pair<int, int> pos10 = make_pair(204, -404);
pair<int, int> pos11 = make_pair(310, -360);
pair<int, int> pos12 = make_pair(-320, 0); //Coordinates mapped on the box, not the unit coordinates
pair<int, int> pos13 = make_pair(320, 0);
pair<int, int> pos14 = make_pair(0, 407);
pair<int, int> pos15 = make_pair(0, -407);
pair<int, int> pos16 = make_pair(0, -160);
pair<int, int> pos17 = make_pair(0, 160);
pair<int, int> pos18 = make_pair(-160, 0);
pair<int, int> pos19 = make_pair(160, 0);
pair<int, int> pos20 = make_pair(320, -510);
pair<int, int> pos21 = make_pair(280, -510);
pair<int, int> pos22 = make_pair(0, -510);
pair<int, int> pos23 = make_pair(320, -510);

pair<int, int> posWOMD = make_pair(310, -510);
float calculateDistance(pair<int, int> p1, pair<int, int> p2)
{
	return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}
float getDistanceFromWOMD(int pos)
{
	//pos of WOM D: x/y -> 310,-510

	switch (pos)
	{
	case 0:
		return calculateDistance(posWOMD, pos0);
	case 1:
		return calculateDistance(posWOMD, pos1);
	case 2:
		return calculateDistance(posWOMD, pos2);
	case 3:
		return calculateDistance(posWOMD, pos3);
	case 4:
		return calculateDistance(posWOMD, pos4);
	case 5:
		return calculateDistance(posWOMD, pos5);
	case 6:
		return calculateDistance(posWOMD, pos6);
	case 7:
		return calculateDistance(posWOMD, pos7);
	case 8:
		return calculateDistance(posWOMD, pos8);
	case 9:
		return calculateDistance(posWOMD, pos9);
	case 10:
		return calculateDistance(posWOMD, pos10);
	case 11:
		return calculateDistance(posWOMD, pos11);
	case 12:
		return calculateDistance(posWOMD, pos12);
	case 13:
		return calculateDistance(posWOMD, pos13);
	case 14:
		return calculateDistance(posWOMD, pos14);
	case 15:
		return calculateDistance(posWOMD, pos15);
	case 16:
		return calculateDistance(posWOMD, pos16);
	case 17:
		return calculateDistance(posWOMD, pos17);
	case 18:
		return calculateDistance(posWOMD, pos18);
	case 19:
		return calculateDistance(posWOMD, pos19);
	case 20:
		return calculateDistance(posWOMD, pos20);
	case 21:
		return calculateDistance(posWOMD, pos21);
	case 22:
		return calculateDistance(posWOMD, pos22);
	case 23:
		return calculateDistance(posWOMD, pos23);

	default:
		return -1000;
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

	TApplication *myapp = new TApplication("myapp", 0, 0);

	//Output Plot Folder
	string outDir = outputFolder + "/";
	if (!fs::is_directory(outputFolder) || !fs::exists(outputFolder))
	{
		fs::create_directory(outputFolder);
	}
	if (!fs::is_directory(outDir) || !fs::exists(outDir))
	{
		fs::create_directory(outDir);
	}

	std::ifstream inFile(effFile);
	int numberOfLines = std::count(std::istreambuf_iterator<char>(inFile),
								   std::istreambuf_iterator<char>(), '\n');

	ifstream data_store(effFile);

	TCanvas *effCanvas = new TCanvas("effCanvas", "Sum Histogram", 1000, 600);
	effCanvas->SetGrid();

	map<int, map<int, vector<float>>> dataMap; //Key= Energy , Val: Runs<ValueTypes<Values>>>

	TMultiGraph *mg = new TMultiGraph();

	for (string lineWithDataMeasurement; getline(data_store, lineWithDataMeasurement);)
	{

		string effName = lineWithDataMeasurement.substr(0, lineWithDataMeasurement.find("=")); //12_pos3_angle0_e52_ch32
		if (effName.find("dc") != std::string::npos)
		{
			continue;
		}
		vector<string> rawData = split(effName, "_");
		string runPosRaw = rawData[1];

		int runPos = stoi(runPosRaw.substr(3, runPosRaw.size() - 1));
		int runEnergy = stoi(rawData[3].substr(1, rawData[3].size()));
		string dataStr = lineWithDataMeasurement.substr(lineWithDataMeasurement.find("=") + 1, lineWithDataMeasurement.size());
		vector<string> data = split(dataStr, "/");
		float eff = stof(data[0]);
		float errP = stof(data[1]);
		float errM = stof(data[2]);

		vector<float> valueTypes;

		valueTypes.push_back(eff);
		valueTypes.push_back(errP);
		valueTypes.push_back(errM);

		map<int, vector<float>> c = dataMap[runEnergy];
		c[runPos] = valueTypes;
		dataMap[runEnergy] = c;
	}
	TLegend *h_leg = new TLegend(0.12, 0.12, 0.35, 0.30);
	h_leg->SetTextSize(0.04);

	int colorCounter = 1;
	for (std::pair<int, map<int, vector<float>>> element : dataMap)
	{
		int energy = element.first;
		map<int, vector<float>> data = element.second;
		printf("Doing Energy: %d \n", energy);
		vector<float> x;
		vector<float> y;
		vector<float> yErrU;
		vector<float> yErrL;

		for (std::pair<int, vector<float>> posRun : data)
		{
			int runPos = posRun.first;
			vector<float> dat = posRun.second;
			float d = getDistanceFromWOMD(runPos);
			x.push_back(d);
			y.push_back(dat[0]);
			yErrU.push_back(dat[1]);
			yErrL.push_back(dat[2]);
			printf("Doing Position: %d , Eff: %1.4f (%1.4f,%1.4f), Distance: %1.1f\n", runPos, dat[0], dat[1], dat[2], d);
		}

		const Int_t n = x.size();
		TGraphAsymmErrors *gr = new TGraphAsymmErrors(n, &(x[0]), &(y[0]), 0, 0, &(yErrL[0]), &(yErrU[0]));

		if (colorCounter == 1)
		{
			gr->SetLineColor(2);
			gr->SetMarkerStyle(21);
			gr->SetMarkerColor(2);
		}
		else if (colorCounter == 2)
		{
			gr->SetLineColor(4);
			gr->SetMarkerStyle(22);
			gr->SetMarkerColor(4);
		}
		else
		{
			gr->SetLineColor(1);
			gr->SetMarkerStyle(20);
			gr->SetMarkerColor(1);
		}
		gr->SetTitle("Energy");
		gr->SetLineWidth(2.5);
		gr->SetMarkerSize(1.2);
		gr->GetXaxis()->SetLabelSize(0.085);
		gr->GetYaxis()->SetLabelSize(0.085);

		h_leg->AddEntry(gr, Form("Energy: %1.2f GeV", energy / 10.0), "p");

		TF1 *f1 = new TF1("fit", "-[0]*exp(x*[2])+[1]", 160, 2000);
		f1->SetParameter(0, 0.01);
		f1->SetParameter(1, 100);
		f1->SetParameter(2, 0.009);
		//gr->Fit("fit", "R");

		mg->Add(gr);
		gStyle->SetTitleSize(0.2);
		//mg->SetTitle("Efficiency Measurement");

		colorCounter++;
	}
	mg->Draw("AP ");

	TAxis *yaxisP = mg->GetYaxis();
	TAxis *xaxisP = mg->GetXaxis();
	yaxisP->SetLabelSize(0.04);
	yaxisP->SetTitle("efficiency [%]");
	yaxisP->SetTitleSize(0.04);
	xaxisP->SetLabelSize(0.04);
	xaxisP->SetTitle("d [mm]");
	xaxisP->SetTitleSize(0.04);
	h_leg->SetTextFont(42);
	h_leg->Draw();

	effCanvas->Print((outDir + "PosEfficiency.pdf").c_str());
	cout << "This needs to be closed manual, since the ROOT app is running! -> Strg + C" << endl;
	myapp->Run();

	return 0;
}