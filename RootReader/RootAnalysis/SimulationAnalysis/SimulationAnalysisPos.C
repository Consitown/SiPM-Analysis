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

string outputFolder = "./SimulationPlots";
string effFile = "./Data.txt";
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
  //run_pos0_e14={3.495308,1.485189,42.490936,2477664}
		string effName = lineWithDataMeasurement.substr(0, lineWithDataMeasurement.find("=")); 
	
		vector<string> rawData = split(effName, "_");
		string runPosRaw = rawData[1];

		int runPos = stoi(runPosRaw.substr(3, runPosRaw.size() - 1));
		int runEnergy = stoi(rawData[2].substr(1, rawData[2].size()));
		string dataStr = lineWithDataMeasurement.substr(lineWithDataMeasurement.find("=") + 2, lineWithDataMeasurement.size()-1);
		vector<string> data = split(dataStr, ",");
		
		cout<<"DD: "<<data[0]<<"  "<<data[1]<<"  "<<data[2]<<"  "<<data[3] <<endl;

		float pReachedWOM=stof(data[0]); //Overall reached the WOM
		float pDetectedWOM=stof(data[1]); //Overall detected at the WOM
		float pSiPMs=stof(data[2]);
		int nPhotons=stoi(data[3]);

		float eSc=stoi(data[4]);
		float eW=stoi(data[5]);

		vector<float> valueTypes;

		valueTypes.push_back(pReachedWOM);
		valueTypes.push_back(pDetectedWOM);
		valueTypes.push_back(pSiPMs);
		valueTypes.push_back(nPhotons);
		valueTypes.push_back(eSc);
		valueTypes.push_back(eW);

		map<int, vector<float>> c = dataMap[runEnergy];
		c[runPos] = valueTypes;
		dataMap[runEnergy] = c;
	}
	//TLegend *h_leg = new TLegend(0.22, 0.12, 0.45, 0.30);
	TLegend *h_leg = new TLegend(0.6, 0.65, 0.85, 0.85);

	h_leg->SetTextSize(0.04);

	int colorCounter = 1;
	for (std::pair<int, map<int, vector<float>>> element : dataMap)
	{
		int energy = element.first;
		map<int, vector<float>> data = element.second;
		printf("Doing Energy: %d \n", energy);
		vector<float> x;
		vector<float> yDetected;
		vector<float> yReached;

		vector<float> yNumberOfPhotons;

		vector<float> yEsc;
		vector<float> yEW;


		for (std::pair<int, vector<float>> posRun : data)
		{
			int runPos = posRun.first;
			vector<float> dat = posRun.second;
			float d=getDistanceFromWOMD(runPos);
			x.push_back(d);
			yDetected.push_back(dat[1]);
			yReached.push_back(dat[0]);

			yNumberOfPhotons.push_back(dat[3]);

			yEsc.push_back(dat[4]);
			yEW.push_back(dat[5]);


			printf("Doing Position: %d , Detected: %1.4f, NoF: %1.4f, EDep Sc: %1.4f, Distance: %1.1f\n", runPos, dat[0], dat[3],dat[4],d);
		}




		//vector<float> y=yNumberOfPhotons;
		//vector<float> y=yReached;
		vector<float> y=yDetected;




		const Int_t n = x.size();
		TGraph *gr = new TGraph(n, &(x[0]), &(y[0]));
		gr->SetTitle("Energy");
		gr->SetLineColor(colorCounter);
		gr->SetLineWidth(2.5);
		gr->SetMarkerColor(colorCounter);
		gr->SetMarkerSize(1.0);
		gr->SetMarkerStyle(19+colorCounter);
		gr->GetXaxis()->SetLabelSize(0.085);
		gr->GetYaxis()->SetLabelSize(0.085);
		h_leg->AddEntry(gr, Form("Energy: %1.2f GeV", energy / 10.0), "p");
		mg->Add(gr);

		/*TGraph *gr2 = new TGraph(n, &(x[0]), &(yReached[0]));
		gr2->SetTitle("Energy");
		gr2->SetLineColor(colorCounter);
		gr2->SetLineWidth(2.5);
		gr2->SetMarkerColor(colorCounter);
		gr2->SetMarkerSize(1.0);
		gr2->SetFillColorAlpha(3,1);
		gr2->SetMarkerStyle(19+colorCounter);
		gr2->GetXaxis()->SetLabelSize(0.085);
		gr2->GetYaxis()->SetLabelSize(0.085);		
		mg->Add(gr2);
		*/











		gStyle->SetTitleSize(0.2);
		//mg->SetTitle("Efficiency Measurement");
	
		colorCounter++;
	}
	mg->Draw("AP");
	//mg->SetTitle("number of photons created inside the scintillator");
	//mg->SetTitle("number of photons detected at WOM D");
	//mg->SetTitle("number of photons that reached WOM D");
	mg->SetTitle("energy deposition in walls");


	TAxis *yaxisP = mg->GetYaxis();
	TAxis *xaxisP = mg->GetXaxis();
	yaxisP->SetLabelSize(0.04);
	yaxisP->SetTitle("Normalized EDep[%]");
	yaxisP->SetTitleSize(0.04);
	xaxisP->SetLabelSize(0.04);
	xaxisP->SetTitle("distance [mm]");
	xaxisP->SetTitleSize(0.04);
	h_leg->SetTextFont(42);
	h_leg->Draw();

	effCanvas->Print((outDir + "PosSimulation.pdf").c_str());
	

	return 0;
}