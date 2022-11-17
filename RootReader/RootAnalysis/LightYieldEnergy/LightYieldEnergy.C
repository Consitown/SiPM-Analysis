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

string outputFolder = "./LightYieldEnergyPlots";

string vectorToString(vector<string> vec)
{
	std::ostringstream vts;

	if (!vec.empty())
	{
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec.begin(), vec.end() - 1,
				  std::ostream_iterator<string>(vts, ", "));

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

float lightyield1=0;
float lightyield2=0;
float lightyield3=0;
int main(int argc, char *argv[])
{
	gStyle->SetOptStat(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
	gStyle->SetOptTitle(0);
	gStyle->SetLineStyleString(11, "5 5");
	gErrorIgnoreLevel = kWarning;

	vector<string> runNames14 = split((string)argv[1], ",");
	vector<string> runNames26 = split((string)argv[2], ",");
	vector<string> runNames52 = split((string)argv[3], ",");
	vector<string> positions = split((string)argv[4], ",");

	TCanvas *masterCanvas = new TCanvas("all", "Sum Histogram", 2000, 3000);
	masterCanvas->Divide(3, 4);
	masterCanvas->SetTopMargin(0.);
	masterCanvas->SetBottomMargin(0.);
	masterCanvas->SetLeftMargin(0.);
	masterCanvas->SetRightMargin(0.);

	string outDirAll = outputFolder + "/";

	for (size_t i = 0; i < positions.size(); i++)
	{
		masterCanvas->cd(i + 1);

		int position = stoi(positions[i]);
		int angle = 0;

		vector<string> runList;
		//Find energies for positions
		cout << "Position: " << positions[i] << endl;

		string printName = "Position_" + positions[i];

		string outDir = Form((outputFolder + "/%s/").c_str(), printName.c_str());
	/*	if (!fs::is_directory(outputFolder) || !fs::exists(outputFolder))
		{
			fs::create_directory(outputFolder);
		}
		if (!fs::is_directory(outDir) || !fs::exists(outDir))
		{
			fs::create_directory(outDir);
		}*/

		for (std::vector<string>::size_type j = 0; j != runNames14.size(); j++)
		{

			string runtemp = runNames14[j];
			string ser = string("pos" + positions[i]+"_");

			if (runtemp.find(ser) != std::string::npos)
			{
				runList.push_back(runtemp);
				break;
			}
		}
		for (std::vector<string>::size_type j = 0; j != runNames26.size(); j++)
		{

			string runtemp = runNames26[j];
			string ser = string("pos" + positions[i]+"_");

			if (runtemp.find(ser) != std::string::npos)
			{
				runList.push_back(runtemp);
				break;
			}
		}

		for (std::vector<string>::size_type j = 0; j != runNames52.size(); j++)
		{

			string runtemp = runNames52[j];
			string ser = string("pos" + positions[i]+"_");

			if (runtemp.find(ser) != std::string::npos)
			{
				runList.push_back(runtemp);
				break;
			}
		}

		std::cout << "RunList done: " << vectorToString(runList) << endl;
		;
		cout << endl;

		TCanvas *effCanvas = new TCanvas("effCanvas", "Sum Histogram", 1000, 500);
		effCanvas->SetGrid();
		effCanvas->cd();
		THStack *hs = new THStack("hs", "");

		TLegend *h_leg = new TLegend(0.44, 0.72, 1, 1);
		h_leg->SetTextSize(0.04);

		int colorCounter = 1;
		
		for (std::vector<string>::size_type j = 0; j != runList.size(); j++)
		{

			string runName = runList[j];
			string inDir = "../rootfiles/" + runName + ".root";

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

			int channelNumber = 16;
			int entryNumber = 16;
			float energy = 0;
			int position = 0;
			angle = 0;

			file->GetObject("T", tree);
			tree->SetBranchAddress("nCh", &channelNumber);
			tree->SetBranchAddress("runEnergy", &energy);
			tree->SetBranchAddress("runPosition", &position);
			tree->SetBranchAddress("runAngle", &angle);

			tree->GetEntry(1); //sollte alles konstant sein, also kann man irgendeinen Entry nehmen.

			energy = energy * 0.1f;
			cout << "Doing: " << runName << "  Energy: " << energy << "  Pos: " << position << endl;

			int maxX = -900;
			int min = -100;
			int bins = 450;

			TH1D *allHist = new TH1D("allHist", "", bins, min, maxX);

			allHist->SetFillColorAlpha(colorCounter, 0.4);
			allHist->SetLineColorAlpha(1, 0);
			allHist->SetLineColor(colorCounter);
			//Fill
			tree->Draw("chargeChannelSumWOM[3]>>allHist", "");

			float mean = allHist->GetMean();
			//allHist->Draw("hist");
			h_leg->AddEntry(allHist, Form("%1.1f GeV  Mean: %1.0f N_{pe}", energy, mean), "f");

			hs->Add(allHist);


			switch (colorCounter)
			{
			case 1:
					lightyield1=mean;
				break;
			case 2:
					lightyield2=mean;
				break;
				case 3:
					lightyield3=mean;
				break;
		
			}
			colorCounter++;



		}

		FILE *file_listP;
		string list_filenameP = outputFolder + "/Lightyields.txt";
		file_listP = fopen(list_filenameP.c_str(), "a");

		fprintf(file_listP, "%d_%d=%1.2f/%1.2lf/%1.2lf \n", position,angle, lightyield1,lightyield2,lightyield3);
		fclose(file_listP);

		hs->Draw("nostack");
		TAxis *yaxisP = hs->GetYaxis();
		TAxis *xaxisP = hs->GetXaxis();
		yaxisP->SetLabelSize(0.04);
		yaxisP->SetTitle("counts");
		yaxisP->SetTitleSize(0.04);
		xaxisP->SetLabelSize(0.04);
		xaxisP->SetTitle("N_{pe}");
		xaxisP->SetTitleSize(0.04);

		//h_leg->AddEntry(limit, Form("Threshold from: %s", thresholdName.c_str()), "l");
		//h_leg->AddEntry(limit, Form("#Lambda_{thr}: %1.2f (+%1.2f -%1.2f)", threshold, sumCutErrorP,sumCutErrorM), "l");
		//h_leg->AddEntry(limit, Form("P_{thr}: %1.2f (+%1.2f -%1.2f%s)", sumPercentage, sumPercentageErrorP,sumPercentageErrorM, "%"), "l")
		h_leg->AddEntry((TObject *)0, Form("Position: %d", position), "");

		h_leg->Draw();
		masterCanvas->cd(i + 1);
		effCanvas->DrawClonePad();
	}
	masterCanvas->Print((outDirAll + "LYE.pdf").c_str());
	return 0;
}