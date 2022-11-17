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
	while (!ss.eof())
	{

		/* extracting word by word from stream */
		ss >> temp;

		/* Checking the given word is integer or not */
		if (stringstream(temp) >> found)
			cout << found << " ";

		/* To save from space at the end of string */
		temp = "";
	}
}

int main(int argc, char *argv[])
{

	TString c_name1, c_name2, h_name1, h_name2, draw_cmnd1, draw_cmnd2, cut_cmnd2;

	string fileLocation = string(argv[1]);
	vector<string> splitted = split(string(argv[1]), "//");
	string run_name = splitted[1];
	vector<string> runInformations = split(run_name.c_str(), "_");
	cout << "Doing:  " << argv[1] << " runName: " << run_name << endl;

	float energy = 1.4;
	if (runInformations[2].find("5") != std::string::npos)
	{
		energy = 5.2;
	}
	else if (runInformations[2].find("2") != std::string::npos)
	{
		energy = 2.6;
	}
	int position = -1;
	string ps = split(runInformations[3], ".")[0];
	string posraw = ps.substr(8, ps.length());
	//cout<< "Doing:  "<<posraw<<endl;
	position = stoi(posraw);
	cout << "pos: " << position << "  energy: " << energy << endl;
	//wls_e_5.2GeV_position0

	TApplication *ROOTapp;

	TFile *file = new TFile(fileLocation.c_str());
	TTree *treeAll;
	file->GetObject("EventStat", treeAll);

	Int_t sp, cp;
	treeAll->SetBranchAddress("scintillation_photons", &sp);
	treeAll->SetBranchAddress("Cerenkov_photons", &cp);

	TH1F *spH = new TH1F("sp", "sp", 100, 0, -30000000);
	TH1F *cpH = new TH1F("cp", "cp", 100, 0, -30000000);

	Int_t nentries = (Int_t)treeAll->GetEntries();

	for (Int_t i = 0; i < nentries; i++)
	{
		treeAll->GetEntry(i);
		spH->Fill(sp);
		cpH->Fill(cp);
	}

	float meanSP = spH->GetMean();
	float meanCP = cpH->GetMean();
	float totalPhotons = meanSP + meanCP;

	spH->Draw();
	cout << "Means Total:" << totalPhotons << " Scintillator: " << meanSP << "  Cherenkov: " << meanCP << endl;

	TTree *treeDetected;
	file->GetObject("Detected", treeDetected);
	Int_t detected, WOM;
	treeDetected->SetBranchAddress("detection", &detected);
	treeDetected->SetBranchAddress("WOWnumber", &WOM);

	int atWOMD = 0;
	int detectedAtWOMD = 0;

	nentries = (Int_t)treeDetected->GetEntries();

	for (Int_t i = 0; i < nentries; i++)
	{
		treeDetected->GetEntry(i);
		if (WOM == 1)
		{
			//WOM D
			atWOMD++;
			if (detected == 1)
			{
				//Detected at WOM D
				detectedAtWOMD++;
			}
		}
	}

	float pReachedWOM = (atWOMD / totalPhotons) * 100.0;				 //Overall reached the WOM
	float pDetectedWOM = ((float)detectedAtWOMD / totalPhotons) * 100.0; //Overall detected at the WOM
	//float pDetectedWOM = detectedAtWOMD; //Overall detected at the WOM

	float pSiPMs = ((float)detectedAtWOMD / atWOMD) * 100.0; //Overall detected at the WOM

	cout << "Photons at reached WOM D: " << atWOMD << "  Detected: " << detectedAtWOMD << endl;
	cout << "Percentages:  ReachedWOM/Detected/SiPM:  " << pReachedWOM << "/" << pDetectedWOM << "/" << pSiPMs << endl;

	string formattedName = "run_pos" + to_string(position) + "_e" + to_string((int)(energy * 10));

	TTree *treeEnergy;
	file->GetObject("Energy deposition", treeEnergy);
	Int_t volume;
	Double_t Edep;
	treeEnergy->SetBranchAddress("volume", &volume); //1 scintillator, 2 walls
	treeEnergy->SetBranchAddress("Edep", &Edep);





	TCanvas *effCanvas = new TCanvas("effCanvas", "Sum Histogram", 1000, 600);
	effCanvas->SetGrid();
	// canvas style
	gPad->SetRightMargin(0.00);
	gPad->SetLeftMargin(.12);
	int n_bins = 150;
	// int n_bins = (int)(Xmax - Xmin)*1;
	TH1F* hist = new TH1F("sd","Histtt", 400,0,-1);

	hist->GetXaxis()->SetTitle("pulse-height(charge) [mv #times ns]");
	// h_vec[i]->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}]");
	hist->GetXaxis()->SetTitleOffset(1.3);
	treeEnergy->Draw("Edep","volume==1","");
	effCanvas->SaveAs("out.pdf");






	double eDepSc = 0;
	double eDepW = 0;
	double eDepElse = 0;

	double eDepAll = 0; //Only percentages no MeV

	nentries = (Int_t)treeEnergy->GetEntries();

	for (Int_t i = 0; i < nentries; i++)
	{
		treeEnergy->GetEntry(i);
		eDepAll += Edep;
		if (volume == 1)
		{
			eDepSc += Edep;
		}
		if (volume == 2)
		{
			eDepW += Edep;
		}
		if (volume == 0)
		{
			eDepElse += Edep;
		}
	}

	float pESc = ((float)eDepSc / eDepAll) * 100.0; //Edep in scintillator percentage
	float pEW = ((float)eDepW / eDepAll) * 100.0;
	float pEElse = ((float)eDepElse / eDepAll) * 100.0;

	cout << "Energy Deposition ABS:  Scintillator/Walls/All:  " << eDepSc * (1.0 / 500) << "/" << eDepW * (1.0 / 500) << "/" << eDepAll * (1.0 / 500) << endl;
	cout << "Energy Deposition:  Scintillator/Walls/Else:  " << pESc << "/" << pEW << "/" << pEElse << endl;

	FILE *pFile = fopen("Data.txt", "a");
	fprintf(pFile, "%s={%f,%f,%f,%d,%f,%f,%f}\n", string(formattedName).c_str(), pReachedWOM, pDetectedWOM, pSiPMs, (int)totalPhotons, eDepSc, eDepW , eDepElse);
	fclose(pFile);

	cout << "" << endl;

	return 0;
}