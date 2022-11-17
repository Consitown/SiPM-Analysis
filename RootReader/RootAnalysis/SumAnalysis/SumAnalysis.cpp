//Start ohne Compile
//root -l ./MaskAnalysis

/**
 * The point of this is to analyze the light yield changes when using a different amount of sipms per channel.
 * Light Yield is defined as the mean value of the charge/amplitude histograms 
 **/

//root
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TF1.h>
#include <TString.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TText.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TVirtualPad.h>
//C, C++
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <vector>
#include <map>
using namespace std;
struct stat st = {0};

int nBins = 250;
int xmin = -1000;
int xmax = -5000000; //AUTOBINNING XMAX<XMIN

string parent_dir = "./SumAnalysis/";
char path[1000] = "../rootfiles/";
int charge = 0; //AMP-CHARGE

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
int main(int argc, char *argv[])
{
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(1110);
	gStyle->SetOptFit(1);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1);
	gStyle->SetStatW(0.3); 
	gStyle->SetStatH(0.15);

	gErrorIgnoreLevel = kError;

	int array_size = 0;
	DIR *d;
	struct dirent *dirs;
	d = opendir(path);
	char full_path[1000];
	char runName[1000];
	string workingDir;
	string calib_path;
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
	{
		printf("Current working dir: %s\n", cwd);
		workingDir = cwd;
	}
	else
	{
		perror("getcwd() error");
	}

	int file_count = 0;
	DIR *dirp;
	struct dirent *entry;

	dirp = opendir(path); /* There should be error handling after this */
	while ((entry = readdir(dirp)) != NULL)
	{
		if (entry->d_type == DT_REG)
		{ /* If the entry is a regular file */
			file_count++;
		}
	}
	closedir(dirp);

	if (d)
	{
		int rootfile_iteration = 1;

		while ((dirs = readdir(d)) != NULL)
		{
			if (dirs->d_type != DT_REG)
				continue;

			full_path[0] = '\0';
			runName[0] = '\0';

			strcat(full_path, path);
			strcat(full_path, "/");
			strcat(full_path, dirs->d_name);
			strcat(runName, dirs->d_name);

			TString g = full_path;
			TFile *file0 = new TFile(full_path);
			strip_ext(runName);

			TTree *tree;
			file0->GetObject("T", tree);

			TCanvas *parentCanvas = new TCanvas("PC", "Parent", 2000, 1000);
			parentCanvas->Divide(2, 1, 0.03, 0.03);

			TText *title = new TText(0.01, 0.96, runName);
			title->SetNDC();
			title->SetTextSize(0.030);
			title->Draw();
			for (int j = 0; j < 2; j++)
			{
				charge = j;

				int channel = 0;
				int i = 0;

				TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 1000);
				masterCanvas->Divide(2, 2, 0.03, 0.03);
				masterCanvas->SetGrid();

				string type = "charge";
				for (i = 0; i < 4; i++)
				{
					TVirtualPad *pad = masterCanvas->cd(i + 1);
					pad->SetGrid();
					string label = "WOM A";
					switch (i)
					{
					case 1:
						label = "WOM B";
						break;
					case 2:
						label = "WOM C";
						break;
					case 3:
						label = "WOM D";
						break;
					default:
						break;
					}

					TH1F *h = new TH1F("h", label.c_str(), nBins, xmin, xmax);
					h->GetYaxis()->SetLabelOffset(0.015);
					h->GetYaxis()->SetTitleOffset(1.65);

					h->GetXaxis()->SetLabelOffset(0.015);

					h->GetYaxis()->SetTitle("entries");
					if (j == 0)
						h->GetXaxis()->SetTitle("mv");
					else
					{
						h->GetXaxis()->SetTitle("mv*ns");
					}

					h->GetYaxis()->SetLabelSize(0.03);
					h->GetXaxis()->SetLabelSize(0.03);

					TString cut("");

					if (charge == 1)
					{
						tree->Draw(Form("chargeChannelSumWOM[%d]>>h", i), cut); //qS2 integration window width 55ns
						type = "Charge";
					}
					else
					{
						tree->Draw(Form("amplitudeChannelSumWOM[%d]>>h", i), cut); //qS2 integration window width 55ns
						type = "Amplitude";
					}





					double mean = h->GetMean();
					TLine *l = new TLine(mean, 0, mean, 5);
					gStyle->SetLineColor(2);
					//	l->Draw();
					TLegend *leg = new TLegend(0.7, 0.8, 0.90, 0.9);

					char arr[sizeof(mean)];
					memcpy(arr, &mean, sizeof(mean));

					char *valueStr = (char *)malloc(13 * sizeof(char));

					sprintf(valueStr, "Mean: %d ", (int)round(mean));

					gStyle->SetLegendTextSize(0.030);
					leg->AddEntry("Mean", valueStr, "");
					//leg->Draw();
				}




				masterCanvas->cd();

				char *valueStr = (char *)malloc(50 * sizeof(char));
				sprintf(valueStr, "%s", type.c_str());
				TText *t = new TText(0.00, 0.005, valueStr);
				t->SetTextSize(0.025);
				t->Draw();

				parentCanvas->cd(j + 1);
				masterCanvas->DrawClonePad();

				string target_dir = parent_dir + runName + "/";

				if (stat("./SumAnalysis", &st) == -1)
				{
					mkdir("./SumAnalysis", 0700);
				}
				if (stat(parent_dir.c_str(), &st) == -1)
				{
					mkdir(parent_dir.c_str(), 0700);
				}
				if (stat(target_dir.c_str(), &st) == -1)
				{
					mkdir(target_dir.c_str(), 0700);
				}
			}

			cout << "File: " << rootfile_iteration << "/" << file_count << endl;
			if (rootfile_iteration == 1 && rootfile_iteration != file_count) //START
			{
				parentCanvas->Print((parent_dir + "overviewSumAnalysis.pdf(").c_str());
			}
			else if (rootfile_iteration == file_count && file_count > 1) //END
			{
				parentCanvas->Print((parent_dir + "overviewSumAnalysis.pdf)").c_str());
			}
			else
			{
				parentCanvas->Print((parent_dir + "overviewSumAnalysis.pdf").c_str());
			}
			rootfile_iteration++;
		}

		closedir(d);
	}

	return 0;
}