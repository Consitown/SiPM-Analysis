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
int xmin = -80;
int xmax = -1000; //AUTOBINNING XMAX<XMIN
int channels[8] = {0, 1, 2, 3, 4, 5, 6, 7};
string parent_dir = "./MaskAnalysis/";
int charge = 0; //AMP-CHARGE
//vector<vector<double>> values;   //1-Channels -> 2-Sipm count (filename exracted: sipmX) -> 3-Data
map<int, map<int, map<int, double>>> values; //0-> Charge 1->Channel 2->Filev

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

	int array_size = 0;
	DIR *d;
	struct dirent *dirs;
	char path[1000] = "../rootfiles/";
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

	dirp = opendir("../rootfiles/"); /* There should be error handling after this */
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
			for (int j = 0; j < 2; j++)
			{
				charge = j;
				//Condition to check regular file.
				if (dirs->d_type == DT_REG)
				{

					full_path[0] = '\0';
					runName[0] = '\0';

					strcat(full_path, path);
					strcat(full_path, "/");
					strcat(full_path, dirs->d_name);
					strcat(runName, dirs->d_name);

					printf("%s\n", full_path);
					TString g = full_path;
					TFile *file0 = new TFile(full_path);
					strip_ext(runName);

					TTree *tree;
					file0->GetObject("T", tree);
					gStyle->SetOptStat(0);
					gStyle->SetOptFit(1);
					gStyle->SetStatW(0.2);
					gStyle->SetStatH(0.1); // stats box size
					gStyle->SetGridStyle(3);
					gStyle->SetGridWidth(1);
					gStyle->SetGridColor(16);

					int channel = 0;
					int i = 0;

					TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 1000);
					masterCanvas->Divide(3, 3);
					string type = "charge";
					for (i = 0; i < sizeof(channels) / sizeof(int); i++)
					{
						masterCanvas->cd(i + 1);
						channel = channels[i];

						TH1F *h = new TH1F("h", "PE spectrum - integration window 35 ns", nBins, xmin, xmax);

						TString cut("");

						if (charge == 1)
						{
							tree->Draw(Form("Integral_inRange[%d]>>h", channel), cut); //qS2 integration window width 55ns
							type = "charge";
						}
						else
						{
							tree->Draw(Form("amp_inRange[%d]>>h", channel), cut); //qS2 integration window width 55ns
							type = "Amplitude";
						}
						double mean = h->GetMean();
						TLine *l = new TLine(mean, 0, mean, 1000);
						gStyle->SetLineColor(2);
						l->Draw();
						TLegend *leg = new TLegend(0.7, 0.8, 0.95, 0.9);

						char arr[sizeof(mean)];
						memcpy(arr, &mean, sizeof(mean));

						char *valueStr = (char *)malloc(13 * sizeof(char));
						sprintf(valueStr, "Mean: %lf ", mean);

						gStyle->SetLegendTextSize(0.025);
						leg->AddEntry("Mean", valueStr, "");
						leg->Draw();

						//SIPM Information
						//Naming Scheme: calib_vb58_tune90_3SiPMs_sn2_mask2
						//mask not there ->SIPM information

						char copyRunName[100];
						strcpy(copyRunName, runName); //OTHERWISE CONSUMED
						//string copyRunName(str2);
						int init_size = strlen(copyRunName);
						char delim[] = "_";

						char *ptr = strtok(copyRunName, delim);
						int maskNumber = -1;
						int sipmNumber = -1;
						int finalsipmNumber = 0;
						while (ptr != NULL)
						{
							ptr = strtok(NULL, delim);

							if (ptr != NULL)
							{
								string maskInfo(ptr);
								string sipmInfo(ptr);

								if (maskInfo.find("mask") != string::npos)
								{
									sscanf(ptr, "mask%d", &maskNumber);
								}
								if (sipmInfo.find("SiPMs") != string::npos)
								{
									sscanf(ptr, "%dSiPMs", &sipmNumber);
								}
							}
						}

						if (maskNumber == -1)
						{
							//NO MASK INFORMATION, SiPM Number=XSiPM
							finalsipmNumber = sipmNumber;
						}
						else
						{
							finalsipmNumber = maskNumber;
						}

						values[charge][channel][finalsipmNumber] = mean;
						//cout << "PRINT: " << i << " " << finalsipmNumber <<"     "  <<values[i][finalsipmNumber]<< endl;
					}

					masterCanvas->cd(9);

					char *valueStr = (char *)malloc(50 * sizeof(char));
					sprintf(valueStr, "%s | %s ", type.c_str(), runName);
					TText *t = new TText(.1, .5, valueStr);
					t->SetTextSize(0.035);
					t->Draw();

					//	string target_dir = parent_dir + runName + "/";
					string target_dir = parent_dir + runName + "/";
					string overview_filename = target_dir + "overview_" + runName + ".pdf";

					if (stat("./MaskAnalysis", &st) == -1)
					{
						mkdir("./MaskAnalysis", 0700);
					}
					if (stat(parent_dir.c_str(), &st) == -1)
					{
						mkdir(parent_dir.c_str(), 0700);
					}
					if (stat(target_dir.c_str(), &st) == -1)
					{
						mkdir(target_dir.c_str(), 0700);
					}

					if (rootfile_iteration == 1 && j == 0) //START
					{
						masterCanvas->Print((parent_dir + "overviewMaskAnalysis.pdf(").c_str());
					}
					else if (rootfile_iteration == file_count && j == 1) //END
					{
						masterCanvas->Print((parent_dir + "overviewMaskAnalysis.pdf)").c_str());
					}
					else
					{
						masterCanvas->Print((parent_dir + "overviewMaskAnalysis.pdf").c_str());
					}

					if (j == 1)
					{ //Directories also increment -> call it here inside the charge/amp loop but only once
						rootfile_iteration++;
					}
				}
			}
		}

		closedir(d);
	} //FILELOOP END

	//Summaries
	int charge = 0;
	string type = "charge";
	for (int j = 0; j < 2; j++)
	{
		charge = j;
		if (j == 0)
		{
			type = "amplitude";
		}
		else
			type = "charge";

		TCanvas *summaryCanvas = new TCanvas("SUMMARY", "Overview", 1000, 1000);
		summaryCanvas->Divide(3, 3);
		int i = 0;
		for (i = 0; i < sizeof(channels) / sizeof(int); i++)
		{
			summaryCanvas->cd(i + 1);
			int channel = channels[i];

			map<int, double> infoInChannel = values[charge][channel];

			vector<int> keys;
			vector<double> data;

			for (map<int, double>::iterator it = infoInChannel.begin(); it != infoInChannel.end(); ++it)
			{
				keys.push_back(it->first);
				data.push_back(it->second);
			}

			float scaleBy = data[0];
			for (std::vector<double>::iterator iter = data.begin(); iter != data.end(); ++iter)
			{
				int index = iter - data.begin();
				float toScale = data[index];
				data[index] = toScale / scaleBy;
			}

			vector<float> x(keys.begin(), keys.end());
			vector<float> y(data.begin(), data.end());
			const Int_t n = x.size();

			TGraphErrors *gr = new TGraphErrors(n, &(x[0]), &(y[0]), 0, 0);
			gr->SetTitle("");
			gr->SetLineColor(4);
			gr->SetLineWidth(1);
			gr->SetMarkerColor(1);
			//gr->SetMarkerColor(4);
			gr->SetMarkerSize(1.0);
			gr->SetMarkerStyle(20);
			gr->GetXaxis()->SetLabelSize(0.035);
			gr->GetYaxis()->SetLabelSize(0.035);
			gr->GetXaxis()->SetTitle("Number of SiPMs");
			gr->GetYaxis()->SetTitle("Light Yield");
			//	legend->AddEntry(gr, ("Position: " + to_string(it->first) + " WOM: " + womIterator->first).data(), "lep");
			gStyle->SetOptFit(1111);

			TF1 *func = new TF1("fit", "pol1");
			func->SetParName(1, "Slope");
			func->SetParName(0, "Offset");
			gr->Fit("fit");
			gr->Draw();
			gPad->Update();
			auto stat = dynamic_cast<TPaveStats *>(gr->FindObject("stats"));
			if (stat)
			{
				stat->SetX1NDC(0.15);
				stat->SetX2NDC(0.45);
				stat->Draw();
			}
			else
			{
				cout << "No stats box found!\n";
			}
		}

		map<int, double> infoInChannel = values[charge][0];
		vector<int> keys;
		vector<float> sumArr; // KEY: SiPM Number DATA: Sum of all Values

		for (map<int, double>::iterator it = infoInChannel.begin(); it != infoInChannel.end(); ++it)
		{
			keys.push_back(it->first);
			sumArr.push_back(0.0); //only need to make sumArr the right size to use sumArr[2]=...
		}
		sumArr.push_back(0.0); //only need to make sumArr the right size to use sumArr[2]=...

		//FINAL SUM CHANNEL
		summaryCanvas->cd(9);

		//TO GET Y VALUES
		for (map<int, map<int, double>>::iterator it = values[charge].begin(); it != values[charge].end(); ++it)
		{
			//CHANNEL DATA
			map<int, double> temp = it->second;

			for (map<int, double>::iterator iti = temp.begin(); iti != temp.end(); ++iti)
			{

				float add = iti->second;
				//	cout << "SUM " << sumArr[iti->first] << "      " << iti->first << "       " << sumArr.size() << endl;

				sumArr[iti->first] = sumArr[iti->first] + add;
			}
		}

		vector<float> x(keys.begin(), keys.end());
		const Int_t n = x.size();

		float scaleBy = sumArr[1];
		for (std::vector<float>::iterator iter = sumArr.begin() + 1; iter != sumArr.end() + 1; ++iter)
		{
			int index = iter - sumArr.begin();
			float toScale = sumArr[index];
			sumArr[index] = toScale / scaleBy;
		}

		TGraphErrors *gr = new TGraphErrors(n, &(x[0]), &(sumArr[1]), 0, 0); //Has to be a 1 cause something is broken
		gr->SetTitle("");
		gr->SetLineColor(4);
		gr->SetLineWidth(1);
		gr->SetMarkerColor(1);
		//gr->SetMarkerColor(4);
		gr->SetMarkerSize(1.0);
		gr->SetMarkerStyle(20);
		gr->GetXaxis()->SetLabelSize(0.035);
		gr->GetYaxis()->SetLabelSize(0.035);
		gr->GetXaxis()->SetTitle("Number of SiPMs");
		gr->GetYaxis()->SetTitle("Light Yield Sum");
		//	legend->AddEntry(gr, ("Position: " + to_string(it->first) + " WOM: " + womIterator->first).data(), "lep");
		gStyle->SetOptFit(1111);

		TF1 *func = new TF1("fit", "pol1");
		func->SetParName(1, "Slope");
		func->SetParName(0, "Offset");
		gr->Fit("fit");
		gr->Draw();

		gPad->Update();
		auto stat = dynamic_cast<TPaveStats *>(gr->FindObject("stats"));
		if (stat)
		{
			stat->SetX1NDC(0.15);
			stat->SetX2NDC(0.45);
			stat->Draw();
		}
		else
		{
			cout << "No stats box found!\n";
		}

		char *valueStr = (char *)malloc(50 * sizeof(char));
		sprintf(valueStr, "%s ", type.c_str());
		TText *t = new TText(.1, .15, valueStr);
		t->SetTextColor(2);
		t->SetTextSize(0.035);
		t->Draw();

		if (charge == 0)
			summaryCanvas->Print((parent_dir + "summaryMaskAnalysis.pdf(").c_str());
		else
			summaryCanvas->Print((parent_dir + "summaryMaskAnalysis.pdf)").c_str());
	}
	return 0;
}