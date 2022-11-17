#include <TApplication.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2C.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPaveLabel.h>
#include <TGraphErrors.h>

#include <TTreeViewer.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <windows.h>
#include <stdio.h>
#include <filesystem>
#include <experimental/filesystem>
#include <iostream>
#include <numeric>  
#include <sys/types.h>
#include <sys/stat.h>

struct stat info;

/**
Skript um die Mean Werte von tSiPM zu bestimmen (T_trigger-t_signal)
Man muss immer angeben welche Channel, weil tSiPM immer für 2 WOMS ist
Die Histogramme von tSipM werden gefittet bis zu "percentage" einer Gaußglocke
Manuell zusammenfügen der PDFs ist notwendig.
Es nimmt alle ROOTFILES im Ordner rootfiles und erzeugt Dateien im Ordner "Plots"
Es muss richtig benannt werden, sonst funktioniert die automatische Legende nicht


SETUP: ROOT VISUAL STUDIO runterladen und installieren-> Umgebungsvariablen setzen: 
ROOTSYS: D:\root\root_v6.16.00
PATH: D:\root\bin
LD_LIBRARY_PATH_  D:\root\lib oder so

Im Projekt in VS:
Linker->Eingabge->Zusätzliche Abhängigkeiten: D:\root\lib\*.lib;
VC++ Verzeichnisse-> Include/Libbibliothek jeweils: D:\root\include; bzw. D:\root\lib;

**/



//https://www.pp.rhul.ac.uk/~cowan/stat/root/hist/rootTest.cc
using namespace std;
namespace fs = std::experimental::filesystem;



float percentage = 0.5f;
float referenceChannel = 0;
vector<int> useChannels = { 7,8,9,10,11,12,13,14 };
int channelOffset = 7;  //FOR XAXIS TICK Labels
int channelNr = useChannels.size();
Double_t rangeMin = -0.5; //Mit OFFSET!!! Aus Channel 7 wird channel 0 -> hier -0.5 eintrage
Double_t rangeMax = 7.5;


void printi(std::vector<float> const& input)
{
	std::copy(input.begin(),
		input.end(),
		std::ostream_iterator<float>(std::cout, " "));
}

vector<string> split(const string& str, const string& delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}

int main(int argc, char** argv)

{
	
	TApplication* theApp = new TApplication("App", &argc, argv);

	system("cls"); //Remove ROOT Warnings
	char buf[256];
	GetCurrentDirectoryA(256, buf);
	string dir = string(buf) + '\\';
	string rootfolder = "rootfiles\\";
	string savefolder = "Plots\\";
	string muonDir = dir + rootfolder+"muon\\";	
	string pionDir = dir + rootfolder + "pion\\";
	string electronDir = dir + rootfolder + "electron\\";

	string generalPrintDir = dir + savefolder ;

	string muonPrintDir = dir + savefolder + "muon\\";
	string pionPrintDir = dir + savefolder + "pion\\";
	string electronPrintDir = dir + savefolder + "electron\\";



	vector<string> paths = {};
	vector<string> pathNames = {};


	if (stat(muonDir.c_str(), &info) != 0)
		printf("cannot access %s\n", muonDir.c_str());
	else if (info.st_mode & S_IFDIR)
	{
		paths.push_back(muonDir);
		pathNames.push_back("muon");
	}
	else
		printf("%s is no directory\n", muonDir.c_str());

	if (stat(pionDir.c_str(), &info) != 0)
		printf("cannot access %s\n", pionDir.c_str());
	else if (info.st_mode & S_IFDIR)
	{
		paths.push_back(pionDir);
		pathNames.push_back("pion");

	}
	else
		printf("%s is no directory\n", pionDir.c_str());

	if (stat(electronDir.c_str(), &info) != 0)
		printf("cannot access %s\n", electronDir.c_str());
	else if (info.st_mode & S_IFDIR) {
		paths.push_back(electronDir);
		pathNames.push_back("electron");

	}
	else
		printf("%s is no directory\n", electronDir.c_str());


	if (!fs::is_directory(generalPrintDir) || !fs::exists(generalPrintDir)) { // Check if src folder exists
		fs::create_directory(generalPrintDir); // create src folder
	}
	if (!fs::is_directory(muonPrintDir) || !fs::exists(muonPrintDir)) { // Check if src folder exists
		fs::create_directory(muonPrintDir); // create src folder
	}
	if (!fs::is_directory(pionPrintDir) || !fs::exists(pionPrintDir)) { // Check if src folder exists
		fs::create_directory(pionPrintDir); // create src folder
	}
	if (!fs::is_directory(electronPrintDir) || !fs::exists(electronPrintDir)) { // Check if src folder exists
		fs::create_directory(electronPrintDir); // create src folder
	}











	map<string,map<int, vector<float>>> womMap; //Key: WOMNAME: A,B,C,D..., Value: Mean Value Vector for all Channels for the different runs
	map<string, map<int, vector<float>>> womMapSpecialShifted; //Key: WOMNAME: A,B,C,D..., Value: Mean Value Vector for all Channels for the different runs  SHIFTED um Channel 0 Wert

	map<string, map<string, vector<float>>> finalMeanChannelMapY; //KEY: WOM < Key=PARTICLE TYPE Value=mittlere Kanäle>>>
	map<string, map<string, vector<float>>> finalMeanChannelMapX; //KEY: WOM < Key=PARTICLE TYPE Value=X>>>
	map<string, map<string, vector<float>>> finalMeanChannelMapYError; //KEY: WOM < Key=PARTICLE TYPE Value=yError>>>


	for (std::vector<string>::iterator iti = paths.begin(); iti != paths.end(); ++iti) {

		//Overview
		//Show Mean for all Channels for different ROOT Files
		map<int, vector<float>>meanMap; //Key: RunNR=0,1,2,3,4,..., Value: Mean Value Vector for all Channels
		map<int, vector<float>>meanErrorMap; //Key: RunNR=0,1,2,3,4,..., Value: Mean Value Errors Vector for all Channels

		map<string, map<int, vector<float>>> WOMRunMap; //Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>
		map<string, map<int, vector<float>>> WOMRunMapErrors; //Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>  Gauß Errors


		int index = distance(paths.begin(), iti);
		std::cout << *iti;
		string path = *iti;

		int counter = 0;
		float channels[8] = {};

		//PRINT OR NOT PRINT
		gROOT->SetBatch(true);

		TCanvas* c1;

		vector<string> posNames, womNames;
		for (const auto& entry : fs::directory_iterator(path)) {
			float stdFits[8] = {};
			vector<float> meanFits;
			vector<float> meanErrors;

			vector<float> meanFitsShifted;		//All Shifted um Pol0 Fit durch Werte
			vector<float> meanFitsShiftedSpecialChannel; //ALL Shifted um Channel 0 Wert

			float chi2Fits[8] = {};


			std::cout << entry.path() << std::endl;
			string file = entry.path().u8string();

			vector<string> namecontent = split(file, "_"); //ONLY FOR NAMESHEME: 65_muon6_pos5_AB
			string posData = namecontent[2];
			string pos = posData.erase(0, 3); //EXTRACT 5 for legend
			string wom = namecontent[3].erase(0, 1);
			wom = wom.erase(1, wom.length());

			posNames.push_back(pos);
			womNames.push_back(wom);

			cout << pos << wom;
			Double_t referenceError;
			char const* c = file.data();
			TFile* f = new TFile(c);
			TH1F* h;
			if (f != nullptr) {
				TTree* tree = (TTree*)f->Get("T");
				c1 = new TCanvas("c", "tSiPM", 1080, 1080);
				c1->SetTitle("TTTT");

				c1->Divide(3, 3);
				for (int a = 1; a <= channelNr; a++) {
					c1->cd(a);
					h = new TH1F("h", "", 1000, -85, -15);
					//h->SetStats(false);

					int selectedChannel = useChannels[a - 1];

					tree->Draw(("tSiPM[" + to_string(selectedChannel) + "]>>h").c_str());

					gStyle->SetOptFit();


					int binmax = h->GetMaximumBin();
					float xbinmax = h->GetXaxis()->GetBinCenter(binmax);
					float bin1 = h->FindFirstBinAbove(h->GetMaximum() * percentage);
					float bin2 = h->FindLastBinAbove(h->GetMaximum() * percentage);
					double fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
					float from = h->GetXaxis()->GetBinCenter(bin1);
					float to = h->GetXaxis()->GetBinCenter(bin2);





					h->Draw();
					h->Fit("gaus", "", "", from, to);



					TF1* fit = h->GetFunction("gaus");
					Double_t p1 = fit->GetParameter(2);
					Double_t mean = fit->GetParameter(1);	

					Double_t error = fit->GetParError(1); //wenn die Differenzen angezeigt werden (z.B. zu Referenzkanal), dann nimmt man hier auch die Fehlerfortpflanzung

					if (a - 1 == 0) {
						referenceError = error;
					}
					else {
						double pot1 = pow(error, 2);
						double pot2 = pow(referenceError, 2);
						error = sqrt(pot1 + pot2);
					}





					printf("TTTTTT: %g und %g und %g \n", fwhm, from, p1);

					TLegend* leg = new TLegend(0.1, 0.7, 0.48, 0.9);
					leg->AddEntry(("e", " ", "Std: " + to_string(p1)).c_str());
					//leg->Draw();

					stdFits[a - 1] = p1;
					meanFits.push_back(mean);
					meanErrors.push_back(error);
					channels[a - 1] = selectedChannel - channelOffset;
					chi2Fits[a - 1] = (fit->GetChisquare()) / (fit->GetNDF());
					c1->Update();
				}
				c1->cd();
				TPaveLabel* title = new TPaveLabel(.85, .01, 0.99, .07, ("POS: " + pos + " WOM: " + wom).data(), "brndc");
				title->Draw();
				//	c1->Close();
				TCanvas* c2 = new TCanvas("c2", "multigraph", 1000, 1000);
				gPad->Modified();
				gPad->Update();
				c2->Divide(1, 2);

				TGraph* gr = new TGraph(channelNr, channels, meanFits.data());
				TGraph* grchi2 = new TGraph(channelNr, channels, chi2Fits);
				gr->SetTitle("Mean per Channel");
				//	gr->SetLineColor(0);
				gr->SetLineWidth(0);
				gr->SetMarkerColor(4);
				gr->SetMarkerSize(1.5);
				gr->SetMarkerStyle(8);
				gr->GetXaxis()->SetTitle("Channel");
				gr->GetYaxis()->SetTitle("Mean");

				grchi2->SetTitle("Chi^{2}/NDF per Channel");
				//grchi2->SetLineColor(0);
				grchi2->SetLineWidth(0);
				grchi2->SetMarkerColor(4);
				grchi2->SetMarkerSize(1.5);
				grchi2->SetMarkerStyle(8);
				grchi2->GetXaxis()->SetTitle("Channel");
				grchi2->GetYaxis()->SetTitle("Chi^{2}/NDF");


				c2->cd(1);

				gr->Draw();
				gr->Fit("pol0");

				gPad->SetGrid(1, 1);
				gPad->Modified();
				gPad->Update();

				c2->cd(2);

				grchi2->Draw();
				gPad->SetGrid(1, 1);
				gPad->Modified();
				gPad->Update();
				c2->cd();
				title->Draw();



				meanMap[counter] = meanFits;
				meanErrorMap[counter] = meanErrors;

				map<int, vector<float>> tempMap = womMap[wom];
				map<int, vector<float>> tempMapSpecialShifted = womMapSpecialShifted[wom];

				TF1* gfit = gr->GetFunction("pol0");
					
				//SHIFT-> Damit man Mean werte für verschiedene Positionen miteinander vergleichen kann, müssen die alle erstmal auf das selbe Niveau gebracht werden-> Um constante nach unten verschieben
				int cc = 0;
				float reference = 0; //Channel 0 Value
				
				for (std::vector<float>::iterator it = meanFits.begin(); it != meanFits.end(); ++it) {
					reference = meanFits[referenceChannel];
					meanFitsShifted.push_back(meanFits[cc] - (gfit->GetParameter(0)));
					meanFitsShiftedSpecialChannel.push_back(meanFits[cc] -reference );

					cc++;
				}

				
				tempMap[counter] = meanFitsShifted;
				tempMapSpecialShifted[counter] = meanFitsShiftedSpecialChannel;

				
				
				womMap[wom] = tempMap;
				womMapSpecialShifted[wom] = tempMapSpecialShifted;


				map<int, vector<float>> runs = WOMRunMap[wom];
				runs[stoi(pos)] = meanFitsShiftedSpecialChannel;
				WOMRunMap[wom] = runs;


				map<int, vector<float>> runsError = WOMRunMapErrors[wom];
				runsError[stoi(pos)] = meanErrors;
				WOMRunMapErrors[wom] = runsError;


				string savefolderType = savefolder + pathNames[index];

				string saveHist = savefolderType + "//histogram_" + to_string(counter) + "_POS" + pos + "_WOM" + wom + ".pdf";
				string saveResult = savefolderType + "//results_" + to_string(counter) + "_POS" + pos + "_WOM" + wom + ".pdf";
				c1->SaveAs(saveHist.data());

				c2->SaveAs(saveResult.data());



				counter++;





			}

			gPad->Modified();
			gPad->Update();

		}

		//SAVE OVERVIEW

		//gROOT->SetBatch(false);
		cout << "\n";
		map<int, vector<float>>::iterator temp_iti = meanMap.begin();
		map<int, vector<float>>::iterator temp_iti_error = meanErrorMap.begin();

		int tempi = 0;
		while (temp_iti != meanMap.end())
		{
			string pos = posNames[tempi];
			string wom = womNames[tempi];
			cout << ("Position: " + pos + " WOM: " + wom +" referenceChannel: "+to_string(referenceChannel) ).data() << endl;
			printi(temp_iti->second);
			cout << "\n" << endl;

			printi(temp_iti_error->second);

			cout << "\n" << endl;

			printi(womMapSpecialShifted[wom][temp_iti->first]);
			cout << "\n" << endl;

			tempi++;
			temp_iti++;
			temp_iti_error++;

		}
		cout << "\n" << endl;













		TCanvas* c3 = new TCanvas("c3", "overview", 1500, 1200);
		//c3->DrawFrame(-10, -10, 10, 10);
		//c3->SetRightMargin(0.3);

		c3->cd();
		TMultiGraph* mg = new TMultiGraph();

		gPad->Modified();
		gPad->Update();

		auto legend = new TLegend(0.1, 0.91, 0.5, 0.99);
		legend->SetNColumns(3);
		legend->SetHeader(("Arrival Time Analysis: " + pathNames[index]).data(), "L"); // option "C" allows to center the header
		map<int, vector<float>>::iterator it = meanMap.begin();
		int counti = 0;
		while (it != meanMap.end())
		{

			int co = counti + 1;
			if (co == 5)co = 12;
			TGraph* gr = new TGraph(channelNr, channels, it->second.data());
			gr->SetTitle("Mean per Channel");
			gr->SetLineColor(co);
			gr->SetLineWidth(1);
			gr->SetMarkerColor(co);
			//gr->SetMarkerColor(4);
			gr->SetMarkerSize(1.5);
			gr->SetMarkerStyle(20);
			gr->GetXaxis()->SetLabelSize(0.035);
			gr->GetYaxis()->SetLabelSize(0.035);

			//gr->Draw();
			//gPad->Modified();
			//gPad->Update();

			string pos = posNames[counti];
			string wom = womNames[counti];

			legend->AddEntry(gr, ("Position: " + pos + " WOM: " + wom).data(), "lep");

			gStyle->SetOptFit(0000);
			gr->Fit("pol0");
			gr->GetFunction("pol0")->SetLineColor(co);
			mg->Add(gr);

			counti++;
			it++;
		}
		gStyle->SetLegendBorderSize(1);
		legend->SetMargin(0.25);

		gStyle->SetLegendTextSize(0.015);
		//mg->Draw("c p A pmc plc");
		mg->GetXaxis()->SetTitle("Channel");
		mg->GetYaxis()->SetTitle("Mean of tSiPM");
		mg->GetYaxis()->SetTitleSize(0.02);
		mg->GetXaxis()->SetTitleSize(0.02);

		mg->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
		mg->GetXaxis()->SetTickLength(0.015);
		mg->GetYaxis()->SetTickLength(0.015);
		mg->GetXaxis()->SetLabelSize(0.018);
		mg->GetYaxis()->SetLabelSize(0.018);
		mg->Draw("c p A");
		legend->Draw();

		gPad->SetGrid(1, 1);
		gPad->Modified();
		gPad->Update();

		string savefolderType = savefolder + pathNames[index];

		string saveFinal = savefolderType + "//overview.pdf";
		c3->SaveAs(saveFinal.data());



		//Abweichungen: IDEE: Die Verläufe der Mean Kurven für verschiedene Strahlpositionen sind für die selben WOMs relativ gleich. Hier geht es jetzt darum zu untersuchen, wie die Abweichung von
		//einem Mittelwert ist. Dazu filtert man die Runs/MeanMap nach WOM, bildet dann für jeden Channel für jeden Wom den Mittelwert. Die Fehlerbalken sind dann entweder Standardabweichung oder Min/Max
		//gROOT->SetBatch(false);

		map<string, map<int, vector<float>>>::iterator itu = womMap.begin();
		TCanvas* c10= new TCanvas("c10", "Overview", 1080, 1080); // problems here
		c10->Divide(1, 2, 1E-11, 1E-11);
		c10->cd(1);

		int womcount = 1;
		while (itu != womMap.end()) {
			//WOMS
			map<int, vector<float>> channelMap; //KEY: CHANNEL, Value: Means for Runs for WOM
			map<int, vector<float>> temp = itu->second; //RUNS FOR ONE WOM
			map<int, vector<float>>::iterator iti = temp.begin();
			int c = 0;
			while (iti != temp.end())
			{
				//ALL RUNS
				vector<float> channelData = iti->second; //0,1,2,3,4,5,6,7 Means for Channels
				int d = 0;

				for (auto i : channelData) {
					channelMap[d].push_back(i);
					d++;
				}
				c++;	//INDEX
				iti++; //ITERATOR
			}

			//MEAN FOR CHANNELMAPS
			map<int, vector<float>>::iterator itis = channelMap.begin();
			int e = 0;

			//3 Lists for Mean, Distance Min, Distance Max
			vector<double> meansForChannels;
			vector<double> minsForChannels;
			vector<double> maxsForChannels;
			vector<double> yError;

			vector<double> x = { 0,1,2,3,4,5,6,7 };


			
			while (itis != channelMap.end())  //FOR CHANNELS
			{
				float average = accumulate(itis->second.begin(), itis->second.end(), 0.0) / itis->second.size();
				
				
				auto max_ = max_element(begin(itis->second), end(itis->second)); // c++11
				auto min_ = min_element(begin(itis->second), end(itis->second)); // c++11

				float max = *max_;
				float min = *min_;

				meansForChannels.push_back(average);
				minsForChannels.push_back(min);
				maxsForChannels.push_back(max);

				float stdDeviation = 0.0, variance = 0.0, sum = 0.0, mean = 0.0;
				int i;
				for (i = 0; i < itis->second.size(); ++i)
					sum += itis->second[i];
				mean = sum / itis->second.size();
				for (i = 0; i < itis->second.size(); ++i)
					variance += pow(itis->second[i] - mean, 2);
				variance = variance / itis->second.size();
				stdDeviation = sqrt(variance);


				yError.push_back(stdDeviation); //USE AS Y ERROR STDEVIATION
				itis++;
				e++;
			}
			c10->cd(womcount);
			
		//	gROOT->SetBatch(false);
			TCanvas* c5 = new TCanvas("c5", "Fluctuation", 1080, 1080);
			c5->cd();
			auto legend = new TLegend(0.1, 0.87, 0.4, 0.99);

			c5->SetGrid();
			gStyle->SetOptFit(1111);
			const Int_t n = meansForChannels.size();
			vector<double> xErr = vector<double>(n);
			auto test = &(x[0]);
			TGraphErrors* gr = new TGraphErrors(n, &(x[0]), &(meansForChannels[0]), 0, &(yError[0]));

			float yErrorMean = accumulate(yError.begin(), yError.end(), 0.0) / yError.size();
			float yMean = accumulate(meansForChannels.begin(), meansForChannels.end(), 0.0) / meansForChannels.size();
			legend->SetHeader(("Fluctuations for WOM: " + itu->first + " and " + pathNames[index]).data(), "L"); // option "C" allows to center the header

			legend->AddEntry((TObject*)0, ("Mean: " + to_string(yMean)).data(), "");
			legend->AddEntry((TObject*)0, ("Std Deviation Mean: " + to_string(yErrorMean)).data(), "");
			gr->GetXaxis()->SetRangeUser(-0.1, 7.1);
			gr->GetXaxis()->SetRangeUser(-0.1, 7.1);
			gr->SetLineColor(4);
			gr->SetMarkerSize(1.5);
			gr->SetMarkerColor(1);
			gr->SetMarkerStyle(20);
			gStyle->SetLegendBorderSize(1.0);
			legend->SetMargin(0.25);
			gr->SetTitle("");
			gStyle->SetLegendTextSize(0.025);
			//mg->Draw("c p A pmc plc");
			gr->GetXaxis()->SetTitle("Channel");
			gr->GetYaxis()->SetTitle("Mean for Runs");
			gr->GetYaxis()->SetTitleSize(0.02);
			gr->GetXaxis()->SetTitleSize(0.02);

			gr->GetXaxis()->SetTickLength(0.015);
			gr->GetYaxis()->SetTickLength(0.015);
			gr->GetXaxis()->SetLabelSize(0.018);
			gr->GetYaxis()->SetLabelSize(0.018);
			gr->Draw("c p A");
			legend->Draw();

			

			string savefolderType = savefolder + pathNames[index];

		

			c10->cd(womcount);
			c5->DrawClonePad();
			

			womcount++;
			itu++;
		}
		c10->Draw();
		saveFinal = savefolderType + "//fluctuations.pdf";
		c10->SaveAs(saveFinal.data());
	








		/*
		Mit den Fluctuations hatte man gesehen, wie stark die Daten schwanken gemittelt über die Runs, jetzt:
		Wie stark schwankt ein Kanal zum anderen je nach Run.
		Dazu: Mean Kurven für jeden Run nach unten verschoben um ReferenzKanalwert (z.B. Channel 0 )

		
		*/
		gROOT->SetBatch(true);
		map<string, map<int, vector<float>>>::iterator womIter =WOMRunMap.begin();//Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>
		 womcount = 1;
		 TCanvas* c11 = new TCanvas("c11", "Overview", 1080, 1080); // problems here
		 c11->Divide(1, 2, 1E-11, 1E-11);
		 vector<float> x = { 0,1,2,3,4,5,6,7 };

		while (womIter != WOMRunMap.end()) {

			c11->cd(womcount);
			TMultiGraph* mg = new TMultiGraph();

			gPad->Modified();
			gPad->Update();

			auto legend = new TLegend(0.1, 0.91, 0.5, 0.99);
			legend->SetNColumns(3);
			legend->SetHeader(("Arrival Time Analysis: " + pathNames[index]).data(), "L"); // option "C" allows to center the header
			map<int, vector<float>>::iterator it = womIter->second.begin();
			int counti = 0;
			while (it != womIter->second.end())
			{

				int co = counti + 1;
				if (co == 5)co = 12;
				
				vector<float> channeldata = (it->second);
				const Int_t n = channeldata.size();
				vector<float> yError=WOMRunMapErrors[womIter->first][it->first];

				TGraphErrors* gr = new TGraphErrors(n, &(x[0]), &(channeldata[0]), 0, &(yError[0]));
				gr->SetTitle("Mean per Channel");
				gr->SetLineColor(co);
				gr->SetLineWidth(1);
				gr->SetMarkerColor(co);
				//gr->SetMarkerColor(4);
				gr->SetMarkerSize(1.5);
				gr->SetMarkerStyle(20);
				gr->GetXaxis()->SetLabelSize(0.035);
				gr->GetYaxis()->SetLabelSize(0.035);

				legend->AddEntry(gr, ("Position: " + to_string(it->first) + " WOM: " + womIter->first).data(), "lep");

				gStyle->SetOptFit(0000);
		
				mg->Add(gr);

				counti++;
				it++;
			}
			gStyle->SetLegendBorderSize(1);
			legend->SetMargin(0.25);

			gStyle->SetLegendTextSize(0.02);
			//mg->Draw("c p A pmc plc");
			mg->GetXaxis()->SetTitle("Channel");
			mg->GetYaxis()->SetTitle("Mean of tSiPM");
			mg->GetYaxis()->SetTitleSize(0.02);
			mg->GetXaxis()->SetTitleSize(0.02);

			mg->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
			mg->GetXaxis()->SetTickLength(0.015);
			mg->GetYaxis()->SetTickLength(0.015);
			mg->GetXaxis()->SetLabelSize(0.018);
			mg->GetYaxis()->SetLabelSize(0.018);
			mg->Draw("c p A");
			legend->Draw();

			gPad->SetGrid(1, 1);
			gPad->Modified();
			gPad->Update();

			womcount++;
			womIter++;

		}
		saveFinal = savefolderType + "//runFluctuations.pdf";
		c11->SaveAs(saveFinal.data());





		/*
		Jetzt: für jeden Run: gewichtetes Mittel der Kanäle bilden und diese Mittel zu Strahlpositionen auftragen für jeden WOM
		*/
		gROOT->SetBatch(true);
		map<string, map<int, vector<float>>>::iterator womIterator = WOMRunMap.begin();//Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>
		womcount = 1;
		
		TCanvas* c12 = new TCanvas("c12", "Overview", 1080, 1080); // problems here
		c12->Divide(1, 2, 1E-11, 1E-11);

		while (womIterator != WOMRunMap.end()) {

			c12->cd(womcount);
			TMultiGraph* mg = new TMultiGraph();

			gPad->Modified();
			gPad->Update();

			auto legend = new TLegend(0.1, 0.91, 0.5, 0.99);
			legend->SetNColumns(3);
			string title = "Particle: " + pathNames[index] + " WOM:  " + womIterator->first;
			legend->SetHeader(title.c_str(), "L"); // option "C" allows to center the header
			map<int, vector<float>>::iterator it = womIterator->second.begin();


			vector<float> runs;
			vector<float> meanDeviations;
			vector<float> standardDeviation;

			while (it != womIterator->second.end())
			{ //OVER RUNS

				vector<float> abweichung = (it->second);
				float abweichungSumme = accumulate(abweichung.begin() + 1, abweichung.end(), 0.0);
				//gewichteter Mittelwert des Kanals
				float gewichteterMittelwert = 0.0;
				float gewichteterFehler = 0.0; //http://people.physik.hu-berlin.de/~julien/sub/Kurzeinfuerung%20-%20Fehlerrechnung.pdf  Problem: Gewichte haben Fehler-> komplizierte Fehlerfortpflanzung

				for (int j = 1; j < abweichung.size(); ++j) {
					float gewicht = abweichung[j];
					float unsicherheit = WOMRunMapErrors[womIterator->first][it->first][j]; //Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>  Gauß Errors

					gewichteterMittelwert = gewichteterMittelwert + (gewicht * j);
				}
				for (int j = 1; j < abweichung.size(); ++j) {
					float gewicht = abweichung[j];
					float unsicherheit = WOMRunMapErrors[womIterator->first][it->first][j]; //Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>  Gauß Errors

					float zaehlerAbleitung = ((j * abweichungSumme) - gewichteterMittelwert);
					float nennerAbleitung = pow(abweichungSumme, 2);
					float zusammen =pow( (zaehlerAbleitung / nennerAbleitung) * unsicherheit,2); //Aus Gauß Fehlerfortpflanzung
					gewichteterFehler = gewichteterFehler + zusammen;

				}



				gewichteterMittelwert = (gewichteterMittelwert / abweichungSumme);
				gewichteterFehler = sqrt(gewichteterFehler);

				float y = gewichteterMittelwert; //DROP Channel 0
				float yError = gewichteterFehler;
			
				
				

				runs.push_back(it->first);
				meanDeviations.push_back(y);
				standardDeviation.push_back(yError);
				
				it++;
			}
			



			vector<float> yError = standardDeviation;
			vector<float> x = runs;
			vector<float> y = meanDeviations;
			const Int_t n = x.size();


			finalMeanChannelMapY[womIterator->first][pathNames[index]]=y;
			finalMeanChannelMapX[womIterator->first][pathNames[index]] = x;
			finalMeanChannelMapYError[womIterator->first][pathNames[index]] = yError;

			TGraphErrors* gr = new TGraphErrors(n, &(x[0]), &(y[0]), 0, &(yError[0]));
				
				
				gr->SetTitle("Mean per Channel");
				gr->SetLineColor(4);
				gr->SetLineWidth(1);
				gr->SetMarkerColor(1);
				//gr->SetMarkerColor(4);
				gr->SetMarkerSize(1.0);
				gr->SetMarkerStyle(20);
				gr->GetXaxis()->SetLabelSize(0.035);
				gr->GetYaxis()->SetLabelSize(0.035);

			//	legend->AddEntry(gr, ("Position: " + to_string(it->first) + " WOM: " + womIterator->first).data(), "lep");

				gStyle->SetOptFit(0000);

				mg->Add(gr);

			
			gStyle->SetLegendBorderSize(1);
			legend->SetMargin(0.25);

			gStyle->SetLegendTextSize(0.02);
			//mg->Draw("c p A pmc plc");
			mg->GetXaxis()->SetTitle("Beam Position");
			mg->GetYaxis()->SetTitle("Weighted channel mean");
			mg->GetYaxis()->SetTitleSize(0.02);
			mg->GetXaxis()->SetTitleSize(0.02);

			mg->GetXaxis()->SetRangeUser(0, 15);
			mg->GetYaxis()->SetRangeUser(2, 9);

			mg->GetXaxis()->SetTickLength(0.015);
			mg->GetYaxis()->SetTickLength(0.015);
			mg->GetXaxis()->SetLabelSize(0.018);
			mg->GetYaxis()->SetLabelSize(0.018);
			mg->Draw("c p A");
			legend->Draw();

			gPad->SetGrid(1, 1);
			gPad->Modified();
			gPad->Update();

			womcount++;
			womIterator++;

		}
		saveFinal = savefolderType + "//runMeanFluctuations.pdf";
		c12->SaveAs(saveFinal.data());




}

		gROOT->SetBatch(false);



		int womcount = 1;

		map<string, map<string, vector<float>>>::iterator womIterator =finalMeanChannelMapY.begin();//Key: WOMNAME, VALUE: map<Runnummer,Channeldaten>

		TCanvas* c20 = new TCanvas("c20", "FinalOverview", 1080, 1080); // problems here
		c20->cd();
		c20->Divide(1, 2, 1E-11, 1E-11);

	

		

		while (womIterator != finalMeanChannelMapY.end())
		{ //OVER WOMS
			c20->cd(womcount);
			TMultiGraph* mg = new TMultiGraph();

			auto legend = new TLegend(0.1, 0.91, 0.5, 0.99);
			legend->SetNColumns(3);
			string title = " WOM:  " + womIterator->first;
			legend->SetHeader(title.c_str(), "L"); // option "C" allows to center the header
			legend->SetMargin(0.25);





			int color = 1;
			map<string, vector<float>>::iterator particleIterator = womIterator->second.begin();
			while (particleIterator != womIterator->second.end()) {
				vector<float> y = particleIterator->second;
				vector<float> x = finalMeanChannelMapX[womIterator->first][particleIterator->first];
				vector<float> yError = finalMeanChannelMapYError[womIterator->first][particleIterator->first];;
				const Int_t n = x.size();

				TGraphErrors* gr = new TGraphErrors(n, &(x[0]), &(y[0]), 0, &(yError[0]));
				gr->SetTitle("Mean per Channel");
				gr->SetLineColor(color);
				gr->SetLineWidth(1);
				gr->SetMarkerColor(color);
				//gr->SetMarkerColor(4);
				gr->SetMarkerSize(1.0);
				gr->SetMarkerStyle(20);
				gr->GetXaxis()->SetLabelSize(0.035);
				gr->GetYaxis()->SetLabelSize(0.035);
				gStyle->SetOptFit(0000);

				legend->AddEntry(gr, ("Particle: " + particleIterator->first).data(), "lep");

				mg->Add(gr);

				gPad->Modified();
				gPad->Update();

				particleIterator++;
				color++;
			}

			

			gStyle->SetLegendBorderSize(1);

			gStyle->SetLegendTextSize(0.02);
			mg->GetXaxis()->SetTitle("Beam Position");
			mg->GetYaxis()->SetTitle("Weighted channel mean");
			mg->GetYaxis()->SetTitleSize(0.02);
			mg->GetXaxis()->SetTitleSize(0.02);

			mg->GetXaxis()->SetRangeUser(0, 15);
			mg->GetYaxis()->SetRangeUser(2, 9);

			mg->GetXaxis()->SetTickLength(0.015);
			mg->GetYaxis()->SetTickLength(0.015);
			mg->GetXaxis()->SetLabelSize(0.018);
			mg->GetYaxis()->SetLabelSize(0.018);
			mg->Draw("c p A");
			legend->Draw();

			//mg->Draw("c p A pmc plc");

			gPad->SetGrid(1, 1);
			gPad->Modified();
			gPad->Update();


			womcount++;
			womIterator++;
		}

		string saveFinal = savefolder + "//meanChannelOverview.pdf";
		c20->SaveAs(saveFinal.data());
	





	

	theApp->Run();

	return 0;


}

