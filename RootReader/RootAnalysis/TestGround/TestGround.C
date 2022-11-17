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
#include <TColor.h>
#include <TGraph2D.h>
#include <TPaletteAxis.h>

namespace fs = std::experimental::filesystem;
using namespace std;
struct stat info;

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

void fitConvolution()
{
	TH1F *h_ExpGauss = new TH1F("h_ExpGauss", "Exponential convoluted by gaussian", 100, 50., 10.);
	for (int i = 0; i < 1e6; i++)
	{
		Double_t x = gRandom->Exp(20); //gives a alpha of -0.3 in the exp
		x += gRandom->Gaus(20, 10);
		h_ExpGauss->Fill(x); //probability density function of the addition of two variables is the convolution of 2 dens. functions
	}
	TF1Convolution *f_conv = new TF1Convolution("expo", "gaus", -1, 10, true);

	f_conv->SetRange(-5., 6.);
	f_conv->SetNofPointsFFT(1000);
	TF1 *f = new TF1("f", *f_conv, -5., 5., f_conv->GetNpar());
	f->SetParameters(1., -0.3, 0., 1.);
	TCanvas *c1 = new TCanvas("c", "c", 800, 1000);
	h_ExpGauss->Fit("f");
	h_ExpGauss->Draw();
	c1->Print("test.pdf");
}

/*******************
__ FUNCTIONS ______
*******************/
float calculateDistance(pair<int, int> p1, pair<int, int> p2)
{
	return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}

int main(int argc, char *argv[])
{
	TApplication *myapp = new TApplication("myapp", 0, 0);

	gStyle->SetOptStat(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);



	TFile *file = new TFile("../rootfiles/58_pos6_angle30_e14_ch32.root");
	if (file->IsZombie())
	{
		cout << "PROBLEM with the initialization of the output ROOT ntuple "
			 << ": check that the path is correct!!!"
			 << endl;

		exit(-1);
	}
	TTree *tree;
	file->GetObject("T", tree);

		int bins=100;
		int maxX=8000;
		int min=-10;

		TCanvas *c = new TCanvas("c", "Graph2D example", 0, 0, 1800, 1000);

		THStack *hs = new THStack("hs","");

		TH1D *allHist = new TH1D("allHist", "", bins, min, maxX);
		allHist->SetFillColorAlpha(kBlue, 0.4);
		allHist->SetLineColorAlpha(1, 0);


		TH1D *allHistErrorP = new TH1D("allHistErrorP", "", bins, min, maxX);
		TH1D *allHistErrorM = new TH1D("allHistErrorM", "", bins, min, maxX);

		TLegend *h_leg = new TLegend(0.50, 0.62, 0.90, 0.90);
		h_leg->SetTextSize(0.03);
	

		
		tree->Draw("Integral[2]>>allHist", "");
		tree->Draw("IntegralErrorP[2]>>allHistErrorP", "");
		tree->Draw("IntegralErrorM[2]>>allHistErrorM", "");


		

		float mean=allHist->GetMean();
		float meanErrorP=allHistErrorP->GetMean();
		float meanErrorM=allHistErrorM->GetMean();
		float cf=1.16818;
		float dcf=0.077463;
		

		h_leg->AddEntry(allHist, Form("default PCS, mean: %1.0f",mean), "f");
		h_leg->AddEntry(allHistErrorP, Form("up-scaled PCS, mean: %1.0f",meanErrorP), "l");
		h_leg->AddEntry(allHistErrorM, Form("down-scaled PCS, mean: %1.0f",meanErrorM), "l");
		h_leg->AddEntry((TObject *)0, Form("CF: %1.2f #pm %1.2f",cf,dcf), "");



		allHistErrorP->SetLineWidth(3);
		allHistErrorM->SetLineWidth(3);
		
		allHistErrorP->SetLineColorAlpha(4, 1);
		allHistErrorM->SetLineColorAlpha(2, 1);

		


 		hs->Add(allHist);
 		hs->Add(allHistErrorP);
 		hs->Add(allHistErrorM);

		hs->Draw("nostack");

		TAxis *yaxisP = hs->GetYaxis();
		TAxis *xaxisP =hs->GetXaxis();
		yaxisP->SetLabelSize(0.04);
		yaxisP->SetTitle("counts");
		yaxisP->SetTitleSize(0.04);
		xaxisP->SetLabelSize(0.04);
		xaxisP->SetTitle("#Lambda [mV #times ns]");
		xaxisP->SetTitleSize(0.04);

		h_leg->Draw();
		c->Print("test.pdf");
		myapp->Run();





	/*TCanvas *c = new TCanvas("c", "Graph2D example", 0, 0, 1000, 1000);
	Double_t x, y, z, P = 6.;
	Int_t np = 200;
	TGraph2D *dt = new TGraph2D();
	TRandom *r = new TRandom();
	int xMin = -400;
	int xMax = 400;
	int yMin = -600;
	int yMax = 600;
	int pixelSize = 5;
	pair<int, int> posWOMD = make_pair(310, -510);

	int numberRows = 1200 / pixelSize;
	int numberColumns = 800 / pixelSize;
	int index = 0;

	TF1 *low = new TF1("fit", "-[0]*exp(x*[2])+[1]", 0, 5000);
	low->SetParameter(0, 3.08772e-04);
	low->SetParameter(1, 9.98528e+01);
	low->SetParameter(2, 9.79445e-03);

	TF1 *mid = new TF1("fit", "-[0]*exp(x*[2])+[1]", 0, 5000);
	mid->SetParameter(0, 7.85506e-04);
	mid->SetParameter(1, 9.99429e+01);
	mid->SetParameter(2, 7.63927e-03);

	TF1 *high = new TF1("fit", "-[0]*exp(x*[2])+[1]", 0, 5000);
	high->SetParameter(0, 1.54281e-08);
	high->SetParameter(1, 9.99753e+01);
	high->SetParameter(2, 1.97434e-02);

	for (Int_t rowIndex = yMin; rowIndex < yMax; rowIndex = rowIndex + pixelSize)
	{

		for (Int_t columnIndex = xMin; columnIndex < xMax; columnIndex = columnIndex + pixelSize)
		{
			pair<int, int> point = make_pair(columnIndex, rowIndex);
			float distance = calculateDistance(posWOMD, point);
			float value = mid->Eval(distance);
			if (value < 99.0)
			{
				value = 99.0;
			}
			dt->SetPoint(index, columnIndex, rowIndex, value);

			index++;
		}
	}

	dt->SetPoint(index, 0, 0, 100.00);
	index++;

	dt->SetPoint(index, 100, 100, 100.00);

	//gStyle->SetPalette(1);

	/*UInt_t Number = 6;
   Double_t Red[Number]    = {1,1,1,1,1,1};
   Double_t Green[Number]  = {1,1,0.7,0.5,0,0};
   Double_t Blue[Number]   = {1,1,0,0,0,0};
   Double_t Length[Number] = {0,0.99000,0.99800,0.99850,0.99870,0.99875};
   

	UInt_t Number = 12;
	Double_t Red[Number] = {1, 1, 1, 0, 0, 0, 0, 0.333, 0.666, 1, 1, 1, 1};
	Double_t Green[Number] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0.666, 0.333, 0, 0};
	Double_t Blue[Number] = {1, 1, 1, 1, 0.666, 0.333, 0, 0, 0, 0, 0, 0, 1};
	Double_t Length[Number] = {0.00000, 0.99000, 0.99002, 0.99100, 0.99200, 0.99300, 0.99400, 0.99500, 0.99600, 0.99700, 0.99800, 0.99825, 0.99850};

	Int_t nb = 30900;
	//	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
	gStyle->SetNumberContours(nb);

	dt->GetXaxis()->SetLabelOffset(0.01);
	dt->GetXaxis()->SetLabelSize(0.1);

	//dt->Draw("TRI2 Z ");
	dt->Draw("TRI2 Z");

	//c->SetTheta(90.0-0.001);
	//c->SetPhi(0.0+0.001);

	dt->Print("test.pdf");

	myapp->Run();

	//
	/*	int bins = 150;
	int xmin = -2;
	int xmax = 16;

	int smooth = 1;
	int n = 100000;
	float labelTextSize = 0.045;

	TCanvas *c1;
	c1 = new TCanvas("c1", "Exponential distribution", 1000, 500);
	TH1F *distribution1;
	TH1F *distribution2;
	TH1F *distribution3;
	TH1F *distribution4;
	distribution1 = new TH1F("", "", bins, xmin, xmax);
	distribution2 = new TH1F("", "", bins, xmin, xmax);
	distribution3 = new TH1F("", "", bins, xmin, xmax);
	distribution4 = new TH1F("", "", bins, xmin, xmax);

	float multi = 0.5;

	TRandom *random1 = new TRandom;
	Double_t data1;
	Double_t data2;
	Double_t data3;
	Double_t data4;

	TF1 *fb2 = new TF1("fa3", "TMath::Poisson(x,[0])", xmin, 50);
	fb2->SetParameter(0, 2);
	//fb2->SetParameter(1, 0.4);

	for (Int_t i = 0; i < n; i++)
	{

		//float randomNumber;
		//randomNumber = ((rand() % 500) + 1)/100.0;
		Double_t rand = fb2->GetRandom();

		data2 = rand+rand*0.1;
		data3 = rand-rand*0.1;
		data1 = rand;
	//	data4 = data1 * (5 * multi);

		distribution1->Fill(data1);
		distribution2->Fill(data2);
		distribution3->Fill(data3);
	//	distribution4->Fill(data4);
	}

	c1->cd(1);
	//fb2->Draw();

	distribution1->Smooth(smooth);
	distribution2->Smooth(smooth);
	distribution3->Smooth(smooth);
	distribution4->Smooth(smooth);
	Int_t color= TColor::GetColor("#30c18f");
	Int_t colorStretched= kRed;
	Int_t colorCompressed= kAzure;




	distribution1->SetFillColorAlpha(color, 0.5);
	//distribution3->SetFillColorAlpha(colorStretched, 0.3);
	//distribution2->SetFillColorAlpha(colorCompressed, 0.3);

	distribution2->SetLineColor(colorStretched);
	distribution1->SetLineColor(1);
	distribution3->SetLineColor(colorCompressed);
	distribution1->SetLineWidth(0);
	distribution2->SetLineWidth(4);
	distribution3->SetLineWidth(4);


	TAxis *yaxis = distribution1->GetYaxis();
	TAxis *xaxis = distribution1->GetXaxis();
	yaxis->SetLabelSize(labelTextSize);
	xaxis->SetLabelSize(labelTextSize);
	yaxis->SetTitle("Counts");
	xaxis->SetLabelSize(labelTextSize);
	xaxis->SetTitleSize(labelTextSize);
	yaxis->SetTitleSize(labelTextSize);
	xaxis->SetTitle(Form("arb.units"));


	TLegend *h_leg = new TLegend(0.65, 0.70, 0.89, 0.88);
	h_leg->SetTextSize(labelTextSize);
	h_leg->SetTextFont(42);
	h_leg->AddEntry(distribution1, Form("unshifted PCS"), "f");
	h_leg->AddEntry(distribution2, Form("up-shifted PCS"), "l");
	h_leg->AddEntry(distribution3, Form("down-shifted PCS"), "l");

















	distribution1->Draw("hist");
	distribution2->Draw("same hist");
	distribution3->Draw("same hist");

	distribution1->GetYaxis()->SetRangeUser(0., 4500);

	h_leg->Draw();

	//	distribution4->Draw("same hist");
	c1->Print("test.pdf");
*/
	//{
	// Generation of an exponential distribution

	/*	fitConvolution();
		gROOT->Reset();
		gROOT->SetStyle("Plain");
		TCanvas *c1;
		TCanvas *c2;
		TH1F *distribution1;
		TH1F *distribution2;
		TH1F *distribution3;
		TH1F *distribution4;
		TF1 *exp1;

		c1 = new TCanvas("c1", "Exponential distribution", 200, 10, 600, 800);
		//   c1->Divide();
		c2 = new TCanvas("c2", "Exponential distribution", 200, 10, 600, 800);
		c2->Divide(1, 3);

		// Create some histograms.
		distribution1 = new TH1F("Distribution", "Distribution", 100, 0, 1);
		distribution2 = new TH1F("Distribution", "Distribution", 100, 0, 1);
		Int_t m = 100;
		Double_t a = -15;
		Double_t b = 15;
		distribution3 = new TH1F("Exp-Distribution", "Exp-Distribution", m, a, b);
		distribution4 = new TH1F("Exp-Distribution1", "Exp-Distribution2", m, a, b);

		// Use two pseudo-random number generator
		TRandom *random1 = new TRandom;
		TRandom *random2 = new TRandom;
		TRandom *random3 = new TRandom;
		// Random seeds for the random generators are set
		// These are the starting values for the algorithm
		// producing the pseudo-random numbers
		// seed=0: actual computer clock time in seconds
		random1->SetSeed(0);
		random2->SetSeed(5000);
		random3->SetSeed(1);
		// number of generated numbers
		int n = 50000;

		Double_t expo1 = 0;
		Double_t x1 = 0;
		Double_t x2 = 0;
		Double_t lambda = 0.75;
		Double_t Pi = 3.1415;
		Double_t sigma = 1.0;

		Double_t max = lambda * exp(0);
		//  cout << "max: " << max  << endl;

		// Weighting factor to obtain a normalized distribution
		int success = 0;

		Double_t data1;
		Double_t data2;
		Double_t data3;

		Double_t par[2];
		exp1 = new TF1("exp1", "expo", 8, 10);

		for (Int_t i = 0; i < n; i++)
		{
			// Produce n-times uniform random numbers data1 and data2
			data1 = random1->Rndm();
			data2 = random2->Rndm();
			data3 = random3->Rndm();

			x1 = sqrt(-2 * log(data2)) * cos(2 * Pi * data3);
			expo1 = (-1 / (lambda)) * log(data1);
			distribution1->Fill(data1);
			distribution3->Fill(expo1);
			expo1 = expo1 + x1 * sigma;
			distribution4->Fill(expo1);
		}

		c2->cd(1);
		distribution3->Draw();
		c2->cd(2);
		distribution4->Draw("e1p");

		// Plot next figures into canvas c2
		TPad *pad1 = new TPad("pad1", "pad", 0.00, 0.00, 1.00, 1.00);
		c2->cd(3);
		pad1->Draw();
		pad1->SetLogy(1);
		pad1->cd();
		distribution4->Fit("exp1", "R");
		distribution4->Draw("e1p");
		c2->Update();

		// Produce a eps-file
		c2->Print("exponentialgaussian.gif", "gif");
	}
*/

	return 0;
}