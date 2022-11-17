//Including root functionalities:
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPad.h>
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
#include <TNtuple.h>
#include <TImage.h>
#include <TAttImage.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <numeric>
#include <tuple>
#include <map>

using namespace std;

void plasticScintAmps(TTree* tree, string rootFilename, string saveFolder, float* positionInfo);

void plasticScintAmps(TTree* tree, string rootFilename, string saveFolder, float* positionInfo) {

	int runNumber = positionInfo[0];
	float posTopLeftLimit = positionInfo[1];
	float posTopRightLimit = positionInfo[2];
	float posBotLeftLimit = positionInfo[3];
	float posBotRightLimit  = positionInfo[4];
		
	//Style Settings:
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);
  gStyle->SetLineScalePS(1);

	// gStyle->SetOptStat(210);
	// gStyle->SetStatFontSize(0.32);


  std::string locationAmps = saveFolder + "/plasticScintAmps.pdf";

/**
	TFile file(rootFilename.c_str());
	if (file.IsZombie())
    {
        std::cout << "Problem with file " << rootFilename << " in plasticScintPlot; check if file path is correct!" << endl;
        exit(-1);
    }
  // std::cout << file.IsZombie() << endl;

	TTree* tree = new TTree;
  file.GetObject("T", tree);
  tree->GetEntry(1); **/

	string histNames[] = {"minimum_ch10", "minimum_ch11", "minimum_ch12", "minimum_ch13"}; //, "Incidence_time_ch10", "Incidence_time_ch11", "Incidence_time_ch12", "Incidence_time_ch13"};
	string histTitles[] = {"Amplitude ch10", "Amplitude ch11", "Amplitude ch12", "Amplitude ch13"}; // "Incidence Time Ch10", "Incidence Time Ch11", "Incidence Time Ch12", "Incidence Time Ch13"};
  int nrPlots = 4;
	TString histDraw[nrPlots];
	TString xlabel = "Amplitude [mV]";
	TString ylabel = "Number of Entries";
  Int_t xMin[] = {-120, -120, -120, -120};
  Int_t xMax[] = {5, 5, 5, 5};

	TH1F* histVec[nrPlots]; 

  TCanvas canvas("canvas", "Plastic Scintillators", 1557, 2000);
  TPaveLabel title(0.1, 0.96, 0.9, 0.99, Form("Incidence Time; Run %d, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot", runNumber, posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit));
  TPaveLabel xTitle(0, 0.01, 1, 0.03, xlabel);
  TPaveLabel yTitle(0.01, 0, 0.03, 1, ylabel);
  title.SetTextSize(.7);
  xTitle.SetTextSize(.7);
  yTitle.SetTextAngle(90);
  yTitle.SetTextSize(.017);
  xTitle.SetLineColor(0);
  yTitle.SetLineColor(0);
  title.SetBorderSize(0);
  xTitle.SetBorderSize(0);
  yTitle.SetBorderSize(0);
  title.SetFillColor(0);
  xTitle.SetFillColor(0);
  yTitle.SetFillColor(0);
  title.Draw();
  xTitle.Draw();
  yTitle.Draw();
  TPad graphPad("Graphs", "Graphs", 0.03, 0.03, 1, 0.96);
  graphPad.Draw();
  graphPad.cd();
  int row = 4;
  int col = 2;
  graphPad.Divide(col, row);
  //graphPad.SetLeftMargin(0.2);
  //graphPad.SetBottomMargin(0.15);

  // std::cout << "3" << endl;

  

  for (int nPad = 1; nPad < nrPlots+1; nPad++) {
		graphPad.cd(nPad);
		gPad->SetLeftMargin(.08); //.18
		gPad->SetBottomMargin(.052); //.15
		gPad->SetRightMargin(0.065);
		gPad->SetGrid();
  }

	for (int i = 0; i < nrPlots; i++) {
    graphPad.cd(i + 1);

		histDraw[i] = Form("%s>>%s", histNames[i].c_str(), histNames[i].c_str());

    histVec[i] = new TH1F(histNames[i].c_str(), histTitles[i].c_str(), (xMax[i] - xMin[i])*5, xMin[i], xMax[i]);

		gStyle->SetTitleSize(0.08, "t"); 
		
    histVec[i]->GetXaxis()->SetLabelSize(.06);
    histVec[i]->GetYaxis()->SetLabelSize(.06);
    histVec[i]->SetLineColorAlpha(kBlack, 0.7);
    histVec[i]->SetFillColorAlpha(kBlack, 0.5);
    histVec[i]->SetMarkerStyle(8);
    histVec[i]->SetMarkerSize(0.2);
    histVec[i]->SetMarkerColorAlpha(kBlack, 0.6);



    tree->Draw(histDraw[i], "", "HIST");

    // Customizing the legend:
    TLegend* histLeg = new TLegend(0.64, 0.58, 0.92, 0.87);
    // histLeg->SetFillColorAlpha(kWhite, 0); //translucent legend
    histLeg->SetBorderSize(1);
    histLeg->AddEntry((TObject*)0, Form("Entries = %d", (int) histVec[i]->GetEntries()));
    histLeg->AddEntry((TObject*)0, Form("Mean = %.2f", histVec[i]->GetMean()));
    histLeg->AddEntry((TObject*)0, Form("Mean = %.2f", histVec[i]->GetStdDev()));
    // histLeg->AddEntry((TObject*)0, Form("Mean = %1.2f #pm %1.2f", histMeanVec[i], histMeanErrVec[i]), "");
    gStyle->SetLegendTextSize(0.05);
		histLeg->Draw();



	}

	canvas.SaveAs(locationAmps.c_str());



}