// reads finished Rootfile and extracts Phi_ew histograms.
// fits histogram with no channel omitted to gauss and saves it
// makes pdf with histograms where each time one channel gets omitted

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
#include <TFitResult.h>
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

//andrea
#include "meanAngleFuncs.h"

using namespace std;

int main(int argc, char *argv[]) {
  // string rootfileStr;
  string rootfileStr = argv[1]; // takes in format e.g. ../RootAnalysis/finishedRootfiles/18_cosmics_vb58_50k_2808_PARTS.root
  std::cout << "1 " << rootfileStr << endl;

  // string rootfileStr = "/mnt/d/Programme/RootAnalysis/RootAnalysis/finishedRootfiles/17_cosmics_vb58_50k_2808.root";
  
  string runNameRaw;
  // cout << rootfileStr.size()-6 << endl;
  for (int i = rootfileStr.size()-6; i >= 0; i--) {
    if (rootfileStr.at(i) == '/' || rootfileStr.at(i) == '\\') {
        if (i == rootfileStr.size()-6) {
          runNameRaw = "EGNAHC_eman_tluafed"; //"default_name_CHANGE";
          // cout << 21 << endl;
        }
        // cout << 23 << endl;
        break;
    }
    // cout << rootfileStr[i] << endl;
    // cout << rootfileStr.at(i) << endl;
    runNameRaw += rootfileStr.at(i);
    // cout << runNameRaw << endl;
    // cout << 22 << endl;
  }
  // printf("Hello world");
  // std::cout << "3 " << runNameRaw << "\n" << endl;

  string runName = "";
  for (int i=0; i<runNameRaw.size(); i++) {
    // cout << runNameRaw.at(runNameRaw.size()-1-i) << endl;
    runName += runNameRaw.at(runNameRaw.size()-1-i);
  }
  // cout << "finished runName " << runName << endl;

  string runNumberString = "";
  for (int i=0; i<runName.size(); i++) {
    if (runName.at(i) == '_') break;
    runNumberString += runName.at(i);
  }
  int runNumber = stoi(runNumberString);

  // if not otherwise specified, files will be save to ../RootAnalysis/integralAnalysis/runName
  // HARDCODED PATH
  string saveFolder;
  if (argc < 3) { // || (string(argv[2])).empty() ) {
    saveFolder = "/mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/integralAnalysis/" + runName;
    // std::cout << " argc < 3 " << saveFolder <<"DUB" << endl;
  } else {
    // std::cout << "argc not <3  " << argc << endl;
    string saveFolder = argv[2]; // path specified when analysis is called (makeIntegralsPDF.sh)
    // std::cout << saveFolder << endl;
  }

  string allChannelsParam_loc = saveFolder + "/allChannelsParam.txt";
  FILE* allChannelsParam = fopen(allChannelsParam_loc.c_str(), "w");

  string omitChannelsParam_loc = saveFolder + "/omitChannelsParam.txt";
  FILE* omitChannelsParam = fopen(omitChannelsParam_loc.c_str(), "w");


  // open rootfile:
  TFile* file = new TFile(rootfileStr.c_str());
 
  if (file->IsZombie())
    {
        std::cout << "Problem with file " << rootfileStr << "; check if file path is correct!" << endl;
        exit(-1);
    }
  // std::cout << file.IsZombie() << endl;

  TTree* tree = new TTree;
  file->GetObject("T", tree);
  tree->GetEntry(1);

  //find out run parameters
  float rawAngleLimitLower, rawAngleLimitUpper;
  float rawTopPositionLimitLower, rawTopPositionLimitUpper;
  float rawBotPositionLimitLower, rawBotPositionLimitUpper;
  int centerChannel;

  tree->SetBranchAddress("angleIntTop", &rawAngleLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("angleIntBot", &rawAngleLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("posTopIntTop", &rawTopPositionLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("posTopIntBot", &rawTopPositionLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("posBotIntTop", &rawBotPositionLimitUpper);
  tree->GetEntry(1);
  tree->SetBranchAddress("posBotIntBot", &rawBotPositionLimitLower);
  tree->GetEntry(1);
  tree->SetBranchAddress("CenterChannelPhiEw", &centerChannel);
  tree->GetEntry(1);


  float angleUpperLimit; // angle interval right (from normal) deg
  float angleLowerLimit; // angle interval left (from normal) deg
  float posTopLeftLimit; // left position limit upper plastic scint, in cm from 0
  float posTopRightLimit; // right position limit upper plastic scint, in cm from 0
  float posBotLeftLimit; // left position limit lower plastic scint, in cm from 0
  float posBotRightLimit; // right position limit lower plastic scint, in cm from 0
  
  if (rawAngleLimitUpper == 999) {
    angleUpperLimit = 90.0;
    angleLowerLimit = -90.0;
  } 
  else {
    angleUpperLimit = 90.0 - TMath::ATan(0.66 * 2 / (0.2998 * rawAngleLimitUpper)) / TMath::Pi() * 180.0;
    angleLowerLimit = -90.0 - TMath::ATan(0.66 * 2 / (0.2998 * rawAngleLimitLower)) / TMath::Pi() * 180.0;
  }

  if (rawTopPositionLimitLower == 999) {
    posTopLeftLimit = -100;
    posTopRightLimit = 100;
  }
  else {
    posTopLeftLimit = 0.2998 * 0.5 * 0.5 * rawTopPositionLimitLower * 100.0;
    posTopRightLimit = 0.2998 * 0.5 * 0.5 * rawTopPositionLimitUpper * 100.0;
  }

  if (rawBotPositionLimitLower == 999) {
    posBotLeftLimit = -100;
    posBotRightLimit = 100;
  }
  else {
    posBotLeftLimit = 0.2998 * 0.5 * 0.5 * rawBotPositionLimitLower * 100.0;
    posBotRightLimit = 0.2998 * 0.5 * 0.5 * rawBotPositionLimitUpper * 100.0;
  }

  float positionInfo[] = {(float) runNumber, posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit};


  // declare stuff to access data
  int nCh = 8;
  TH1F* omitChannels[nCh];

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);
  gStyle->SetLineScalePS(1);

  TCanvas canvas("canvas", "Phi_ew for channels", 1557, 2000);
  TPaveLabel title(0.1, 0.96, 0.9, 0.99, Form("Phi_ew omitting different channels; centered opposite omitted channel, Run %s, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot", runNumberString.c_str(), posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit));
  // TPaveLabel title(0.1, 0.96, 0.9, 0.99, "Phi_ew omitting different channels; centered opposite omitted channel");
  TPaveLabel xTitle(0, 0.01, 1, 0.03, "Phi_ew (deg.)");
  TPaveLabel yTitle(0.01, 0, 0.03, 1, "Number of Entries");
  title.SetTextSize(.6);
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
  graphPad.Divide(2, 4);

  for (int nPad = 1; nPad < 9; nPad++) {
      graphPad.cd(nPad);
      gPad->SetLeftMargin(.065); //.18
      gPad->SetBottomMargin(.052); //.15
      gPad->SetRightMargin(0.065);
      gPad->SetGrid();
  }

  TString histTitle;
  TString histName;
  TString histDraw;

  Int_t xMax = 200;
  Int_t xMin = -200;


  for (int i=0; i<nCh; i++) {
    // omitChannels[i] = (TH1F*) file.Get(Form("Phi_ew_omit_ch%d", i));
    graphPad.cd(i + 1);

    //Naming the histogram and setting upt the draw command:
    histTitle.Form("Omitting Channel %d", i);
    histName.Form("Hist%d", i);
    histDraw.Form("Phi_ew_omit_ch%d>>Hist%d", i, i);

    omitChannels[i] = new TH1F(histName, histTitle, (xMax - xMin), xMin, xMax);

    gStyle->SetTitleSize(0.08, "t"); 
    //histVec[i]->SetTitleSize(0.25, "t");
    //histVec[i]->GetXaxis()->SetTitle("Number of photoelectrons");
    //histVec[i]->GetXaxis()->SetTitleSize(.07);
    omitChannels[i]->GetXaxis()->SetLabelSize(.06);
    //histVec[i]->GetYaxis()->SetTitle("Number of entries");
    //histVec[i]->GetYaxis()->SetTitleSize(.07);
    omitChannels[i]->GetYaxis()->SetLabelSize(.06);
    omitChannels[i]->SetLineColorAlpha(kBlack, 0.7);
    omitChannels[i]->SetFillColorAlpha(kBlack, 0.5);
    omitChannels[i]->SetMarkerStyle(8);
    omitChannels[i]->SetMarkerSize(0.2);
    omitChannels[i]->SetMarkerColorAlpha(kBlack, 0.6);

    tree->Draw(histDraw, "", "HIST");

    TLegend* histLeg = new TLegend(0.64, 0.68, 0.92, 0.85);
    histLeg->SetFillColorAlpha(kWhite, 0); //translucent legend
    histLeg->SetBorderSize(1);
    histLeg->AddEntry((TObject*)0, Form("Entries = %d", (int) omitChannels[i]->GetEntries()));
    // histLeg->AddEntry((TObject*)0, Form("Mean = %1.2f #pm %1.2f", histMeanVec[i], histMeanErrVec[i]), "");
    gStyle->SetLegendTextSize(0.05);

    histLeg->Draw();

  }

title.SetLabel(Form("#splitline{Phi_ew omitting different channels; centered opposite omitted channel,}{Run %s, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot}", runNumberString.c_str(), posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit));

string loc_name = saveFolder + "/phi_ew_omitting.pdf";
// std::cout << loc_name << endl;
canvas.SaveAs(loc_name.c_str());
std::cout << "- phi_ew_omitting.pdf" << endl;

canvas.Clear();



// save phi ew histogram with no omissions too

//Naming the histogram and setting upt the draw command:
TString histTitle2Shift;
TString histName2Shift;
TString histDraw2Shift;
histTitle2Shift.Form("#splitline{Phi_ew from all channels, centered around channel %d}{Shifted by + and -360 deg.}", centerChannel);
histName2Shift.Form("HistAllShift");
histDraw2Shift.Form("Phi_ew_shifted>>HistAllShift");
TString histTitle2;
TString histName2;
TString histDraw2;
histTitle2.Form("#phi_{ew} centered around channel %d", centerChannel);
histName2.Form("HistAll");
histDraw2.Form("Phi_ew_all_ch>>HistAll");
TString histTitle3;
TString histName3;
TString histDraw3;
histTitle3.Form("Standard deviation Phi_ew");
histName3.Form("HistStd");
histDraw3.Form("Std_Phi_ew_all>>HistStd");

TString histNamesPhi[] = {histName2, histName3, histName2Shift};
TString histTitlesPhi[] = {histTitle2, histTitle3, histTitle2Shift};
TString histDrawsPhi[] = {histDraw2, histDraw3, histDraw2Shift};
// cout << 2 << endl;
TH1F* phiShifted = new TH1F(histName2Shift, histTitle2Shift, 3*(xMax - xMin)/5, 3*xMin, 3*xMax);
TH1F* allChannels = new TH1F(histName2, histTitle2, 180, xMin, xMax);
TH1F* stdPhiew = new TH1F(histName3, histTitle3, (1200)/10, 0, 250);

TH1F* phiEwStuff[3] = {allChannels, stdPhiew, phiShifted};
// cout << 3 << endl;

// gStyle->SetTitleFont(28,'t');

// cout << 1 << endl;
TCanvas canvas2("canvas", "Phi_ew from all channels", 1500, 1000);
TPaveLabel title2(0.1, 0.96, 0.9, 0.99, Form("Phi_ew centered opposite omitted channel, Run %s, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot", runNumberString.c_str(), posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit));
TPaveLabel xTitle2(0, 0.01, 1, 0.03, "Phi_ew (deg.)");
TPaveLabel yTitle2(0.03, 0, 0.03, 1, "#splitline{Number of Entries}{ }");

xTitle2.SetTextSize(1.2);
yTitle2.SetTextAngle(90);
yTitle2.SetTextSize(.03);
// xTitle2.SetLineColor(0);
// yTitle2.SetLineColor(0);

xTitle2.SetBorderSize(0);
yTitle2.SetBorderSize(0);

xTitle2.SetFillColor(0);
yTitle2.SetFillColor(0);

// xTitle2.Draw();
// yTitle2.Draw();
title2.Draw();
TPad graphPad2("Graphs", "Graphs", 0.01, 0.03, 1, 0.96);
graphPad2.Draw();
graphPad2.Divide(2, 2);
graphPad2.cd();
TString xAxesPhiAll[] = {"#phi_{ew} (deg.)", "#phi_{ew} (deg.)", "#phi_{ew} (deg.)"};


// fit stuff
// TF1* g1 = new TF1("m1", "gaus", -540, -180);
// TF1* g2 = new TF1("m2", "gaus", -180, 180);
// TF1* g3 = new TF1("m3", "gaus", 180, 540);

TF1 *func = new TF1("tripleGauss", fitf, -540, 540, 4);


// TF1* total = new TF1("mstotal", "gaus(0)+gaus(3)+gaus(6)", -540,540);
// Double_t parameters[9];




//histVec[i]->SetTitleSize(0.25, "t");
//histVec[i]->GetXaxis()->SetTitle("Number of photoelectrons");
//histVec[i]->GetXaxis()->SetTitleSize(.07);
// allChannels->GetXaxis()->SetLabelSize(.06);
//histVec[i]->GetYaxis()->SetTitle("Number of entries");
//histVec[i]->GetYaxis()->SetTitleSize(.07);
// allChannels->GetYaxis()->SetLabelSize(.06);
int is[] = {1, 2, 3};
for (int i=0; i<3; i++) {
  graphPad2.cd(is[i]);

  gPad->SetLeftMargin(.12); //.18
  gPad->SetBottomMargin(.152); //.15
  gPad->SetRightMargin(0.065);
  gPad->SetTopMargin(.152);
  gPad->SetGrid();

  // gStyle->SetTitleSize(0.05, "t"); 
  // gStyle->SetTitleFontSize(13);
  // gStyle->SetTitleFont(70, "t");
  gStyle->SetTitleOffset(3, "t");
  gStyle->SetTitleY(1.01);
  // gStyle->SetTitleFont(30,'t');

  phiEwStuff[i]->SetLineColorAlpha(kBlack, 0.7);
  phiEwStuff[i]->SetFillColorAlpha(kBlack, 0.5);
  phiEwStuff[i]->SetMarkerStyle(8);
  phiEwStuff[i]->SetMarkerSize(0.2);
  phiEwStuff[i]->SetMarkerColorAlpha(kBlack, 0.6);
  phiEwStuff[i]->GetXaxis()->SetTitle(xAxesPhiAll[i]);
  phiEwStuff[i]->GetXaxis()->SetTitleOffset(1);
  phiEwStuff[i]->GetXaxis()->SetTitleSize(0.04);
  phiEwStuff[i]->GetYaxis()->SetTitle("Counts");
  phiEwStuff[i]->GetYaxis()->SetTitleOffset(1.3);
  phiEwStuff[i]->GetYaxis()->SetTitleSize(0.04);


  

  // fit triple gauss
  if (i==2) {
    string fitLog_loc = saveFolder + "/.." + "/fitlog.txt";
    FILE* fitLog = fopen(fitLog_loc.c_str(), "a");
    tree->Draw(histDrawsPhi[i]);
    func->SetParameters(phiEwStuff[i]->GetMaximum(), phiEwStuff[0]->GetMean(), phiEwStuff[0]->GetRMS(), phiEwStuff[0]->GetMinimum());
    // cout << "here" << endl;
    // cout << phiEwStuff[0]->GetMaximum() << endl;
    // cout << phiEwStuff[0]->GetMean() << endl;
    // cout << phiEwStuff[0]->GetRMS() << endl;
    // cout << phiEwStuff[0]->GetMinimum() << endl;
    func->SetParNames("Maximum", "Mean", "Sigma", "Base");
    func->SetParLimits(0, 0, 10000);   //Making sure the gauss function is positive
    func->SetParLimits(1, -180, 180); //Making sure the central fit peak is in the histogram range
    func->SetParLimits(2, 0, 10000);   //Making sure the standard deviation is positive (could also just use absolute value); upper limits are arbitrary but should cover everything
    func->SetLineWidth(5);
    phiEwStuff[i]->Fit("tripleGauss");
    TFitResultPtr resultPtr = phiEwStuff[i]->Fit(func,"S");
    TMatrixDSym cov = resultPtr->GetCovarianceMatrix();
    cout << "covariance" << endl;
    cov.Print();
    // cov.Write((saveFolder+"/covarianceMatrix.txt").c_str());
    // resultPtr->Print("V") >> fitLog;
    // resultPtr->Write();

    


    tree->Draw(histDrawsPhi[i]); //, "", "HIST");
    func->Draw("same");

    // TF1* finishedFit = phiEwStuff[i]->GetFunction("fit");
    Double_t chi2 = func->GetChisquare();
    Double_t ndf = func->GetNDF();
    Double_t fitParams[4];
    Double_t fitParamErrs[4];
    for (int k=0; k<4; k++) {
      fitParams[k] = func->GetParameter(k);
      // cout << func->GetParameter(k) << endl;
      fitParamErrs[k] = func->GetParError(k);
      // printArray(fitParams, 4);
      // printArray(fitParamErrs, 4);
    }

    Double_t covMatrixSave[4][4];
    
    for (int m=0; m<4; m++) {
      for (int n=0; n<4; n++) {
        covMatrixSave[m][n] = cov(m,n);
      }
    }


    // fprintf(fitLog, "#Max\tPeak\tSigma\tOffset\n");
    fprintf(fitLog, "\n##Run %s, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot", runNumberString.c_str(), posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit);
    fprintf(fitLog, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", fitParams[0], fitParamErrs[0], fitParams[1], fitParamErrs[1], fitParams[2], fitParamErrs[2], fitParams[3], fitParamErrs[3], chi2/ndf);
    fprintf(fitLog, "#Covariance Matrix\n");
    for (int m=0; m<4; m++) {
      fprintf(fitLog, "%f\t%f\t%f\t%f\n", cov(m,0), cov(m,1), cov(m,2), cov(m,3));
    }

    //Ratio of the y values of the central fit maximum and the minimum between the central and the left maximum:
    float fitMax = func->GetMaximum(func->GetParameter(1) - 1, func->GetParameter(1) + 1); //Searches for maximum in the region around the central peak
    float fitMin = func->Eval(func->GetParameter(1) - 180.0); //Two maxima are 360 deg apart and the minimum is in the middle
    float peakValleyRatio = fitMax / fitMin;

    float sigma = fitParams[2];
    float scaling = fitParams[0];
    float offset = fitParams[3];
    float sigma2 = sigma*sigma;
    float sigmaError = fitParamErrs[2];
    float scalingError = fitParamErrs[0];
    float offsetError = fitParamErrs[3];

    // yanked this from joscha, will check over later
    //Derivatives of f_max / f_min with respect to parameters (calculated by Mathematica, which is great for symbolic calculations but nothing else:
    float scalingDerivative = (1 + 2 * TMath::Exp(-32400.0 / sigma2)) / ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset)
        - ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset)) / TMath::Power(((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset), 2);
    float sigmaDerivative = 129600.0 * TMath::Exp(-32400.0 / sigma2) * scaling / (TMath::Power(sigma, 3) * ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset))
        - scaling * ((145800 * TMath::Exp(-72900.0 / sigma2) + 32400.0 * TMath::Exp(-8100.0 / sigma2)) / TMath::Power(sigma, 3)) * ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset)
        / TMath::Power((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset, 2);
    float offsetDerivative = 1 / ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset)
        - ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset) / TMath::Power((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100 / sigma2)) * scaling + offset, 2);

    float peakValleyRatioErrors = TMath::Sqrt(TMath::Power(scalingDerivative * scalingError, 2) + TMath::Power(sigmaDerivative * sigmaError, 2) + TMath::Power(offsetDerivative * offsetError, 2)); // ADD COVARIANCE STUFF
    
    fprintf(fitLog, "#Peak Valley Ratio\n");
    fprintf(fitLog, "%f\t%f\n", peakValleyRatio, peakValleyRatioErrors);

    TLegend* histLeg2 = new TLegend(0.83, 0.5, 1.0, 0.83);
    histLeg2->SetFillColorAlpha(kWhite, 0.91); //translucent legend
    histLeg2->SetBorderSize(1);

    histLeg2->AddEntry((TObject*)0, Form("Entries = %d", (int) phiEwStuff[i]->GetEntries()), "");
    histLeg2->AddEntry((TObject*)0, Form("#bar{#Phi_{ew}} = %.2f", fitParams[1]), "");
    histLeg2->AddEntry((TObject*)0, Form("#bar{#Phi_{ew}} err. = %.2f", fitParamErrs[1]), "");
    histLeg2->AddEntry((TObject*)0, Form("#sigma = %.2f", fitParams[2]), "");
    histLeg2->AddEntry((TObject*)0, Form("#sigma err. = %.2f", fitParamErrs[2]), "");
    histLeg2->AddEntry((TObject*)0, Form("p_{3} = %.2f", fitParams[3]), "");
    histLeg2->AddEntry((TObject*)0, Form("p_{3} err. = %.2f", fitParamErrs[3]), "");
    histLeg2->AddEntry((TObject*)0, Form("#Chi^{2}/dof = %.2f", chi2/ndf), "");
    // histLeg->AddEntry((TObject*)0, Form("Mean = %1.2f #pm %1.2f", histMeanVec[i], histMeanErrVec[i]), "");
    gStyle->SetLegendTextSize(0.03);
    histLeg2->Draw();
  }
  else {
    // gStyle->SetPalette(kOcean);
    tree->Draw(histDrawsPhi[i]); //, "", "HIST");
    TLegend* histLeg2 = new TLegend(0.9, 0.8, 0.68, 0.65);
    histLeg2->SetFillColorAlpha(kWhite, 0.7); //translucent legend
    histLeg2->SetBorderSize(1);
    histLeg2->AddEntry((TObject*)0, Form("Entries = %d", (int) phiEwStuff[i]->GetEntries()));
    histLeg2->AddEntry((TObject*)0, Form("Mean = %1.2f #pm %1.2f", phiEwStuff[i]->GetMean(), phiEwStuff[i]->GetMeanError()), "");
    gStyle->SetLegendTextSize(0.03);
    histLeg2->Draw();
  }
}

title2.SetLabel(Form("Phi_ew centered opposite omitted channel, Run %s, [%.1f, %.1f]cm top., [%.1f, %.1f]cm bot", runNumberString.c_str(), posTopLeftLimit, posTopRightLimit, posBotLeftLimit, posBotRightLimit));





string loc_name2 = saveFolder + "/phi_ew_all.pdf";
canvas2.SaveAs(loc_name2.c_str());
std::cout << "- phi_ew_all.pdf" << endl;

return 0;

}