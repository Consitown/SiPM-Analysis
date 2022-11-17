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

//Including standard C/C++-Libraries:
#include <sys/stat.h>
#include <sys/types.h>
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

int main(int argc, char* argv[])
{
    const int maxFileNumber = 100;        //Arbitrary; change if needed
    string outPath = "/mnt/d/RootAnalysis/RootAnalysis/integralAnalysis/joscha_test";
    TTree* tree[maxFileNumber];
    const Int_t nCh = 8;                    //eight channels for WOM D (sum is not needed for this analysis at the moment)
    
    //Style Settings:
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(16);
    gStyle->SetLineScalePS(1);

    //Creating variables to read data from the TTree stored in the .file into our own TTree:
    int runNumber[maxFileNumber];
    int runPosition[maxFileNumber];
    float runEnergy[maxFileNumber];
    float photonCountPerEvent[maxFileNumber];
    int nEntries[maxFileNumber];

    TLine* meanLineVec[maxFileNumber][nCh];                //array of lines for mean values
    TLine* leftEdgeLineVec[maxFileNumber][nCh];            //array of lines for left edges during mean calculation
    Float_t histMaxY[maxFileNumber][nCh];                  //y values for maximum amplitude

    //Opening a file to write photon counts into:
    FILE* photonCountFile = fopen(Form("%sphotonCounts.txt", outPath.c_str()), "w");

    //Creating an array to calculate reduced mean values for the three energies (dividing light yields for pos 9 and 11 by those for pos 10):
    float normValues14[nCh][6];     //Three values and three errors for each channel
    float normValues26[nCh][6];
    float normValues52[nCh][6];

    //Mapping switches to their respective angles:
    Int_t angles[nCh] = { 0,45,90,135,180,225,270,315 };
    Int_t channelOrder[nCh] = { 6,5,4,0,1,2,3,7 }; //0 deg --> channel 6, 45 deg --> channel 5, 90 deg --> channel 4, ... , 315 deg --> channel 7
    Int_t channelLabels[nCh] = { 135,180,225,270,90,45,0,315 };

    //Loop over files:
    for (int arg = 1; arg < argc; arg++)
    {
        tree[arg] = new TTree;
        string filePath = "/mnt/d/Programme/RootAnalysis/RootAnalysis/finishedRootfiles/18_cosmics_vb58_50k_2808_PARTS.root";
        printf("Analyzing file %s\n", filePath.c_str());

        //Open .root file:
        TFile file(filePath.c_str());
        if (file.IsZombie())
        {
            cout << "Problem with file " << filePath << "; check if file path is correct!" << endl;
            exit(-1);
        }

        //Reading Data from file into our own TTree:
        file.GetObject("T", tree[arg]);
        tree[arg]->SetBranchAddress("runNumber", &runNumber[arg]);
        tree[arg]->SetBranchAddress("runPosition", &runPosition[arg]);
        tree[arg]->SetBranchAddress("runEnergy", &runEnergy[arg]);
        tree[arg]->GetEntry(1);

        //Drawing range for histograms:
        Int_t xMin, xMax;
        xMin = -10;
        xMax = 500;

        //Retrieve light yield data for all channels. IMPORTANT: This ENTIRE script is currently only set up for the analysis of WOM D!
        TH1F* histVec[nCh];                     //array of histograms
        Float_t histMeanVec[nCh];               //array of mean values
        Float_t histMeanErrVec[nCh];            //array of mean value errors
        Float_t histStdDevVec[nCh];             //array of standard deviations
        Float_t histStdDevErrVec[nCh];          //array of standard deviation errors

        TString histName, histTitle, histDraw;  //name and title for histogram; as well as the command to draw it which will be formatted later
        Int_t histEntries;                      //number of entries in each histogram
        int histLeftEdge[nCh];                //array of minimum values for mean calculation

        //Int_t photonCut[nCh];                    //specifies at which value the photon cut for MeanAngles.cpp will take place

        //Setting up the canvas:
        TCanvas canvas("canvas", "Light yield for different Channels", 1557, 2000);
        TPaveLabel title(0.1, 0.96, 0.9, 0.99, Form("Light yield for different Channels, Run %d, Position %d, %1.1f GeV", runNumber[arg], runPosition[arg], runEnergy[arg] / 10.0));
        TPaveLabel xTitle(0, 0.01, 1, 0.03, "Number of Photoelectrons N_{pe}");
        TPaveLabel yTitle(0.01, 0, 0.03, 1, "Number of Entries");
        title.SetTextSize(.7);
        xTitle.SetTextSize(.7);
        yTitle.SetTextAngle(90);
        yTitle.SetTextSize(.017);
        title.SetLineColor(0);
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
        //graphPad.SetLeftMargin(0.2);
        //graphPad.SetBottomMargin(0.15);
        for (int nPad = 1; nPad < 9; nPad++)
        {
            graphPad.cd(nPad);
            gPad->SetLeftMargin(.065); //.18
            gPad->SetBottomMargin(.052); //.15
            gPad->SetRightMargin(0.065);
            gPad->SetGrid();
        }

        photonCountPerEvent[arg] = 0;
        //Loop over all channels of WOM D (in Jans read.C script at line 670, it is declared that of the 32 overall channels, WOM D is read out in the first eight (0-7).
        for (int i = 0; i < nCh; i++)
        {
            graphPad.cd(i + 1);

            //Naming the histogram and setting upt the draw command:
            histTitle.Form("Channel %d (%d degrees)", i, channelLabels[i]);
            histName.Form("Hist%d", i);
            histDraw.Form("Integral[%d]>>Hist%d", i, i);

            //Drawing the Histogram:
            histVec[i] = new TH1F(histName, histTitle, (xMax - xMin), xMin, xMax);
            gStyle->SetTitleSize(0.08, "t"); //histVec[i]->SetTitleSize(0.25, "t");
            //histVec[i]->GetXaxis()->SetTitle("Number of photoelectrons");
            //histVec[i]->GetXaxis()->SetTitleSize(.07);
            histVec[i]->GetXaxis()->SetLabelSize(.06);
            //histVec[i]->GetYaxis()->SetTitle("Number of entries");
            //histVec[i]->GetYaxis()->SetTitleSize(.07);
            histVec[i]->GetYaxis()->SetLabelSize(.06);
            histVec[i]->SetLineColorAlpha(kBlack, 0.7);
            histVec[i]->SetFillColorAlpha(kBlack, 0.5);
            histVec[i]->SetMarkerStyle(8);
            histVec[i]->SetMarkerSize(0.2);
            histVec[i]->SetMarkerColorAlpha(kBlack, 0.6);
            tree[arg]->Draw(histDraw, "", "HIST");

            //Getting histogram maximum (how high to draw the lines):
            histMaxY[arg][i] = histVec[i]->GetBinContent(histVec[i]->GetMaximumBin());

            /*
            //Finding the left edge for mean calculation:
            histLeftEdge[i] = 10 + TMath::Floor(0.1 * histVec[i]->GetMean());
            //Start out at a bin right of the leftmost bin (ignore 0 --> Underflow, 1-9 --> below 0), depending on where the mean is (not starting out at the first bin makes the algorithm more stable;
            //e.g. if the first bin is randomly very small, the condition below will instantly be fulfilled, which we do not want
            float leftSum = 0;
            float rightSum = 0;
            float leftMean = 0;
            float rightMean = 0;
            Double_t binContent;
            int counter;

            while (rightMean <= 1.3 * leftMean)
            {
                if (counter >= 1000)    //Sometimes the condition is not met at all -- in this case, terminate after 1000 (arbitrary number) tries and continue with next histogram
                {
                    histLeftEdge[i] == 0;
                    break;
                }
                counter += 1;
                leftSum = 0.0;
                rightSum = 0.0;
                for (int bin = 10; bin <= histLeftEdge[i]; bin++)
                {
                    binContent = histVec[i]->GetBinContent(bin);
                    leftSum += binContent;    //Sum of all bins left of histLeftEdge[i] (including histLeftEdge[i])
                }
                leftMean = (float)leftSum / (histLeftEdge[i] - 9);
                for (int bin = 1; bin <= 10; bin++)
                {
                    rightSum += histVec[i]->GetBinContent(histLeftEdge[i] + bin);
                }
                rightMean = (float)rightSum / 10.0;     //Mean of five bins right of histLeftEdge[i] (excluding histLeftEdge[i])
                histLeftEdge[i] += 1;                   //Move histLeftEdge[i] to the right until the condition is met
            }
            photonCut[i] = 0;
            if (histLeftEdge[i] - 11 > 0 && histLeftEdge[i] > 10 + TMath::Ceil(0.1 * histVec[i]->GetMean()))
            {
                photonCut[i] = histLeftEdge[i] - 11;
            }

            //Declaring the leftmost bin for the calculation of the histogram mean and the total light yield:
            histVec[i]->GetXaxis()->SetRange(histLeftEdge[i], xMax); //sets axis range to start at the specified lower bin, since only bins inside axis range are included in GetMean-calculation
            */
            //Obtaining results:
            histMeanVec[i] = histVec[i]->GetMean();
            histMeanErrVec[i] = histVec[i]->GetMeanError();
            histStdDevVec[i] = histVec[i]->GetStdDev();
            histStdDevErrVec[i] = histVec[i]->GetStdDevError();

            //Getting histogram sum (total number of photons, for light yield vs f_{max}/f_{min} plot):
            int nCells = histVec[i]->GetNcells();
            nEntries[arg] = histVec[i]->GetEntries();
            for (int j = 1; j <= nCells; j++)
            {
                if (histVec[i]->GetBinContent(j) > 0)
                {
                    photonCountPerEvent[arg] += histVec[i]->GetBinContent(j) * (float)j;
                }
            }

            //Drawing mean lines:
            meanLineVec[arg][i] = new TLine(histMeanVec[i], 0, histMeanVec[i], histMaxY[arg][i] * 1.01);
            meanLineVec[arg][i]->SetLineColor(4);
            meanLineVec[arg][i]->SetLineStyle(2);
            meanLineVec[arg][i]->SetLineWidth(4);
            meanLineVec[arg][i]->Draw("same");

            //Drawing lines for leftmost bin included in mean calculation:
            /*if (photonCut[i] > 0)
            {
                leftEdgeLineVec[arg][i] = new TLine(photonCut[i], 0, photonCut[i], histMaxY[arg][i] * 1.01);
                leftEdgeLineVec[arg][i]->SetLineColor(46);
                leftEdgeLineVec[arg][i]->SetLineStyle(5);
                leftEdgeLineVec[arg][i]->SetLineWidth(4);
                leftEdgeLineVec[arg][i]->Draw("same");
            }
            */

            //Resetting Drawing Range:
            histVec[i]->GetXaxis()->SetCanExtend(kTRUE);
            histVec[i]->ExtendAxis(histMeanVec[i] * 5, histVec[i]->GetXaxis());
            histVec[i]->GetXaxis()->SetRangeUser(histMeanVec[i] * (-0.2), TMath::Min(histMeanVec[i] * 4.99, (double)histVec[i]->GetNbinsX()));

            //Customizing the legend:
            TLegend* histLeg = new TLegend(0.45, 0.56, 0.92, 0.92);
            histLeg->SetFillColorAlpha(kWhite, 0); //translucent legend
            histLeg->SetBorderSize(0);
            histLeg->AddEntry(meanLineVec[arg][i], "#bar{N_{pe}} = ", "l");
            histLeg->AddEntry((TObject*)0, Form("%1.2f #pm %1.2f", histMeanVec[i], histMeanErrVec[i]), "");
            /*if (photonCut[i] > 0)
            {
                histLeg->AddEntry(leftEdgeLineVec[arg][i], Form("Cutoff Value: %d", photonCut[i]), "l");
            }
            */
            //histLeg->AddEntry((TObject*)0, Form("#sigma = %1.2f #pm %1.2f", histStdDevVec[i], histStdDevErrVec[i]), "");
            histLeg->Draw();
        }   //END OF LOOP OVER CHANNELS
        title.SetLabel(Form("Light yield for different Channels, Run %d, Pos. %d, %1.1f GeV (%d Entries)", runNumber[arg], runPosition[arg], runEnergy[arg] / 10.0, nEntries[arg]));

        //Average photon count per Event:
        photonCountPerEvent[arg] /= nEntries[arg];

        //Exporting canvas to pdf:
        string pdfName = Form("run%d_pos%d_energy%2.0f.pdf", runNumber[arg], runPosition[arg], runEnergy[arg]);
        canvas.SaveAs(Form("%s%s", outPath.c_str(), pdfName.c_str()));

        cout << "Run " << runNumber[arg] << " complete!" << endl;

        //Creating/Opening a file for mean values and errors, and one for Cutoff values:
        FILE* MeanList = fopen(Form("%sMeanList_run%d.txt", outPath.c_str(), runNumber[arg]), "w");
        //FILE* CutoffValues = fopen(Form("%sCutoffValues_run%d.txt", outPath.c_str(), runNumber[arg]), "w");

        //Computing the sum of all means (so that all light yields can be divided by the total light yield for that run for better comparability):
        Float_t histMeanSum = 0;
        Float_t histMeanErrSum = 0;
        Float_t propagatedError[nCh];
        for (int i = 0; i < nCh; i++)
        {
            histMeanSum += histMeanVec[i];
            histMeanErrSum += TMath::Power(histMeanErrVec[i], 2);
        }

        //Computing the propagated error for every value divided by the sum of values:
        for (int i = 0; i < nCh; i++)
        {
            Float_t x_i = histMeanVec[i];
            Float_t dx_i = histMeanErrVec[i];
            Float_t x_sum = histMeanSum;
            Float_t dx_sum = histMeanErrSum - dx_i * dx_i;
            propagatedError[i] = TMath::Sqrt(TMath::Power(((1 / x_sum) - (x_i / (x_sum * x_sum))) * dx_i, 2) - TMath::Power((x_i / (x_sum * x_sum)), 2) * dx_sum);
        }

        //Creating seperate .txt files for each run. 
        //In each MeanList file, the angle corresponding to the channel is written in the first column, mean values in the second and mean errors into the third, separated by tabulators.
        
        for (int i = 0; i < nCh; i++)
        {
            fprintf(MeanList, "%d\t", angles[i]);
            fprintf(MeanList, "%f\t", histMeanVec[channelOrder[i]] / histMeanSum);
            fprintf(MeanList, "%f\n", propagatedError[channelOrder[i]]);
        }
        fclose(MeanList);
        
        if (runNumber[arg] == 30)   //Pos 9, 2.6 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues26[i][0] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues26[i][1] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 31)   //Pos 9, 1.4 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues14[i][0] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues14[i][1] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 32)   //Pos 9, 5.2 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues52[i][0] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues52[i][1] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 33)   //Pos 10, 5.2 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues52[i][2] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues52[i][3] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 34)   //Pos 10, 2.6 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues26[i][2] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues26[i][3] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 35)   //Pos 10, 1.4 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues14[i][2] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues14[i][3] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 36)   //Pos 11, 1.4 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues14[i][4] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues14[i][5] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 37)   //Pos 11, 2.6 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues26[i][4] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues26[i][5] = propagatedError[channelOrder[i]];
            }
        }
        if (runNumber[arg] == 38)   //Pos 11, 5.2 GeV
        {
            for (int i = 0; i < nCh; i++)
            {
                normValues52[i][4] = histMeanVec[channelOrder[i]] / histMeanSum;
                normValues52[i][5] = propagatedError[channelOrder[i]];
            }
        }
        
        fprintf(photonCountFile, "%f\n", photonCountPerEvent[arg]);
    } //END OF LOOP OVER FILES

    FILE* ReducedValues14 = fopen(Form("%sreducedValues14.txt", outPath.c_str()), "w");
    FILE* ReducedValues26 = fopen(Form("%sreducedValues26.txt", outPath.c_str()), "w");
    FILE* ReducedValues52 = fopen(Form("%sreducedValues52.txt", outPath.c_str()), "w");
    for (int i = 0; i < nCh; i++)
    {
        //Calculate propagated errors for reduced values:
        Float_t x1_14 = normValues14[i][0];     //Values for pos 9
        Float_t x2_14 = normValues14[i][2];     //Values for pos 10
        Float_t x3_14 = normValues14[i][4];     //Values for pos 11
        Float_t dx1_14 = normValues14[i][1];    //Errors for pos 9
        Float_t dx2_14 = normValues14[i][3];    //Errors for pos 10
        Float_t dx3_14 = normValues14[i][5];    //Errors for pos 11
        Float_t propagatedErrorNorm9_14 = TMath::Sqrt(TMath::Power(1 / x2_14 * dx1_14, 2) + TMath::Power(-x1_14 / TMath::Power(x2_14, 2) * dx2_14, 2));
        Float_t propagatedErrorNorm11_14 = TMath::Sqrt(TMath::Power(1 / x2_14 * dx3_14, 2) + TMath::Power(-x3_14 / TMath::Power(x2_14, 2) * dx2_14, 2));
        Float_t x1_26 = normValues26[i][0];     //Values for pos 9
        Float_t x2_26 = normValues26[i][2];     //Values for pos 10
        Float_t x3_26 = normValues26[i][4];     //Values for pos 11
        Float_t dx1_26 = normValues26[i][1];    //Errors for pos 9
        Float_t dx2_26 = normValues26[i][3];    //Errors for pos 10
        Float_t dx3_26 = normValues26[i][5];    //Errors for pos 11
        Float_t propagatedErrorNorm9_26 = TMath::Sqrt(TMath::Power(1 / x2_26 * dx1_26, 2) + TMath::Power(-x1_26 / TMath::Power(x2_26, 2) * dx2_26, 2));
        Float_t propagatedErrorNorm11_26 = TMath::Sqrt(TMath::Power(1 / x2_26 * dx3_26, 2) + TMath::Power(-x3_26 / TMath::Power(x2_26, 2) * dx2_26, 2));
        Float_t x1_52 = normValues52[i][0];     //Values for pos 9
        Float_t x2_52 = normValues52[i][2];     //Values for pos 10
        Float_t x3_52 = normValues52[i][4];     //Values for pos 11
        Float_t dx1_52 = normValues52[i][1];    //Errors for pos 9
        Float_t dx2_52 = normValues52[i][3];    //Errors for pos 10
        Float_t dx3_52 = normValues52[i][5];    //Errors for pos 11
        Float_t propagatedErrorNorm9_52 = TMath::Sqrt(TMath::Power(1 / x2_52 * dx1_52, 2) + TMath::Power(-x1_52 / TMath::Power(x2_52, 2) * dx2_52, 2));
        Float_t propagatedErrorNorm11_52 = TMath::Sqrt(TMath::Power(1 / x2_52 * dx3_52, 2) + TMath::Power(-x3_52 / TMath::Power(x2_52, 2) * dx2_52, 2));

        //Print angle in first column:
        fprintf(ReducedValues14, "%d\t", angles[i]);
        fprintf(ReducedValues26, "%d\t", angles[i]);
        fprintf(ReducedValues52, "%d\t", angles[i]);
        //Second and third column are reduced values and errors for pos 9:
        fprintf(ReducedValues14, "%f\t%f\t", normValues14[i][0] / normValues14[i][2], propagatedErrorNorm9_14);
        fprintf(ReducedValues26, "%f\t%f\t", normValues26[i][0] / normValues26[i][2], propagatedErrorNorm9_26);
        fprintf(ReducedValues52, "%f\t%f\t", normValues52[i][0] / normValues52[i][2], propagatedErrorNorm9_52);
        //Fourth and fifth column are reduced values and errors for pos 10 (trivial):
        fprintf(ReducedValues14, "%f\t%f\t", 1.0, 0.0);
        fprintf(ReducedValues26, "%f\t%f\t", 1.0, 0.0);
        fprintf(ReducedValues52, "%f\t%f\t", 1.0, 0.0);
        //Sixth and seventh column are reduced values and errors for pos 11:
        fprintf(ReducedValues14, "%f\t%f\n", normValues14[i][4] / normValues14[i][2], propagatedErrorNorm11_14);
        fprintf(ReducedValues26, "%f\t%f\n", normValues26[i][4] / normValues26[i][2], propagatedErrorNorm11_26);
        fprintf(ReducedValues52, "%f\t%f\n", normValues52[i][4] / normValues52[i][2], propagatedErrorNorm11_52);
    }//Now we have three files (for 3 energies) with 7 columns and 8 rows each: 1st column is angles, 2-3 are pos 9, ..., 6-7 are pos 11 reduced values and errors.
     //This file is now accessed by AngularDistribution.cpp.
    fclose(ReducedValues14);
    fclose(ReducedValues26);
    fclose(ReducedValues52);

    fclose(photonCountFile);
}