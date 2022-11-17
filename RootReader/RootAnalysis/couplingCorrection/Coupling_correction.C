// reads finished Rootfile from (input) folder. Prints generalized x_i and y_i
// lists all these in .txt file generalized_coordinates.txt

//Including root functionalities:
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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
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

int main(int argc, char *argv[]) {
    //defining the directories
    string runName = "1333_1307to0408_day";

    string file_path = "/home/alex/Dokumente/Studium/RootReader/runs/Sicherung_0/";

    string in_filename = file_path + runName + ".root";

    string saveFolder = file_path + "couplingCorrection/" + runName;

    string data_filename = file_path +  "rootfiles/"; //file_path + "integralAnalysis/" + runName;

    //opening the root file
    cout << "Opening root file at " << in_filename << endl;
    TFile file(in_filename.c_str());

    if (file.IsZombie())
    {
        cout << "Problem with file " << in_filename << "; check if file path is correct!" << endl;
        exit(-1);
    }

    TTree* tree = new TTree;
    file.GetObject("T", tree);
    tree->GetEntry(1);

    //initializing the variables
    int n_ch = 8; //number of channels
    int channel_list[8];
    for (int i=0; i<n_ch; i++) channel_list[i] = i;
    int n_events; //number of events
    n_events = tree->GetEntries();
    
    double x_bar; //final x_bar (?)
    double y_bar; //final y_bar (?)

    float x_bar_channel[n_ch] = {0}; //x_bar for each channel
    float y_bar_channel[n_ch] = {0}; //y_bar for each channel

    float ly_event_channel[n_ch][n_events]; //lightyield of each event, so the integral value of one event that is saved in the root file
    float ly_event_res[n_events] = {0}; //lightyield from each event, summed over all channels

    double phi_list[n_ch] = {0, 315, 270, 225, 45,  90, 135, 180}; //angles corresponding to the channels in degree
    //for (int i=0; i<n_ch; i++) phi_list[i] = i * 45;

    //reading the individual integral values from the root files into the ly_event list
    for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < n_events; event_counter++)
        {
            TString integral_name;
            integral_name.Form("integral_hist_%d", channel_counter);
            tree->SetBranchAddress(integral_name, &ly_event_channel[channel_counter][event_counter]);
            tree->GetEntry(event_counter);

            ly_event_res[event_counter] += ly_event_channel[channel_counter][event_counter];
        }
    }

    cout << "Finished opening and reading." << endl;

    //calculate the x and y bar

    for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < n_events; event_counter++)
        {
            x_bar_channel[channel_counter] += TMath::Cos(phi_list[channel_counter] / 180. * TMath::Pi()) * ly_event_channel[channel_counter][event_counter] / ly_event_res[event_counter] / n_events;
            y_bar_channel[channel_counter] += TMath::Sin(phi_list[channel_counter] / 180. * TMath::Pi()) * ly_event_channel[channel_counter][event_counter] / ly_event_res[event_counter] / n_events;
            //cout << ly_event_channel[channel_counter][event_counter] << "\t" << ly_event_res[event_counter] << "\t" << n_events << endl;
        }
    }

    

    //calculate x bar sum and y bar sum and mean 
    double x_bar_sum, y_bar_sum;
    for (int i = 0; i < n_ch; i++)
    {
        x_bar_sum += x_bar_channel[i];
        y_bar_sum += y_bar_channel[i];
    }

    double x_bar_mean = x_bar_sum / 8;
    double y_bar_mean = y_bar_sum / 8;

    //save results in file couplingCorrection.txt

    ofstream couplingCorrection_file ("couplingCorrection.txt");//, ios::app);

    for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        couplingCorrection_file << std::noshowpos << channel_counter << "\t";
        couplingCorrection_file << std::fixed << std::showpos << std::setw(7) << std::setprecision(4) << x_bar_channel[channel_counter] << "\t\t";
        couplingCorrection_file << std::fixed << y_bar_channel[channel_counter] << "\n";
    } 

    couplingCorrection_file.close();

    //save mean in another file
    ofstream couplingCorrection_mean_file ("couplingCorrection_mean.txt");//, ios::app);


    couplingCorrection_mean_file << x_bar_mean << "\t" << y_bar_mean;

    couplingCorrection_mean_file.close();

    //Creating angles list
    Double_t pi = TMath::Pi();

    ofstream couplingCorrection_angles_file ("couplingCorrection_angles.txt");//, ios::app);

    double phi_ew_list[n_ch];
    for (int i = 0; i < n_ch; i++)
    {
        phi_ew_list[i] = atan(y_bar_channel[i] / x_bar_channel[i]) * 180 / pi;
        couplingCorrection_angles_file << i << "\t" << phi_ew_list[i] << "\n";
    }

    couplingCorrection_angles_file.close();

    cout << "Created .txt files." << endl;



     /*
    ___________GRAPHICS EXPORT ______________________
    */
    //style options
    /*gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1); // stats box size
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(16);
    gStyle->SetMarkerStyle(7);
    gStyle->SetMarkerSize(7);
    gStyle->SetMarkerColorAlpha(kBlack, 1.6);


    //the canvas in which all plots are saved
    TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 1000);
    masterCanvas->Divide(1, 1);
    masterCanvas->SetGrid();

    TGraph *angle_plot = new TGraph(8, channel_list, phi_ew_list);
    masterCanvas->cd(1);
    angle_plot->Draw("APE");
    masterCanvas->Print("Phi_ew_bar.pdf");*/
}