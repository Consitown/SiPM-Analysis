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
#include <TLatex.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TPave.h>
#include <TPaveText.h>


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

int main()
{
    //defining lists for variables (length 24 for my first measurement without mirror)
    double run_number_list[24];
    double dist_max_list[24];
    double dist_mean_list[24];
    double gauss_max_list[24];
    double mpv_list[24];
    double mpv_err_list[24];

    //storing the mean temperature of all three thermometers over the period of 8 hours
    double temperature_list[24];

    //storing the Gaussian max error
    double gauss_max_error_y[24];

    //for the plots, we need shorter lists
    double run_number_list_part[6];
    double dist_max_list_part[6];
    double dist_mean_list_part[6];
    double gauss_max_list_part[6];
    double mpv_list_part[6];
    double mpv_err_list_part[6];
    double temperature_list_part[6];

    //List for errors (just bin width)
    double dist_max_error_x[6] {0, 0, 0, 0, 0, 0};
    double dist_max_error_y[6] {2.67, 2.67, 2.67, 2.67, 2.67, 2.67};
    double gauss_max_error_y_part[6];

    //Variable for temperature correction
    double correction_factor;

    //variable to read single lines from the file
    string line;

    //variable to store numbers as strings first
    string help_str;

    //variable to convert from help_str to number
    double help_num;

    int counter = 0;

    //opening the file Temperatures.txt and storing the results in the list temperature_list
    ifstream temp_file ("Temperature dependence/no_mirror_jan/Temperature.txt");
    if (temp_file.is_open())
    {
        while (getline (temp_file, line))
        {
            help_str = "";

            for (int i=4; i<=7; i++) help_str = help_str + line[i];

            help_num = atof(help_str.c_str());
            temperature_list[counter] = help_num;

            counter++;
        }
    }
    else cout << "Could not open temperature file." << endl;

    temp_file.close();
    
    //just for debugging, write Temperature.txt
    //cout << "What I have read as Temperature.txt: \n" << endl;
    //for (int i = 0; i<=23; i++) cout << "temperature_list[" << i << "] = " << temperature_list[i] << endl;

    //opening the file Values.txt and storing the values in the respective lists
    counter = 0;

    ifstream myfile ("Temperature dependence/no_mirror_jan/sum_values.txt");
    if (myfile.is_open())
    {
        while (getline (myfile, line))
        { 

            //0 to 2 are run_number
            help_str = "";
            for (int i=0; i<=2; i++) help_str = help_str + line[i];
            
            help_num = atof(help_str.c_str());
            run_number_list[counter] = help_num;

            //4 to 8 are dist_max
            help_str = "";
            for (int i=4; i<=8; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            dist_max_list[counter] = help_num;

            //10 to 14 are dist_mean
            help_str = "";
            for (int i=10; i<=14; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            dist_mean_list[counter] = help_num;

            //22 to 26 are langau MPV
            help_str = "";
            for (int i=22; i<=26; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            mpv_list[counter] = help_num;

            //39 to 42 are langau MPV
            help_str = "";
            for (int i=28; i<=31; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            mpv_err_list[counter] = help_num;

            counter++;           
        }
    }
    else cout << "Could not open values file." << endl;

    myfile.close();
    for (int i = 0; i<=23; i++) 
    {
        cout << run_number_list[i] << "\t" << mpv_list[i] << "\t" << mpv_err_list[i] << endl;
    }

    //just for debugging, write Values.txt
    /*cout << "\n \n \nWhat I have read as sum_values.txt:\n" << endl;
    for (int i = 0; i<=23; i++) 
    {
        cout << run_number_list[i] << "\t" << dist_max_list[i] << "\t" << dist_mean_list[i] << "\t" << mpv_list[i] << "\t" << mpv_err_list[i] << endl;
    }*/


    //plotting the results

    //style options
    gStyle->SetOptStat(0);
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
    masterCanvas->Divide(2, 4);
    masterCanvas->SetGrid();


    /*
    ----------------PCB C-------------------------------
    */
    for (int i = 0; i<=5; i++)
    {
        run_number_list_part[i] = run_number_list[i];
        dist_max_list_part[i] = dist_max_list[i];
        dist_mean_list_part[i] = dist_mean_list[i];
        temperature_list_part[i] = temperature_list[i];
        mpv_list_part[i] = mpv_list[i];
        mpv_err_list_part[i] = mpv_err_list[i];
    } 
 


    //Distribution maximum
    TGraph* dist_max_plot_c = new TGraphErrors(6, temperature_list_part, dist_max_list_part, dist_max_error_x, dist_max_error_y);
    dist_max_plot_c->SetTitle("Distribution maximum for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(1);
    dist_max_plot_c->Draw("APE");

    //Distribution mean
    /*TGraph* dist_mean_plot_c = new TGraphErrors(6, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_c->SetTitle("Distribution mean for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    dist_mean_plot_c->Draw("APE");*/

    //Gaussian maximum
    /*TGraph* mpv_plot_c = new TGraphErrors(6, temperature_list_part, gauss_max_list_part, dist_max_error_x, gauss_max_error_y_part);
    mpv_plot_c->SetTitle("Gaussian maximum for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    mpv_plot_c->Draw("APE");*/

    //Langau MPV
    TGraph* mpv_plot_c = new TGraphErrors(6, temperature_list_part, mpv_list_part, dist_max_error_x, mpv_err_list_part);
    mpv_plot_c->SetTitle("Landau #times Gauss MPV for PCB C; Temperature in ^{#circ}C; MPV [mV x ns]");
    masterCanvas->cd(2);
    mpv_plot_c->Draw("APE");


    /*
    ----------------PCB J-------------------------------
    */
    
    for (int i = 0; i<=5; i++)
    {
        run_number_list_part[i] = run_number_list[i + 6];
        dist_max_list_part[i] = dist_max_list[i + 6];
        dist_mean_list_part[i] = dist_mean_list[i + 6];
        temperature_list_part[i] = temperature_list[i + 6];
        mpv_list_part[i] = mpv_list[i + 6];
        mpv_err_list_part[i] = mpv_err_list[i + 6];
    } 

    
    //distribution maximum
    TGraph* dist_max_plot_j = new TGraphErrors(6, temperature_list_part, dist_max_list_part, dist_max_error_x, dist_max_error_y);
    dist_max_plot_j->SetTitle("Distribution maximum for PCB J; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(3);
    dist_max_plot_j->Draw("APE");

    //Distribution mean
    /*TGraph* dist_mean_plot_j = new TGraphErrors(6, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_j->SetTitle("Distribution mean for Juelich PCB (J); Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    dist_mean_plot_j->Draw("APE");*/

    //Gaussian maximum/MPV
    /*TGraph* mpv_plot_j = new TGraphErrors(6, temperature_list_part, gauss_max_list_part, dist_max_error_x, gauss_max_error_y_part);
    mpv_plot_j->SetTitle("Gaussian maximum for Juelich PCB (J); Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    mpv_plot_j->Draw("APE");*/

    //Langau MPV
    TGraph* mpv_plot_j = new TGraphErrors(6, temperature_list_part, mpv_list_part, dist_max_error_x, mpv_err_list_part);
    mpv_plot_j->SetTitle("Landau #times Gauss MPV for PCB J; Temperature in ^{#circ}C; MPV [mV x ns]");
    masterCanvas->cd(4);
    mpv_plot_j->Draw("APE");


    /*
    ----------------PCB D-------------------------------
    */

    for (int i = 0; i<=5; i++)
    {
        run_number_list_part[i] = run_number_list[i + 12];
        dist_max_list_part[i] = dist_max_list[i + 12];
        dist_mean_list_part[i] = dist_mean_list[i + 12];
        temperature_list_part[i] = temperature_list[i + 12];
        mpv_list_part[i] = mpv_list[i + 12];
        mpv_err_list_part[i] = mpv_err_list[i + 12];
    } 

    
    //distribution maximum
    TGraph* dist_max_plot_d = new TGraphErrors(6, temperature_list_part, dist_max_list_part, dist_max_error_x, dist_max_error_y);
    dist_max_plot_d->SetTitle("Distribution maximum for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    dist_max_plot_d->Draw("APE");

    //Distribution mean
    /*TGraph* dist_mean_plot_d = new TGraphErrors(6, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_d->SetTitle("Distribution mean for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    dist_mean_plot_d->Draw("APE");*/

    //Gaussian maximum/MPV
    /*TGraph* mpv_plot_d = new TGraphErrors(6, temperature_list_part, gauss_max_list_part, dist_max_error_x, gauss_max_error_y_part);
    mpv_plot_d->SetTitle("Gaussian maximum for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    mpv_plot_d->Draw("APE");*/

    //Langau MPV
    TGraph* mpv_plot_d = new TGraphErrors(6, temperature_list_part, mpv_list_part, dist_max_error_x, mpv_err_list_part);
    mpv_plot_d->SetTitle("Landau #times Gauss MPV for PCB D; Temperature in ^{#circ}C; MPV [mV x ns]");
    masterCanvas->cd(6);
    mpv_plot_d->Draw("APE");



    /*
    ----------------PCB D but without new optical gel---------------------
    */

    double run_number_list_part2[6];
    double dist_max_list_part2[6];
    double dist_mean_list_part2[6];
    double temperature_list_part2[6];
    double mpv_list_part2[6];
    double mpv_err_list_part2[6];

    double dist_max_error_x2[6] {0, 0, 0, 0, 0, 0};
    double dist_max_error_y2[6] {2.67, 2.67, 2.67, 2.67, 2.67, 2.67};

    for (int i = 0; i<=5; i++)
    {
        run_number_list_part2[i] = run_number_list[i + 18];
        dist_max_list_part2[i] = dist_max_list[i + 18];
        dist_mean_list_part2[i] = dist_mean_list[i + 18];
        temperature_list_part2[i] = temperature_list[i + 18];
        mpv_list_part2[i] = mpv_list[i + 18];
        mpv_err_list_part2[i] = mpv_err_list[i + 18];
    } 

    
    //distribution maximum
    TGraph* dist_max_plot_d2 = new TGraphErrors(6, temperature_list_part2, dist_max_list_part2, dist_max_error_x2, dist_max_error_y2);
    dist_max_plot_d2->SetTitle("Distribution maximum for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    dist_max_plot_d2->Draw("APE");

    //Distribution mean
    /*TGraph* dist_mean_plot_d2 = new TGraphErrors(4, temperature_list_part2, dist_mean_list_part2);
    dist_mean_plot_d2->SetTitle("Distribution mean for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    dist_mean_plot_d2->Draw("APE");*/

    //Gaussian maximum
    /*TGraph* mpv_plot_d2 = new TGraphErrors(4, temperature_list_part2, gauss_max_list_part2, dist_max_error_x2, gauss_max_error_y_part2);
    mpv_plot_d2->SetTitle("Gaussian maximum for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    mpv_plot_d2->Draw("APE");*/

     //Langau MPV
    TGraph* mpv_plot_d2 = new TGraphErrors(6, temperature_list_part2, mpv_list_part2, dist_max_error_x2, mpv_err_list_part2);
    mpv_plot_d2->SetTitle("Landau #times Gauss MPV for PCB D without new optical gel; Temperature in ^{#circ}C; MPV [mV x ns]");
    masterCanvas->cd(8);
    mpv_plot_d2->Draw("APE");


    //print everything to a pdf

    masterCanvas->Print("Temperature dependence/no_mirror_jan/Correlation plot.pdf");
}

