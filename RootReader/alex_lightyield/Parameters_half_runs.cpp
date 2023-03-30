/*
THIS IS FOR THE HALF RUNS!!!
*/


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
    //defining lists for variables (length 24 for my first measurement without mirror, 48 for half measurement)
    double run_number_list[48];
    double dist_max_list[48];
    double dist_mean_list[48];
    double mpv_list[48];
    double gauss_mean_list[48];

    //storing the mean temperature of all three thermometers over the period of 8 hours
    double temperature_list[48];

    //storing the Gaussian max error
    double gauss_max_error_y[48];

    //for the plots, we need shorter lists
    double run_number_list_part1[6];
    double dist_max_list_part1[6];
    double dist_mean_list_part1[6];
    double mpv_list_part1[6];
    double gauss_mean_list_part1[6];
    double temperature_list_part1[6];
    double gauss_max_error_y_part1[6];


    double run_number_list_part2[6];
    double dist_max_list_part2[6];
    double dist_mean_list_part2[6];
    double mpv_list_part2[6];
    double gauss_mean_list_part2[6];
    double temperature_list_part2[6];
    double gauss_max_error_y_part2[6];

    //for the errors
    double error_x[6] {0, 0, 0, 0, 0, 0};
    double dist_max_error_y[6] {2.67, 2.67, 2.67, 2.67, 2.67, 2.67};
    double gauss_max_error_y_part[12];

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
    ifstream temp_file ("Temperature.txt");

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

    /*
    //just for debugging, write Temperature.txt
    cout << "What I have read as Temperature.txt: \n" << endl;
    
    for (int i = 0; i<48; i++) cout << "temperature_list[" << i << "] = " << temperature_list[i] << endl;
    */

    //opening the file Values.txt and storing the values in the respective lists
    counter = 0;

    ifstream myfile ("Values.txt");
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

            //10 to 14 are dist_mean which is mpv
            help_str = "";
            for (int i=10; i<=14; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            dist_mean_list[counter] = help_num;

            //16 to 20 are gauss_max which is mpv
            help_str = "";
            for (int i=16; i<=20; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            mpv_list[counter] = help_num;

            //22 to 26 are gauss_mean
            help_str = "";
            for (int i=22; i<=26; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            gauss_mean_list[counter] = help_num;

            //28 to 31 are gauss_max_error
            help_str = "";
            for (int i=28; i<=31; i++) help_str=help_str+line[i];
            
            help_num = atof(help_str.c_str());
            gauss_max_error_y[counter] = help_num;


            counter++;           
        }
    }
    else cout << "Could not open values file." << endl;

    myfile.close();


    //just for debugging, write Values.txt
    /* cout << "\n \n \nWhat I have read as Values.txt:\n" << endl;
    for (int i = 0; i<48; i++) 
    {
        cout << run_number_list[i] << "\t" << dist_max_list[i] << "\t" << dist_mean_list[i] << "\t" << mpv_list[i] << "\t" << gauss_mean_list[i] << endl;
    }
    */

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
    masterCanvas->Divide(4, 3);
    masterCanvas->SetGrid();


    /*
    ----------------PCB C-------------------------------
    */  
    for (int i = 0; i<=5; i++)
    {
        run_number_list_part1[i] = run_number_list[2 * i];
        dist_max_list_part1[i] = dist_max_list[2 * i];
        dist_mean_list_part1[i] = dist_mean_list[2 * i];
        mpv_list_part1[i] = mpv_list[2 * i];
        gauss_mean_list_part1[i] = gauss_mean_list[2 * i];
        temperature_list_part1[i] = temperature_list[2 * i];
        gauss_max_error_y_part1[i] = gauss_max_error_y[2 * i];

        run_number_list_part2[i] = run_number_list[2 * i + 1];
        dist_max_list_part2[i] = dist_max_list[2 * i + 1];
        dist_mean_list_part2[i] = dist_mean_list[2 * i + 1];
        mpv_list_part2[i] = mpv_list[2 * i + 1];
        gauss_mean_list_part2[i] = gauss_mean_list[2 * i + 1];
        temperature_list_part2[i] = temperature_list[2 * i + 1];
        gauss_max_error_y_part2[i] = gauss_max_error_y[2 * i + 1];
    } 

    //just for debugging
    //for (int i = 0; i<6; i++)   cout << "dist_mean_list_part1[" << i << "] = " << dist_mean_list_part1[i] << endl;

    //correct for individual temperature values
    //cout << "\n\n\nWhat I have calculated as Correction factors:\n" << endl;
    /*for (int i = 0; i<=5; i++)
    {
        correction_factor = (4e6 - 0.2e6 * (temperature_list_part1[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part1[0] - 25));
        // cout << correction_factor << endl;

        dist_max_list_part1[i] = dist_max_list_part1[i] * correction_factor;
        dist_mean_list_part1[i] = dist_mean_list_part1[i] * correction_factor;
        mpv_list_part1[i] = mpv_list_part1[i] * correction_factor;
        gauss_mean_list_part1[i] = gauss_mean_list_part1[i] * correction_factor;



        correction_factor = (4e6 - 0.2e6 * (temperature_list_part2[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part2[0] - 25));
        // cout << correction_factor << endl;

        dist_max_list_part2[i] = dist_max_list_part2[i] * correction_factor;
        dist_mean_list_part2[i] = dist_mean_list_part2[i] * correction_factor;
        mpv_list_part2[i] = mpv_list_part2[i] * correction_factor;
        gauss_mean_list_part2[i] = gauss_mean_list_part2[i] * correction_factor;
    } 
    */


    //Distribution maximum
    TGraph* dist_max_plot_c1 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part1, error_x, dist_max_error_y);
    TGraph* dist_max_plot_c2 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part2, error_x, dist_max_error_y);
    dist_max_plot_c1->SetTitle("Distribution maximum for PCB C; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(1);
    dist_max_plot_c1->Draw("APE");
    dist_max_plot_c2->SetMarkerColorAlpha(kRed, 1.6);
    dist_max_plot_c2->Draw("SAME P");

    //Distribution mean
    /*TGraph* dist_mean_plot_c1 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part1);
    TGraph* dist_mean_plot_c2 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part2);
    dist_mean_plot_c1->SetTitle("Distribution mean for PCB C; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    dist_mean_plot_c1->Draw("APE");
    dist_mean_plot_c2->SetMarkerColorAlpha(kRed, 1.6);
    dist_mean_plot_c2->Draw("SAME P");*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_c1 = new TGraphErrors(6, run_number_list_part1, mpv_list_part1, error_x, gauss_max_error_y_part1);
    TGraph* mpv_plot_c2 = new TGraphErrors(6, run_number_list_part1, mpv_list_part2, error_x, gauss_max_error_y_part2);
    mpv_plot_c1->SetTitle("Gaussian maximum for PCB C; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    mpv_plot_c1->Draw("APE");
    mpv_plot_c2->SetMarkerColorAlpha(kRed, 1.6);
    mpv_plot_c2->Draw("SAME P");


    /*
    ----------------PCB J-------------------------------
    */
    for (int i = 0; i<=5; i++)
    {
        run_number_list_part1[i] = run_number_list[2 * (i + 6)];
        dist_max_list_part1[i] = dist_max_list[2 * (i + 6)];
        dist_mean_list_part1[i] = dist_mean_list[2 * (i + 6)];
        mpv_list_part1[i] = mpv_list[2 * (i + 6)];
        gauss_mean_list_part1[i] = gauss_mean_list[2 * (i + 6)];
        temperature_list_part1[i] = temperature_list[2 * (i + 6)];
        gauss_max_error_y_part1[i] = gauss_max_error_y[2 * (i + 6)];

        run_number_list_part2[i] = run_number_list[2 * (i + 6) + 1];
        dist_max_list_part2[i] = dist_max_list[2 * (i + 6) + 1];
        dist_mean_list_part2[i] = dist_mean_list[2 * (i + 6) + 1];
        mpv_list_part2[i] = mpv_list[2 * (i + 6) + 1];
        gauss_mean_list_part2[i] = gauss_mean_list[2 * (i + 6) + 1];
        temperature_list_part2[i] = temperature_list[2 * (i + 6) + 1];
        gauss_max_error_y_part2[i] = gauss_max_error_y[2 * (i + 6) + 1];
    } 

    /*
    //just for debugging
    for (int i = 0; i<6; i++)   cout << "dist_mean_list_part1[" << i << "] = " << dist_mean_list_part1[i] << endl;
    */


    //correct for individual temperature values
    /*for (int i = 0; i<=5; i++)
    {
        correction_factor = (4e6 - 0.2e6 * (temperature_list_part1[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part1[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part1[i] = dist_max_list_part1[i] * correction_factor;
        dist_mean_list_part1[i] = dist_mean_list_part1[i] * correction_factor;
        mpv_list_part1[i] = mpv_list_part1[i] * correction_factor;
        gauss_mean_list_part1[i] = gauss_mean_list_part1[i] * correction_factor;

        correction_factor = (4e6 - 0.2e6 * (temperature_list_part2[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part2[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part2[i] = dist_max_list_part2[i] * correction_factor;
        dist_mean_list_part2[i] = dist_mean_list_part2[i] * correction_factor;
        mpv_list_part2[i] = mpv_list_part2[i] * correction_factor;
        gauss_mean_list_part2[i] = gauss_mean_list_part2[i] * correction_factor;
    } */
    
    
    //distribution maximum
    TGraph* dist_max_plot_j1 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part1, error_x, dist_max_error_y);
    TGraph* dist_max_plot_j2 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part2, error_x, dist_max_error_y);
    dist_max_plot_j1->SetTitle("Distribution maximum for Juelich PCB (J); Run number; Pulse height [mV x ns]");
    masterCanvas->cd(2);
    dist_max_plot_j1->Draw("APE");
    dist_max_plot_j2->SetMarkerColorAlpha(kRed, 1.6);
    dist_max_plot_j2->Draw("SAME P");

    //Distribution mean
    /*TGraph* dist_mean_plot_j1 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part1);
    TGraph* dist_mean_plot_j2 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part2);
    dist_mean_plot_j1->SetTitle("Distribution mean for Juelich PCB (J); Run number; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    dist_mean_plot_j1->Draw("APE");
    dist_mean_plot_j2->SetMarkerColorAlpha(kRed, 1.6);
    dist_mean_plot_j2->Draw("SAME P");*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_j1 = new TGraphErrors(6, run_number_list_part1, mpv_list_part1, error_x, gauss_max_error_y_part1);
    TGraph* mpv_plot_j2 = new TGraphErrors(6, run_number_list_part1, mpv_list_part2, error_x, gauss_max_error_y_part2);
    mpv_plot_j1->SetTitle("Gaussian maximum for Juelich PCB (J); Run number; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    mpv_plot_j1->Draw("APE");
    mpv_plot_j2->SetMarkerColorAlpha(kRed, 1.6);
    mpv_plot_j2->Draw("SAME P");


    /*
    ----------------PCB D-------------------------------
    */

    for (int i = 0; i<=5; i++)
    {
        run_number_list_part1[i] = run_number_list[2 * (i + 12)];
        dist_max_list_part1[i] = dist_max_list[2 * (i + 12)];
        dist_mean_list_part1[i] = dist_mean_list[2 * (i + 12)];
        mpv_list_part1[i] = mpv_list[2 * (i + 12)];
        gauss_mean_list_part1[i] = gauss_mean_list[2 * (i + 12)];
        temperature_list_part1[i] = temperature_list[2 * (i + 12)];
        gauss_max_error_y_part1[i] = gauss_max_error_y[2 * (i + 12)];

        run_number_list_part2[i] = run_number_list[2 * (i + 12) + 1];
        dist_max_list_part2[i] = dist_max_list[2 * (i + 12) + 1];
        dist_mean_list_part2[i] = dist_mean_list[2 * (i + 12) + 1];
        mpv_list_part2[i] = mpv_list[2 * (i + 12) + 1];
        gauss_mean_list_part2[i] = gauss_mean_list[2 * (i + 12) + 1];
        temperature_list_part2[i] = temperature_list[2 * (i + 12) + 1];
        gauss_max_error_y_part2[i] = gauss_max_error_y[2 * (i + 12) + 1];
    } 

    //correct for individual temperature values
    /*for (int i = 0; i<=5; i++)
    {
        correction_factor = (4e6 - 0.2e6 * (temperature_list_part1[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part1[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part1[i] = dist_max_list_part1[i] * correction_factor;
        dist_mean_list_part1[i] = dist_mean_list_part1[i] * correction_factor;
        mpv_list_part1[i] = mpv_list_part1[i] * correction_factor;
        gauss_mean_list_part1[i] = gauss_mean_list_part1[i] * correction_factor;

        correction_factor = (4e6 - 0.2e6 * (temperature_list_part2[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part2[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part2[i] = dist_max_list_part2[i] * correction_factor;
        dist_mean_list_part2[i] = dist_mean_list_part2[i] * correction_factor;
        mpv_list_part2[i] = mpv_list_part2[i] * correction_factor;
        gauss_mean_list_part2[i] = gauss_mean_list_part2[i] * correction_factor;
    } */
    
    
    //distribution maximum
    TGraph* dist_max_plot_d1 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part1, error_x, dist_max_error_y);
    TGraph* dist_max_plot_d2 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part2, error_x, dist_max_error_y);
    dist_max_plot_d1->SetTitle("Distribution maximum for PCB D; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(3);
    dist_max_plot_d1->Draw("APE");
    dist_max_plot_d2->SetMarkerColorAlpha(kRed, 1.6);
    dist_max_plot_d2->Draw("SAME P");

    //Distribution mean
    /*TGraph* dist_mean_plot_d1 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part1);
    TGraph* dist_mean_plot_d2 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part2);
    dist_mean_plot_d1->SetTitle("Distribution mean for PCB D; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    dist_mean_plot_d1->Draw("APE");
    dist_mean_plot_d2->SetMarkerColorAlpha(kRed, 1.6);
    dist_mean_plot_d2->Draw("SAME P");*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_d1 = new TGraphErrors(6, run_number_list_part1, mpv_list_part1, error_x, gauss_max_error_y_part1);
    TGraph* mpv_plot_d2 = new TGraphErrors(6, run_number_list_part1, mpv_list_part2, error_x, gauss_max_error_y_part2);
    mpv_plot_d1->SetTitle("Gaussian maximum for PCB D; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    mpv_plot_d1->Draw("APE");
    mpv_plot_d2->SetMarkerColorAlpha(kRed, 1.6);
    mpv_plot_d2->Draw("SAME P");


    /*
    ----------------PCB D but without new optical gel---------------------
    */

    for (int i = 0; i<=5; i++)
    {
        run_number_list_part1[i] = run_number_list[2 * (i + 18)];
        dist_max_list_part1[i] = dist_max_list[2 * (i + 18)];
        dist_mean_list_part1[i] = dist_mean_list[2 * (i + 18)];
        mpv_list_part1[i] = mpv_list[2 * (i + 18)];
        gauss_mean_list_part1[i] = gauss_mean_list[2 * (i + 18)];
        temperature_list_part1[i] = temperature_list[2 * (i + 18)];
        gauss_max_error_y_part1[i] = gauss_max_error_y[2 * (i + 18)];

        run_number_list_part2[i] = run_number_list[2 * (i + 18) + 1];
        dist_max_list_part2[i] = dist_max_list[2 * (i + 18) + 1];
        dist_mean_list_part2[i] = dist_mean_list[2 * (i + 18) + 1];
        mpv_list_part2[i] = mpv_list[2 * (i + 18) + 1];
        gauss_mean_list_part2[i] = gauss_mean_list[2 * (i + 18) + 1];
        temperature_list_part2[i] = temperature_list[2 * (i + 18) + 1];
        gauss_max_error_y_part2[i] = gauss_max_error_y[2 * (i + 18) + 1];
    } 
    //just for debugging
    //for (int i = 0; i<6; i++)   cout << "dist_mean_list_part1[" << i << "] = " << dist_mean_list_part1[i] << endl;


    //correct for individual temperature values
    /*for (int i = 0; i<=5; i++)
    {
        correction_factor = (4e6 - 0.2e6 * (temperature_list_part1[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part1[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part1[i] = dist_max_list_part1[i] * correction_factor;
        dist_mean_list_part1[i] = dist_mean_list_part1[i] * correction_factor;
        mpv_list_part1[i] = mpv_list_part1[i] * correction_factor;
        gauss_mean_list_part1[i] = gauss_mean_list_part1[i] * correction_factor;

        correction_factor = (4e6 - 0.2e6 * (temperature_list_part2[i] - 25)) / (4e6 - 0.2e6 * (temperature_list_part2[0] - 25));
        //cout << correction_factor << endl;

        dist_max_list_part2[i] = dist_max_list_part2[i] * correction_factor;
        dist_mean_list_part2[i] = dist_mean_list_part2[i] * correction_factor;
        mpv_list_part2[i] = mpv_list_part2[i] * correction_factor;
        gauss_mean_list_part2[i] = gauss_mean_list_part2[i] * correction_factor;
    } */
    
    //just for debugging
    for (int i = 0; i<6; i++)   cout << "dist_max_list_part2[" << i << "] = " << dist_max_list_part2[i] << endl;

    
    //distribution maximum
    TGraph* dist_max_plot_d21 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part1, error_x, dist_max_error_y);
    TGraph* dist_max_plot_d22 = new TGraphErrors(6, run_number_list_part1, dist_max_list_part2, error_x, dist_max_error_y);
    dist_max_plot_d21->SetTitle("Distribution maximum for PCB D without new optical gel; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(4);
    dist_max_plot_d21->GetYaxis()->SetLimits(38, 47);
    dist_max_plot_d21->Draw("APE");
    dist_max_plot_d22->SetMarkerColorAlpha(kRed, 1.6);
    dist_max_plot_d22->Draw("SAME P");

    //Distribution mean
    /*TGraph* dist_mean_plot_d21 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part1);
    TGraph* dist_mean_plot_d22 = new TGraphErrors(6, run_number_list_part1, dist_mean_list_part2);
    dist_mean_plot_d21->SetTitle("Distribution mean for PCB D without new optical gel; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    dist_mean_plot_d21->Draw("APE");
    dist_mean_plot_d22->SetMarkerColorAlpha(kRed, 1.6);
    dist_mean_plot_d22->Draw("SAME P");*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_d21 = new TGraphErrors(6, run_number_list_part1, mpv_list_part1, error_x, gauss_max_error_y_part1);
    TGraph* mpv_plot_d22 = new TGraphErrors(6, run_number_list_part1, mpv_list_part2, error_x, gauss_max_error_y_part2);
    mpv_plot_d21->SetTitle("Gaussian maximum for PCB D without new optical gel; Run number; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    mpv_plot_d21->Draw("APE");
    mpv_plot_d22->SetMarkerColorAlpha(kRed, 1.6);
    mpv_plot_d22->Draw("SAME P");


    //print everything to a pdf

    masterCanvas->Print("Parameters without temperature correction.pdf");

}

