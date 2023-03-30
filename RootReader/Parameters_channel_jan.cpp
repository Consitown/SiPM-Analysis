//original 
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


//C, C++
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <ctime>

using namespace std;

int main()
{
    // defining some parameters for the program
    char algo = 'j'; //'j' for Jan, 'c' for Christian
    char pcb = 'n'; // for the PCB: 'c', 'd' (D with new optical gel) or 'n' (D without new optical gel) or 'j'
    bool with_black_foil = false; //true for with black foil
    bool plot_mean_mpv = true; //plots the mean value of the MPV per channel

    int n_ch_max; //number of channels, changes as long as there is nothing in channel 0 for Jan
    int n_r = 6; //number of runs
    int n_m = 1; //number of measurements (1 for one, 2 for with and without reflector, 3 for with and without reflector and black end)

    string algo_name;
    if (algo == 'c') algo_name = "Christian";
    if (algo == 'j') algo_name = "Jan";

    string pcb_name;
    if (pcb == 'c') pcb_name = "PCB C";
    if (pcb == 'd') pcb_name = "PCB D";
    if (pcb == 'n') pcb_name = "PCB D without new optical gel";
    if (pcb == 'j') pcb_name = "PCB J";

    n_ch_max = 8;


    cout << "Calculating for " << algo_name << "'s algorithm and " << pcb_name << endl;


    //defining lists for variables

    Double_t run_number_list[n_r][n_ch_max]; 
    Double_t channel_list[n_r][n_ch_max]; 
    Double_t mean_list[n_m][n_r][n_ch_max]; //first index is 0=without, 1=with reflector, 2 for black end, second index is run, third index is channel
    Double_t mean_err_list[n_m][n_r][n_ch_max];
    Double_t trunc_mean_list[n_m][n_r][n_ch_max];
    Double_t trunc_mean_err_list[n_m][n_r][n_ch_max];
    Double_t poisson_exp_list[n_m][n_r][n_ch_max];
    Double_t poisson_exp_err_list[n_m][n_r][n_ch_max];
    Double_t langau_mpv_list[n_m][n_r][n_ch_max];
    Double_t langau_mpv_err_list[n_m][n_r][n_ch_max];
    Double_t error_x[n_ch_max] {0};

    Double_t calib_values[n_ch_max]; //for saving the Lambda_pixel values, so that the light yield results can be converted to PEs
    Double_t mean_mpv_values[n_m][n_ch_max]; //shows the mean MPV value per channel
    Double_t mean_mpv_err_values[n_m][n_ch_max];

    Double_t plot_channel_list[n_ch_max]; //allows the dots and triangles in the plot to be shifted slightly from the original poisition for better visibility

    //variable to read single lines from the file
    string line;

    //variable to store numbers as strings first
    string help_str;

    //variable to convert from help_str to number
    double help_num;

    int line_counter = 0;

    string file_path = "/home/alex/Dokumente/Studium/RootReader/RootAnalysis/LightYield/Channel parameters/cuts/Langau/"; //change if other fit

    if (with_black_foil) file_path += "with black foil/";
    else file_path += pcb_name + "/";
    file_path += algo_name + "/"; //leads to the right folder

    string mean_file_path = file_path + "Channel_values_mean.txt";

    //for mean file
    ifstream mean_file (mean_file_path.c_str());
    if (mean_file.is_open())
    {
        while (getline (mean_file, line))
        { 
            // for mean
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
            {
                help_str = "";
                
                if (algo == 'j') for (int i=6; i<=10; i++) help_str = help_str + line[i + n_ch * 13]; //i=6-10 for Jan, n_ch*13
                if (algo == 'c') for (int i=6; i<=12; i++) help_str = help_str + line[i + n_ch * 16]; //i=6-12 for Christian, n_ch*16
            
                help_num = atof(help_str.c_str());

                mean_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;

                //cout << "Ich packe " << help_num << " nach mean_list[" << line_counter / n_r << "][" << line_counter % n_r << "][" <<  n_ch << "]" << endl;
                
            }

            // for mean_err
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++) 
            {
                help_str = "";
                if (algo == 'j') for (int i=12; i<=15; i++) help_str = help_str + line[i + n_ch * 13]; //i=12-15 for Jan, n_ch*13
                if (algo == 'c') for (int i=14; i<=18; i++) help_str = help_str + line[i + n_ch * 16]; //i=14-18 for Christian, n_ch*16
            
                help_num = atof(help_str.c_str());

                mean_err_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
            }

            line_counter++;    
        }
    }
    else cout << "Could not open truncated mean values file." << endl;

    mean_file.close();


    line_counter = 0;

    //for trunc mean file, only do it for Jan's algorithm
    if (algo == 'j')
    {
        string trunc_mean_file_path = file_path + "Channel_values_trunc_mean.txt";
        ifstream trunc_mean_file (trunc_mean_file_path.c_str());
        if (trunc_mean_file.is_open())
        {
            while (getline (trunc_mean_file, line))
            { 
                // for trunc mean
                for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
                {
                    help_str = "";
                    for (int i=6; i<=10; i++) help_str = help_str + line[i + n_ch * 13];
                
                    help_num = atof(help_str.c_str());
                    
                    trunc_mean_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
                }

                // for trunc mean_err
                for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
                {
                    help_str = "";
                    for (int i=12; i<=15; i++) help_str = help_str + line[i + n_ch * 13];
                
                    help_num = atof(help_str.c_str());

                    trunc_mean_err_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
                }

                line_counter++;    
            }
        }
        else cout << "Could not open mean values file." << endl;

        trunc_mean_file.close();
    }

    //for Poisson file
    /*line_counter = 0;
    string poisson_file_path = file_path + "Channel_values_poisson.txt";
    ifstream poisson_file (poisson_file_path.c_str());
    if (poisson_file.is_open())
    {
        while (getline (poisson_file, line))
        { 
             //0 to 2 are run_number
            help_str = "";
            for (int i=0; i<=2; i++) help_str = help_str + line[i];
            
            help_num = atof(help_str.c_str());
            
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
            {
                if (line_counter < n_r)
                {
                    run_number_list[line_counter % n_r][n_ch] = help_num;
                    channel_list[line_counter % n_r][n_ch] = n_ch;
                }
            }

            // for poisson exp value
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++) 
            {
                help_str = "";
                if (algo == 'j') for (int i=6; i<=10; i++) help_str = help_str + line[i + n_ch * 13]; //i=6-10 for Jan, n_ch*13
                if (algo == 'c') for (int i=6; i<=12; i++) help_str = help_str + line[i + n_ch * 15]; //i=6-12 for Christian, n_ch*15
            
                help_num = atof(help_str.c_str());

                poisson_exp_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
            }

            // for poisson_exp_err
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
            {
                help_str = "";
                if (algo == 'j') for (int i=12; i<=15; i++) help_str = help_str + line[i + n_ch * 13]; //i=12-15 for Jan, n_ch*13
                if (algo == 'c') for (int i=14; i<=17; i++) help_str = help_str + line[i + n_ch * 15]; //i=14-17 for Christian, n_ch*15
            
                help_num = atof(help_str.c_str());
                poisson_exp_err_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
            }

            line_counter++;    
        }
    }
    else cout << "Could not open poisson values file." << endl;

    poisson_file.close();*/

    //for Langau file
    line_counter = 0;
    string langau_file_path = file_path + "Channel_values_langau.txt";
    ifstream langau_file (langau_file_path.c_str());
    if (langau_file.is_open())
    {
        while (getline (langau_file, line))
        { 
             //0 to 2 are run_number
            help_str = "";
            for (int i=0; i<=2; i++) help_str = help_str + line[i];
            
            help_num = atof(help_str.c_str());
            
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
            {
                if (line_counter < n_r)
                {
                    run_number_list[line_counter % n_r][n_ch] = help_num;
                    channel_list[line_counter % n_r][n_ch] = n_ch;
                }
            }

            // for Langau MPV/peak
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++) 
            {
                help_str = "";
                if (algo == 'j') for (int i=6; i<=10; i++) help_str = help_str + line[i + n_ch * 13]; //i=6-10 for Jan, n_ch*13
                if (algo == 'c') for (int i=6; i<=12; i++) help_str = help_str + line[i + n_ch * 15]; //i=6-12 for Christian, n_ch*15
            
                help_num = atof(help_str.c_str());

                langau_mpv_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
            }

            // for Langau MPV error
            for (int n_ch = 0; n_ch<n_ch_max; n_ch++)
            {
                help_str = "";
                if (algo == 'j') for (int i=12; i<=15; i++) help_str = help_str + line[i + n_ch * 13]; //i=12-15 for Jan, n_ch*13
                if (algo == 'c') for (int i=14; i<=17; i++) help_str = help_str + line[i + n_ch * 15]; //i=14-17 for Christian, n_ch*15
            
                help_num = atof(help_str.c_str());
                langau_mpv_err_list[line_counter / n_r][line_counter % n_r][n_ch] = help_num;
            }

            line_counter++;    
        }
    }
    else cout << "Could not open Langau values file." << endl;

    langau_file.close();

    //for calibration values file
    line_counter = 0;
    string calib_file_path = "/home/alex/Dokumente/Studium/Kalibrierungsmessungen/calibration_values ";
    calib_file_path += pcb_name;
    calib_file_path += ".txt";

    ifstream calib_file(calib_file_path.c_str());

    cout << "Opening calibration file at " << calib_file_path << endl;

    if (calib_file.is_open())
    {
        while (getline (calib_file, line))
        { 
            if (line_counter > 0)
            {
                help_str = "";
                for (int i=2; i<=6; i++) help_str = help_str + line[i];
                
                help_num = atof(help_str.c_str());

                calib_values[line_counter - 1] = help_num;
            }
            line_counter++;
        }
    }
    else cout << "Could not open calibration values file." << endl;
    
    //debugging
    //for (int ch = 0; ch<n_ch_max; ch++) cout << calib_values[ch] << endl;

    //dividing the light yield values by the calibration values to get the PE number
    for (int measurement = 0; measurement < n_m; measurement++)
    {
        for (int run = 0; run < n_r; run++)
        {
            for (int channel = 0; channel < n_ch_max; channel++)
            {
                mean_list[measurement][run][channel] /= calib_values[channel];
                trunc_mean_list[measurement][run][channel] /= calib_values[channel];
                langau_mpv_list[measurement][run][channel] /= calib_values[channel];

                mean_err_list[measurement][run][channel] /= calib_values[channel];
                trunc_mean_err_list[measurement][run][channel] /= calib_values[channel];
                langau_mpv_err_list[measurement][run][channel] /= calib_values[channel];
                //cout << mean_list[measurement][run][channel] << "\t";
            }
        }
    }

    //calculating the mean MPV values and the STD for each channel
    double sum, mean, var, sd;
    for (int measurement = 0; measurement < n_m; measurement++)
    {
        for (int channel = 0; channel < n_ch_max; channel++)
        {
            //for the mean
            sum = 0;
            if (with_black_foil && measurement == 1) n_r = 3;
            for (int run = 0; run < n_r; run++)
            {
                sum += langau_mpv_list[measurement][run][channel]; //is then the sum of all MPV of one channel and one measurement
            }

            mean = sum / n_r; //this is the mean of the MPV values
            mean_mpv_values[measurement][channel] = mean; 

            //for the variance
            var = 0;
            for (int run = 0; run < n_r; run++)
            {
                var += (langau_mpv_list[measurement][run][channel] - mean) * (langau_mpv_list[measurement][run][channel] - mean);
            }
            var /= n_r; //this is the variance

            sd = sqrt(var); //this is the standard deviation
            mean_mpv_err_values[measurement][channel] = sd; 
        }
    }
    n_r = 6;
    cout << "\n\n";
    float change, change_mean = 0, change_var = 0, change_std = 0;
    for (int channel = 0; channel < n_ch_max; channel++)
    {
        change = (mean_mpv_values[1][channel] / mean_mpv_values[0][channel] - 1) * 100;
        cout << "Change in light yield in percent (ch. " << channel << "):\t" << change << endl;
        change_mean += change;
    }

    change_mean /= n_ch_max;

    for (int channel = 0; channel < n_ch_max; channel++)
    {
        change = (mean_mpv_values[1][channel] / mean_mpv_values[0][channel] - 1) * 100;
        change_var += (change_mean - change) * (change_mean - change);
    }

    change_var /= n_ch_max;
    change_std = sqrt(change_var);
    cout << "\n\nThis is an average change of "<< change_mean << " pm " << change_std << "\n\n" << endl;
    //for (int i = 0; i<6; i++) cout << langau_mpv_list[1][i][0] << endl;
    //cout << "\n\n mean = " << mean_mpv_values[1][0] << endl;
    //cout << "\n\n std = " << mean_mpv_err_values[1][0] << endl;
    /*
    ___________GRAPHICS EXPORT ______________________
    */

    //style options
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1); // stats box size
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(16);
    gStyle->SetMarkerStyle(8);
    gStyle->SetMarkerSize(1.2);
    gStyle->SetMarkerColorAlpha(kBlack, 1.6);


    //the canvas in which all plots are saved
    TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 1350);
    masterCanvas->Divide(1, 3);
    masterCanvas->SetGrid();

    //the multigraph that contains the mean values
    TMultiGraph *mean_mg = new TMultiGraph();
    mean_mg->SetTitle(Form("Mean values of N_{PE} distribution for %s; Channel number; N_{PE}", pcb_name.c_str()));

    //the multigraph that contains the trunc mean values
    TMultiGraph *trunc_mean_mg = new TMultiGraph();
    trunc_mean_mg->SetTitle(Form("Truncated mean (0 to 150) of N_{PE} distribution for %s; Channel number; N_{PE}", pcb_name.c_str()));

    //the multigraph that contains the poisson exp values
    //TMultiGraph *exp_mg = new TMultiGraph();
    //exp_mg->SetTitle(Form("Poisson exp. values for %s; Channel number; Poisson exp. value", pcb_name.c_str()));

    //the multigraph that contains the Langau MPV values
    TMultiGraph *mpv_mg = new TMultiGraph();
    mpv_mg->SetTitle(Form("Landau*Gauss MPV of N_{PE} distribution for %s; Channel number; N_{PE}", pcb_name.c_str()));

    //the multigraph that contains the poisson*Gauss max values
    TMultiGraph *max_mg = new TMultiGraph();
    max_mg->SetTitle(Form("Extreme point of the Poisson*Gaussfit for %s; Channel number; x-value of maximum", pcb_name.c_str()));

    /*
    ____________LEGEND________________
    */
    TLegend *leg_mean = new TLegend(0.75,0.1,0.9,0.2);
    TLegend *leg_trunc = new TLegend(0.75,0.1,0.9,0.2);
    TLegend *leg_mpv = new TLegend(0.65,0.1,0.9,0.27);
    leg_mean->SetTextSize(0.035);
    leg_trunc->SetTextSize(0.035);
    leg_mpv->SetTextSize(0.035);
    //plots without reflector
    //dots should appear slightly left of the actual channel number
    for (int i = 0; i < 8; i++) plot_channel_list[i] = i - 0.13;

    for (int run = 0; run<n_r; run++)
    {
        TGraph* mean_plot = new TGraphErrors(n_ch_max, plot_channel_list, mean_list[0][run], error_x, mean_err_list[0][run]);
        mean_plot->SetMarkerColor(3 * run + 51);
        if (run == 1)
        {
            if (with_black_foil)
            {
                leg_mean->AddEntry(mean_plot, Form("Without black foil"),"p");
                leg_trunc->AddEntry(mean_plot, Form("Without black foil"),"p");
                leg_mpv->AddEntry(mean_plot, Form("Without black foil"),"p");
            }
            else
            {
                leg_mean->AddEntry(mean_plot, Form("Without reflector"),"p");
                leg_trunc->AddEntry(mean_plot, Form("Without reflector"),"p");
                leg_mpv->AddEntry(mean_plot, Form("Without reflector"),"p");
            }
        }
        mean_mg->Add(mean_plot);
        

        if (algo == 'j')
        {
            TGraph* trunc_mean_plot = new TGraphErrors(n_ch_max, plot_channel_list, trunc_mean_list[0][run], error_x, trunc_mean_err_list[0][run]);
            trunc_mean_plot->SetMarkerColor(3 * run + 51);
            trunc_mean_mg->Add(trunc_mean_plot);

            //for Poisson
            /*TGraph* exp_plot = new TGraphErrors(n_ch_max, channel_list[run], poisson_exp_list[0][run], error_x, poisson_exp_err_list[0][run]);
            exp_plot->SetTitle("Poisson exp. value without reflector");
            exp_plot->SetMarkerColor(3 * run + 51);
            exp_mg->Add(exp_plot);*/

            //for Langau
            TGraph* mpv_plot = new TGraphErrors(n_ch_max, plot_channel_list, langau_mpv_list[0][run], error_x, langau_mpv_err_list[0][run]);
            mpv_plot->SetMarkerColor(3 * run + 51);
            mpv_mg->Add(mpv_plot);
        }

        if (algo == 'c')
        {
            TGraph* max_plot = new TGraphErrors(n_ch_max, plot_channel_list, poisson_exp_list[0][run], error_x, poisson_exp_err_list[0][run]);
            max_plot->SetMarkerColor(3 * run + 50);
            max_plot->SetMarkerSize(2);
            max_mg->Add(max_plot);
        }
    }

    //plots with reflector
    //dots should appear slightly left of the actual channel number
    if (n_m > 1){
        for (int i = 0; i < 8; i++) plot_channel_list[i] = i + 0.13;
        if (with_black_foil) n_r = 3;
        for (int run = 0; run<n_r; run++)
        {
            TGraph* mean_plot = new TGraphErrors(n_ch_max, plot_channel_list, mean_list[1][run], error_x, mean_err_list[1][run]);
            mean_plot->SetMarkerStyle(22);
            mean_plot->SetMarkerColor(3 * run + 80);
            if (run == 1)
            {
                if (with_black_foil)
                {
                    leg_mean->AddEntry(mean_plot, Form("With black foil"),"p");
                    leg_trunc->AddEntry(mean_plot, Form("With black foil"),"p");
                    leg_mpv->AddEntry(mean_plot, Form("With black foil"),"p");
                }
                else
                {
                    leg_mean->AddEntry(mean_plot, Form("With reflector"),"p");
                    leg_trunc->AddEntry(mean_plot, Form("With reflector"),"p");
                    leg_mpv->AddEntry(mean_plot, Form("With reflector"),"p");
                }
            }
            mean_mg->Add(mean_plot);

            if (algo == 'j')
            {
                TGraph* trunc_mean_plot = new TGraphErrors(n_ch_max, plot_channel_list, trunc_mean_list[1][run], error_x, trunc_mean_err_list[1][run]);
                trunc_mean_plot->SetMarkerStyle(22);
                trunc_mean_plot->SetMarkerColor(3 * run + 80);
                trunc_mean_mg->Add(trunc_mean_plot);

                //for Poisson
                /*TGraph* exp_plot = new TGraphErrors(n_ch_max, plot_channel_list, poisson_exp_list[1][run], error_x, poisson_exp_err_list[1][run]);
                exp_plot->SetMarkerStyle(22);
                exp_plot->SetMarkerSize(1);
                exp_plot->SetMarkerColor(run + 41);
                exp_mg->Add(exp_plot);*/

                //for Langau
                TGraph* mpv_plot = new TGraphErrors(n_ch_max, plot_channel_list, langau_mpv_list[1][run], error_x, langau_mpv_err_list[1][run]);
                mpv_plot->SetMarkerStyle(22);
                mpv_plot->SetMarkerColor(3 * run + 80);
                mpv_mg->Add(mpv_plot);
            }

            if (algo == 'c')
            {
                TGraph* max_plot = new TGraphErrors(n_ch_max, plot_channel_list, poisson_exp_list[1][run], error_x, poisson_exp_err_list[1][run]);
                max_plot->SetMarkerStyle(22);
                max_plot->SetMarkerColor(3 * run + 80);
                //leg->Draw("same");
                max_mg->Add(max_plot);
            }
        }
    }


    //drawing the mean MPV values
    for (int i = 0; i < 8; i++) plot_channel_list[i] = i;
    for (int measurement = 0; measurement < n_m; measurement++)
    {
        TGraph* mean_mpv_plot = new TGraphErrors(n_ch_max, plot_channel_list, mean_mpv_values[measurement], error_x, mean_mpv_err_values[measurement]);
        //for plots without reflector
        if (measurement == 0)
        {
            mean_mpv_plot->SetMarkerSize(1.5);
            mean_mpv_plot->SetMarkerStyle(28);
            mean_mpv_plot->SetMarkerColor(4);
            //add the mean MPV to the legend
            if (with_black_foil)
                leg_mpv->AddEntry(mean_mpv_plot, Form("Mean N_{PE} without black foil"),"p");
            else
                leg_mpv->AddEntry(mean_mpv_plot, Form("Mean N_{PE} without reflector"),"p");
        }
        if (measurement == 1)
        {
            mean_mpv_plot->SetMarkerStyle(29);
            mean_mpv_plot->SetMarkerColor(2);
            if (with_black_foil)
                leg_mpv->AddEntry(mean_mpv_plot, Form("Mean N_{PE} with black foil"),"p");
            else
                leg_mpv->AddEntry(mean_mpv_plot, Form("Mean N_{PE} with reflector"),"p");
        }

        //draw the mean MPV plot into the MPV multigraph
        mpv_mg->Add(mean_mpv_plot);
    }

    //mean
    masterCanvas->cd(1);
    mean_mg->Draw("APE");
    leg_mean->Draw("same");

    if (algo == 'j')
    {
        //truncated mean
        masterCanvas->cd(2);
        trunc_mean_mg->Draw("APE");
        leg_trunc->Draw("same");

        //Poisson exp
        //masterCanvas->cd(3);
        //exp_mg->Draw("APE");

        //Langau MPV
        masterCanvas->cd(3);
        mpv_mg->Draw("APE");
        leg_mpv->Draw("same");
        //string langau_out_path = file_path + "Channel parameters " + algo_name + " " + pcb_name + " langau.pdf";
        //masterCanvas->cd(3)->Print(langau_out_path.c_str());
    }

    if (algo == 'c')
    {
        //maximum value of the fit
        masterCanvas->cd(2);
        max_mg->Draw("APE");
    }

    

    string out_path = file_path + "Channel parameters " + algo_name + " " + pcb_name;
    if (with_black_foil) out_path += " with black foil";
    out_path += ".pdf";
    masterCanvas->Print(out_path.c_str());

    string langau_out_path = file_path + "Channel parameters " + algo_name + " " + pcb_name + " langau";
    if (with_black_foil) langau_out_path += " with black foil";
    langau_out_path += ".pdf";
    masterCanvas->cd(3)->Print(langau_out_path.c_str());

    return 0;
}

