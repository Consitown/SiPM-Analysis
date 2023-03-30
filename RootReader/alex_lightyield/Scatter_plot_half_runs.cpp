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
#include <TLatex.h>
#include <TAttMarker.h>
#include <TMarker.h>
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
    double run_number_list_part[12];
    double dist_max_list_part[12];
    double dist_mean_list_part[12];
    double mpv_list_part[12];
    double gauss_mean_list_part[12];
    double temperature_list_part[12];

    //for the errors
    double error_x[12] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double dist_max_error_y[12] {2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67};
    double gauss_max_error_y_part[12];


    //Variable for temperature correction
    double correction_factor;

    //variable to read single lines from the file
    string line;

    //variable to store numbers as strings first
    string help_str;

    //variable to convert from help_str to number
    double help_num;

    //for the marker colour later
    div_t divresult;
    Double_t ax[12],ay[12];

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
    /*
    cout << "\n \n \nWhat I have read as Values.txt:\n" << endl;
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
    for (int i = 0; i<12; i++)
    {
        run_number_list_part[i] = run_number_list[i];
        dist_max_list_part[i] = dist_max_list[i];
        dist_mean_list_part[i] = dist_mean_list[i];
        mpv_list_part[i] = mpv_list[i];
        gauss_mean_list_part[i] = gauss_mean_list[i];
        temperature_list_part[i] = temperature_list[i];
        gauss_max_error_y_part[i] = gauss_max_error_y[i];
    } 
 


    //Distribution maximum
    TGraph* dist_max_plot_c = new TGraphErrors(12, temperature_list_part, dist_max_list_part, error_x, dist_max_error_y);
    dist_max_plot_c->SetTitle("Distribution maximum for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(1);
    dist_max_plot_c->SetMarkerColor(0);
    dist_max_plot_c->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_max_plot_c->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }

    //Distribution mean
    /*TGraph* dist_mean_plot_c = new TGraphErrors(12, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_c->SetTitle("Distribution mean for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    dist_mean_plot_c->SetMarkerColor(0);
    dist_mean_plot_c->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_mean_plot_c->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_c = new TGraphErrors(12, temperature_list_part, mpv_list_part, error_x, gauss_max_error_y_part);
    mpv_plot_c->SetTitle("Gaussian maximum for PCB C; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(5);
    mpv_plot_c->SetMarkerColor(0);
    mpv_plot_c->Draw("APE");

    for (int i=0; i<12; i++)
    {
        mpv_plot_c->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }


    /*
    ----------------PCB J-------------------------------
    */
    
    for (int i = 0; i<12; i++)
    {
        run_number_list_part[i] = run_number_list[i + 12];
        dist_max_list_part[i] = dist_max_list[i + 12];
        dist_mean_list_part[i] = dist_mean_list[i + 12];
        mpv_list_part[i] = mpv_list[i + 12];
        gauss_mean_list_part[i] = gauss_mean_list[i + 12];
        temperature_list_part[i] = temperature_list[i + 12];
        gauss_max_error_y_part[i] = gauss_max_error_y[i + 12];
    } 

    
    //distribution maximum
    TGraph* dist_max_plot_j = new TGraphErrors(12, temperature_list_part, dist_max_list_part, error_x, dist_max_error_y);
    dist_max_plot_j->SetTitle("Distribution maximum for Juelich PCB (J); Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(2);
    dist_max_plot_j->SetMarkerColor(0);
    dist_max_plot_j->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_max_plot_j->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }

    //Distribution mean
    /*TGraph* dist_mean_plot_j = new TGraphErrors(12, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_j->SetTitle("Distribution mean for Juelich PCB (J); Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    dist_mean_plot_j->SetMarkerColor(0);
    dist_mean_plot_j->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_mean_plot_j->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_j = new TGraphErrors(12, temperature_list_part, mpv_list_part, error_x, gauss_max_error_y_part);
    mpv_plot_j->SetTitle("Gaussian maximum for Juelich PCB (J); Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(6);
    mpv_plot_j->SetMarkerColor(0);
    mpv_plot_j->Draw("APE");

    for (int i=0; i<12; i++)
    {
        mpv_plot_j->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }


    /*
    ----------------PCB D-------------------------------
    */

    for (int i = 0; i<12; i++)
    {
        run_number_list_part[i] = run_number_list[i + 24];
        dist_max_list_part[i] = dist_max_list[i + 24];
        dist_mean_list_part[i] = dist_mean_list[i + 24];
        mpv_list_part[i] = mpv_list[i + 24];
        gauss_mean_list_part[i] = gauss_mean_list[i + 24];
        temperature_list_part[i] = temperature_list[i + 24];
        gauss_max_error_y_part[i] = gauss_max_error_y[i + 24];
    } 

    
    //distribution maximum
    TGraph* dist_max_plot_d = new TGraphErrors(12, temperature_list_part, dist_max_list_part, error_x, dist_max_error_y);
    dist_max_plot_d->SetTitle("Distribution maximum for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(3);
    dist_max_plot_d->SetMarkerColor(0);
    dist_max_plot_d->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_max_plot_d->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }

    //Distribution mean
    /*TGraph* dist_mean_plot_d = new TGraphErrors(12, temperature_list_part, dist_mean_list_part);
    dist_mean_plot_d->SetTitle("Distribution mean for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    dist_mean_plot_d->SetMarkerColor(0);
    dist_mean_plot_d->Draw("APE");

    for (int i=0; i<12; i++)
    {
        dist_mean_plot_d->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }*/

    //Gaussian maximum/MPV
    TGraph* mpv_plot_d = new TGraphErrors(12, temperature_list_part, mpv_list_part, error_x, gauss_max_error_y_part);
    mpv_plot_d->SetTitle("Gaussian maximum for PCB D; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(7);
    mpv_plot_d->SetMarkerColor(0);
    mpv_plot_d->Draw("APE");

    for (int i=0; i<12; i++)
    {
        mpv_plot_d->GetPoint(i,ax[i],ay[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax[i],ay[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }



    /*
    ----------------PCB D but without new optical gel---------------------
    */

    //we only have data points, as the temperature measurement stopped
    double run_number_list_part2[8];
    double dist_max_list_part2[8];
    double dist_mean_list_part2[8];
    double mpv_list_part2[8];
    double gauss_mean_list_part2[8];
    double temperature_list_part2[8];

    double error_x2[8] {0, 0, 0, 0, 0, 0, 0, 0};
    double dist_max_error_y2[8] {2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67, 2.67};
    double gauss_max_error_y_part2[8];

    for (int i = 0; i<8; i++)
    {
        run_number_list_part2[i] = run_number_list[i + 36];
        dist_max_list_part2[i] = dist_max_list[i + 36];
        dist_mean_list_part2[i] = dist_mean_list[i + 36];
        mpv_list_part2[i] = mpv_list[i + 36];
        gauss_mean_list_part2[i] = gauss_mean_list[i + 36];
        temperature_list_part2[i] = temperature_list[i + 36];
        gauss_max_error_y_part2[i] = gauss_max_error_y[i + 36];

        //for debugging
        //cout << run_number_list_part2[i] << "\t" << dist_max_list_part2[i] << "\t" << dist_mean_list_part2[i] << "\t" << mpv_list_part2[i] << "\t" << gauss_mean_list_part2[i] << "\t" << temperature_list_part2[i] << endl;
    } 

    Double_t ax2[8],ay2[8];
    //distribution maximum
    TGraph* dist_max_plot_d2 = new TGraphErrors(8, temperature_list_part2, dist_max_list_part2, error_x2, dist_max_error_y2);
    dist_max_plot_d2->SetTitle("Distribution maximum for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(4);
    dist_max_plot_d2->SetMarkerColor(0);
    dist_max_plot_d2->Draw("APE");

    for (int i=0; i<8; i++)
    {
        dist_max_plot_d2->GetPoint(i,ax2[i],ay2[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax2[i],ay2[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }

    //Distribution mean
    /*TGraph* dist_mean_plot_d2 = new TGraphErrors(8, temperature_list_part2, dist_mean_list_part2);
    dist_mean_plot_d2->SetTitle("Distribution mean for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    dist_mean_plot_d2->SetMarkerColor(0);
    dist_mean_plot_d2->Draw("APE");

    for (int i=0; i<8; i++)
    {
        dist_mean_plot_d2->GetPoint(i,ax2[i],ay2[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        


        divresult = div (i,2);

        TMarker* m = new TMarker(ax2[i],ay2[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }*/



    //Gaussian maximum/MPV
    TGraph* mpv_plot_d2 = new TGraphErrors(8, temperature_list_part2, mpv_list_part2, error_x2, gauss_max_error_y_part2);
    mpv_plot_d2->SetTitle("Gaussian maximum for PCB D without new optical gel; Temperature in ^{#circ}C; Pulse height [mV x ns]");
    masterCanvas->cd(8);
    
    mpv_plot_d2->SetMarkerColor(0);
    mpv_plot_d2->Draw("APE");

    //for each two points (=one run), choose one colour
    
    for (int i=0; i<8; i++)
    {
        mpv_plot_d2->GetPoint(i,ax2[i],ay2[i]); 

        /*
        //just for debugging
        cout<<i<<"th element of X array: "<<ax[i]<<endl;
        cout<<i<<"th element of Y array: "<<ay[i]<<endl;
        */


        divresult = div (i,2);

        TMarker* m = new TMarker(ax2[i],ay2[i],7);
        m->SetMarkerColor(divresult.quot + 1); //0 would be white
        m->Draw("SAME");
    }






    //print everything to a pdf

    masterCanvas->Print("Correlation plot half runs.pdf");

}

