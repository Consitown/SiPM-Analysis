// reads finished Rootfile from runs folder. Normalizes the x_i and y_i and recreates the phi_ew values
// corrects using the values given in /RootAnalysis/couplingCorrection/couplingCorrection.txt (old version)
// corrects by weighing with the normalized light yield of the C-0 position (new version)

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

#include "/home/alex/Dokumente/Studium/RootReader/src/meanAngleFuncs.h"


using namespace std;

int main(int argc, char *argv[]) {
    bool correction = true; //whether or not to correct for non-linearity, does not affect the normalization

    //defining the directories

    string pos = "2";
    
    //string runName = "1000_12to18_day";
    //string runName = "1222_02to12_day";
    string runName = "1333_1307to0408_day";
    string in_file_path = "/home/alex/Dokumente/Studium/RootReader/runs/";
    string in_filename = in_file_path + "Sicherung_" + pos + "/" + runName + ".root";

    string couplingCorrection_path = "/home/alex/Dokumente/Studium/RootReader/RootAnalysis/couplingCorrection/";

    string lyc0_filename = couplingCorrection_path + "lyc0.txt";
    string out_filename = couplingCorrection_path + runName + ".root";

    //for the x-&y-plots
    string plot_path = "/home/alex/Dokumente/Studium/RootReader/RootAnalysis/integralAnalysis/" + runName + "/";
    string xy_plot_filename = plot_path + "x-y-plot.pdf";
    string new_xy_plot_filename = plot_path + "new_x-y-plot.pdf";
    string ly_plot_filename = plot_path + "normalized_ly_plot.pdf";
    string new_ly_plot_filename = plot_path + "new_normalized_ly_plot.pdf";

    //opening the old root file
    cout << "Opening old root file at " << in_filename << endl;
    TFile oldfile(in_filename.c_str());

    if (oldfile.IsZombie())
    {
        cout << "Problem with file " << in_filename << "; check if file path is correct!" << endl;
        exit(-1);
    }

    //opening the new root file
    cout << "Opening new root file at " << out_filename << endl;
    TFile newfile(out_filename.c_str(), "RECREATE");

    if (newfile.IsZombie())
    {
        cout << "Problem with file " << out_filename << "; check if file path is correct!" << endl;
        exit(-1);
    }


    //opening the old tree
    TTree* oldtree = new TTree;
    oldfile.GetObject("T", oldtree);
    oldtree->GetEntry(1);

    

    //initializing the variables
    int n_ch = 8; //number of channels
    int channel_list[8];
    for (int i=0; i<n_ch; i++) channel_list[i] = i;
    int n_events; //number of events
    n_events = oldtree->GetEntries();

    float cartXarray[n_ch][n_events], cartYarray[n_ch][n_events];
    float x_i_bar[8];
    float y_i_bar[8];

    float phi_ew;
    float phiStd;
    Float_t threePhiEw[3];
    float sumCartX_list[n_events], sumCartY_list[n_events];
    float sumCartX, sumCartY; //uses the sumCartX values from the root file
    float sumCartX_new, sumCartY_new; //calculates the sumCartX values using the CartX (x_i) values from the root file; should not make a difference to sumCartX in theory!
    float sigmaX, sigmaY;

    float X_bar, Y_bar;

    cout << "n_events = " << n_events << endl;

    float x_i[n_ch][n_events], y_i[n_ch][n_events]; //x_i and y_i for each channel

    float ly_event_channel[n_ch][n_events]; //lightyield of each event, so the integral value of one event that is saved in the root file
    float ly_event_res[n_events] = {0}; //lightyield from each event, summed over all channels
    float ly_channel_av[n_ch] = {0}; //av lightyield in each channel (so resulting light yield divided by n_events)

    float phi_ew_compare[n_events];

    string line; //for .txt readout
    string help_string;
    int line_counter;
    double help_num;

    //reading the individual integral values from the root files into the ly_event list
    //reading the individual x_i and y_i values from the root file into x_i and y_i lists
    for (int event_counter = 0; event_counter < n_events; event_counter++)
    {
        for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
        {
            //for light yield
            TString integral_name;
            integral_name.Form("integral_hist_%d", channel_counter);
            oldtree->SetBranchAddress(integral_name, &ly_event_channel[channel_counter][event_counter]);
            oldtree->GetEntry(event_counter);

            ly_event_res[event_counter] += ly_event_channel[channel_counter][event_counter];

            //for x_i and y_i values
            TString x_name, y_name;
            x_name.Form("cartX%d", channel_counter);
            oldtree->SetBranchAddress(x_name, &x_i[channel_counter][event_counter]);
            oldtree->GetEntry(event_counter);

            y_name.Form("cartY%d", channel_counter);
            oldtree->SetBranchAddress(y_name, &y_i[channel_counter][event_counter]);
            oldtree->GetEntry(event_counter);

        }
        
        oldtree->SetBranchAddress("sumCartX", &sumCartX_list[event_counter]);
        oldtree->GetEntry(event_counter);

        oldtree->SetBranchAddress("sumCartY", &sumCartY_list[event_counter]);
        oldtree->GetEntry(event_counter);

        oldtree->SetBranchAddress("Phi_ew_all_ch", &phi_ew_compare[event_counter]);
        oldtree->GetEntry(event_counter);
    }
    int centerChannel = 0;
    oldtree->SetBranchAddress("CenterChannelPhiEw", &centerChannel);
    oldtree->GetEntry(1);

    cout << "centerChannel = " << centerChannel << endl;

    cout << "Finished opening and reading." << endl;
    //just debugging output
    /*for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < 4; event_counter++)
        {
            cout << std::setprecision(4) << "x_i[" << channel_counter <<"][" << event_counter << "]" << x_i[channel_counter][event_counter] << endl; 
        }
    }*/

    //reading x_i_bar and y_i_bar values from couplingCorrection.txt 
    if (correction == true)
    {
        string values_path = couplingCorrection_path + "couplingCorrection.txt";
        
        line_counter = 0;

        cout << "Opening " << values_path << endl;
        
        ifstream couplingCorrection_file (values_path.c_str());

        if (couplingCorrection_file.is_open())
        {
            while (getline (couplingCorrection_file, line))
            {
            help_string = "";
            for (int char_counter = 2; char_counter <= 8; char_counter++)
            {
                help_string += line[char_counter];
            }
            help_num = atof(help_string.c_str());

            x_i_bar[line_counter] = help_num;

            help_string = "";
            for (int char_counter = 11; char_counter <= 17; char_counter++)
            {
                help_string += line[char_counter];
            }
            help_num = atof(help_string.c_str());

            y_i_bar[line_counter] = help_num;

            line_counter++;
            }
        }
        else cout << "Error while opening " << couplingCorrection_path << endl;

        X_bar = 0;
        Y_bar = 0;
        for (int i=0; i<n_ch; i++)
        {
            X_bar += x_i_bar[i];
            Y_bar += y_i_bar[i];
        }

        couplingCorrection_file.close();
    }


    //reading the average J_i^(C0) values from the .txt file
    if (correction == true)
    {
        cout << "Opening " << lyc0_filename << endl;
        ifstream lyc0_file(lyc0_filename.c_str());

        line_counter = 0;

        if (lyc0_file.is_open())
        {
            while (getline (lyc0_file, line))
            {
            help_string = "";
            for (int char_counter = 0; char_counter <= 4; char_counter++)
            {
                help_string += line[char_counter];
            }
            help_num = atof(help_string.c_str());

            ly_channel_av[line_counter] = help_num;
            line_counter++;
            }
        }
        else cout << "Error while opening " << couplingCorrection_path << endl;

        lyc0_file.close();
        cout << "Read lyc0 values from file." << endl;
    }

    //for (int i=0;i<8;i++) cout << ly_channel_av[i] <<endl;



    //cloning the old tree into a new tree, but without phi, x and y value
    //deactivating branches that shouldn't be copied
    TString x_branch_name, y_branch_name;
    for (int i = 0; i < n_ch; i++)
    {
        x_branch_name.Form("cartX%d", i);
        y_branch_name.Form("cartY%d", i);

        oldtree->SetBranchStatus(x_branch_name, 0);
        oldtree->SetBranchStatus(y_branch_name, 0);
    }

    oldtree->SetBranchStatus("Phi_ew_all_ch", 0);
    oldtree->SetBranchStatus("Std_Phi_ew_all", 0);
    oldtree->SetBranchStatus("Phi_ew_shifted", 0);

    oldtree->SetBranchStatus("sumCartX", 0);
    oldtree->SetBranchStatus("sumCartY", 0);

    oldtree->SetBranchStatus("sigmaSumCartX", 0);
    oldtree->SetBranchStatus("sigmaSumCartY", 0);

    TTree* newtree = oldtree->CloneTree();

    //normalizing
    for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < n_events; event_counter++)
        {        
            x_i[channel_counter][event_counter] /= ly_event_res[event_counter]; //normalise by light yield of the event in all channels
            y_i[channel_counter][event_counter] /= ly_event_res[event_counter]; //normalise by light yield of the event in all channels        
        }
    }

    //normalizing for C-0 as proposed by Heiko
    if (correction == true)
    {  
        for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
        {
            for (int event_counter = 0; event_counter < n_events; event_counter++)
            {        
                x_i[channel_counter][event_counter] /= ly_channel_av[channel_counter]; //normalise by light yield of C-0 measurement
                y_i[channel_counter][event_counter] /= ly_channel_av[channel_counter]; //normalise by light yield of C-0 measurement
            }

        }
    }

    
    //saving the new data to the newtree
   

    for (int i=0; i<8; i++) 
    {
        newtree->Branch(Form("cartX%d", i), &cartXarray[i], Form("cartX%d/F", i));
        newtree->Branch(Form("cartY%d", i), &cartYarray[i], Form("cartY%d/F", i));
    }

    newtree->Branch("sumCartX", &sumCartX_new, "sumCartX/F");
    newtree->Branch("sumCartY", &sumCartY_new, "sumCartY/F");
    newtree->Branch("sigmaSumCartY", &sigmaY, "sigmaSumCartY/F");
    newtree->Branch("sigmaSumCartX", &sigmaX, "sigmaSumCartX/F");

    newtree->Branch("Phi_ew_all_ch", &phi_ew, "Phi_ew_all_ch/F");
    newtree->Branch("Std_Phi_ew_all", &phiStd, "Std_Phi_ew_all/F");
    newtree->Branch("Phi_ew_shifted", threePhiEw, "Phi_ew_shifted[3]/F");


    //normalizing the correction values
    for (int channel_counter = 0; channel_counter < 8; channel_counter++)
    {
        x_i_bar[channel_counter] /= ly_channel_av[channel_counter];
        y_i_bar[channel_counter] /= ly_channel_av[channel_counter];
    }

    for (int event_counter = 0; event_counter < n_events; event_counter++)
    {
        sumCartX_new = 0;
        sumCartY_new = 0;
        //save cartX and cartY
        for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
        { 
            if (correction == true) 
            {             
                cartXarray[channel_counter][event_counter] = x_i[channel_counter][event_counter]; // - x_i_bar[channel_counter]; (old version)
                cartYarray[channel_counter][event_counter] = y_i[channel_counter][event_counter]; // - y_i_bar[channel_counter]; (old version)

                
            }
            if (correction == false)
            {
                cartXarray[channel_counter][event_counter] = x_i[channel_counter][event_counter];
                cartYarray[channel_counter][event_counter] = y_i[channel_counter][event_counter];
            }

            sumCartX_new += cartXarray[channel_counter][event_counter];
            sumCartY_new += cartYarray[channel_counter][event_counter];

            /*sumCartX += 1/8. * cartXarray[channel_counter];
            sumCartY += 1/8. * cartYarray[channel_counter];*/


            newtree->GetBranch(Form("cartX%d", channel_counter))->Fill();
            newtree->GetBranch(Form("cartY%d", channel_counter))->Fill();
        }

        //if (sumCartY_list[event_counter] == sumCartY_list[event_counter + 1]) cout << "Komisches Zeug bei Event " << event_counter << ": \t" << sumCartY_list[event_counter] << endl;

        sumCartX = sumCartX_list[event_counter];
        sumCartY = sumCartY_list[event_counter]; //using the values directly from the run file

        //normalizing
        sumCartX /= ly_event_res[event_counter];
        sumCartY /= ly_event_res[event_counter];

        //correcting
        /*if (correction == true)
        {
            sumCartX -= X_bar;
            sumCartY -= Y_bar;
        }*/
        
        //cout << sumCartX_new << "\t" << sumCartX << endl;

        //cout << std::setprecision(4) << sumCartX << "\t" << sumCartY << "\t" << sumCartY / sumCartX << endl;

        //save sigmaX and sigmaY
        float sigmaX_raw;
        float sigmaY_raw;
        for (int j=0; j<8; j++)
        {
            sigmaX_raw += (sumCartX_new - cartXarray[j][event_counter])*(sumCartX_new - cartXarray[j][event_counter]);
            sigmaX_raw += (sumCartY_new - cartYarray[j][event_counter])*(sumCartY_new - cartYarray[j][event_counter]);
        }
        sigmaX = TMath::Sqrt(sigmaX_raw * 1/8.0);
        sigmaY = TMath::Sqrt(sigmaY_raw * 1/8.0);
        
        newtree->GetBranch("sumCartX")->Fill();
        newtree->GetBranch("sumCartY")->Fill();

        newtree->GetBranch("sigmaSumCartX")->Fill();
        newtree->GetBranch("sigmaSumCartY")->Fill();

        //save phi_ew
        //cout << std::setprecision(4) << sumCartX << "\t" << sumCartY << "\t" << sumCartY / sumCartX << endl;
        phi_ew = cartesianToPolar(sumCartX_new, sumCartY_new);
        //for the threefold shifted plot
        Int_t angles[8] = {0, 315, 270, 225, 45,  90, 135, 180};

        phi_ew = translateAngle(phi_ew, angles[centerChannel]);
        /*if (abs(phi_ew_compare[event_counter] - phi_ew) > 0.000000001)
        {
            cout << phi_ew_compare[event_counter] << "\t" << phi_ew << " at " << event_counter << endl;
        }*/

        //if (event_counter >= 1759) cout << "phi_ew_compare[" << event_counter << "] = " << phi_ew_compare[event_counter] << endl;
        //if (event_counter >= 1759) cout << "phi_ew[" << event_counter << "] = " << phi_ew << endl;
        newtree->GetBranch("Phi_ew_all_ch")->Fill();

        // calculate the standard deviation of phi_ew for each event
        float individualPhi[8];

        float cartX_help[8], cartY_help[8];
        for (int i = 0; i < 8; i++)
        {
            cartX_help[i] = cartXarray[i][event_counter];
            cartY_help[i] = cartYarray[i][event_counter];
        } 
        cartesianToPolar(8, cartX_help, cartY_help, individualPhi);
        float sigma_raw = 0;
        float difference = 0;
        for (int m=0; m<8; m++)
        {
            difference = phi_ew - individualPhi[m];
            if (difference > 180) difference -= 180;
            else if (difference < -180) difference += 180;
            sigma_raw += difference * difference; // XXX
            // cout << sigma_raw << endl;
        }
        phiStd = sqrt(1/8.0 * sigma_raw);
        newtree->GetBranch("Std_Phi_ew_all")->Fill();

        
        threePhiEw[0] = phi_ew - 360.0;
        threePhiEw[1] = phi_ew;
        threePhiEw[2] = phi_ew + 360.0;
        newtree->GetBranch("Phi_ew_shifted")->Fill();
    }


    /*for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < 4; event_counter++)
        {
            cout << "x_i[" << channel_counter <<"][" << event_counter << "]" << x_i[channel_counter][event_counter] << endl; 
        }
    }*/

    gErrorIgnoreLevel = kWarning;
    
    newtree->AutoSave();
    //newtree->Print();

    oldfile.Close();
    newfile.Write();
    newfile.Close();

    cout << "Normalization done!" << endl;

    /*
    _______________________________GRAPHIC EXPORT xy plot______________________________________
    */
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

    //plotting parameters
    float bin_width;
    int n_bins = 100;


    //the canvas in which all xy plots are saved
    TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 3000);
    masterCanvas->Divide(2, 8);
    masterCanvas->SetGrid();

    TH1F *x_plots_cor[n_ch], *x_plots_old[n_ch];
    TH1F *y_plots_cor[n_ch], *y_plots_old[n_ch];

    for (int channel_counter = 0; channel_counter < 8; channel_counter++)
    {
        //x plots
        //set the bin width so that the distributions look the same
        if (abs(x_i_bar[channel_counter]) > 0) bin_width = abs(x_i_bar[channel_counter]) / 10;
        else bin_width = 0.00001;

        
        //corrected
        masterCanvas->cd(2 * channel_counter + 1);
        x_plots_cor[channel_counter] = new TH1F(Form("x_cor%d", channel_counter), Form("x in channel %d", channel_counter), n_bins, -100*bin_width, 100*bin_width);
        for (int i = 0; i<n_events; i++) x_plots_cor[channel_counter]->Fill(cartXarray[channel_counter][i]);
        x_plots_cor[channel_counter]->Draw();

        //uncorrected
        x_plots_old[channel_counter] = new TH1F(Form("x_old%d", channel_counter), Form("x in channel %d", channel_counter), n_bins, -100*bin_width, 100*bin_width);
        for (int i = 0; i<n_events; i++) x_plots_old[channel_counter]->Fill(x_i[channel_counter][i]);
        x_plots_old[channel_counter]->SetLineColorAlpha(kBlack, 0.6);
        x_plots_old[channel_counter]->Draw("same");

        //legend
        TLegend* xLeg = new TLegend(0.6, 0.65, 0.9, 0.9);
        xLeg->SetBorderSize(1);
        xLeg->AddEntry((TObject*)0, Form("Entries = %d", (int) x_plots_cor[channel_counter]->GetEntries()));
        xLeg->AddEntry(Form("x_old%d", channel_counter), Form("Original x values"));
        xLeg->AddEntry(Form("x_cor%d", channel_counter), Form("Corrected x values"));
        gStyle->SetLegendTextSize(0.04);
        xLeg->Draw();


        //y plots

        if (abs(y_i_bar[channel_counter]) > 0) bin_width = abs(y_i_bar[channel_counter]) / 10;
        else bin_width = 0.01;
        //corrected
        masterCanvas->cd(2 * channel_counter + 2);
        y_plots_cor[channel_counter] = new TH1F(Form("y_cor%d", channel_counter), Form("y in channel %d", channel_counter), n_bins, -100*bin_width, 100*bin_width);
        for (int i = 0; i<n_events; i++) y_plots_cor[channel_counter]->Fill(cartYarray[channel_counter][i]);
        y_plots_cor[channel_counter]->Draw();

        //uncorrected
        y_plots_old[channel_counter] = new TH1F(Form("y_old%d", channel_counter), Form("y in channel %d", channel_counter), n_bins, -100*bin_width, 100*bin_width);
        for (int i = 0; i<n_events; i++) y_plots_old[channel_counter]->Fill(y_i[channel_counter][i]);
        y_plots_old[channel_counter]->SetLineColorAlpha(kBlack, 0.6);
        y_plots_old[channel_counter]->Draw("same");

        //legend
        TLegend* yLeg = new TLegend(0.6, 0.65, 0.9, 0.9);
        yLeg->SetBorderSize(1);
        yLeg->AddEntry((TObject*)0, Form("Entries = %d", (int) y_plots_cor[channel_counter]->GetEntries()));
        yLeg->AddEntry(Form("y_old%d", channel_counter), Form("Original y values"));
        yLeg->AddEntry(Form("y_cor%d", channel_counter), Form("Corrected y values"));
        gStyle->SetLegendTextSize(0.04);
        yLeg->Draw();

    }

    TPaveLabel title(0.1, 0.9, 0.9, 0.9, Form("Normalized light yields"));
    title.SetTextSize(.7);
    title.SetBorderSize(0);
    title.Draw("same");
    gPad->SetTitle("Hello");
    gPad->Update();
    masterCanvas->Print(xy_plot_filename.c_str());

    /*
    _______________________________GRAPHIC EXPORT light yield plot______________________________________
    */
    //the canvas in which all light yield plots are saved
    TCanvas *lyCanvas = new TCanvas("LY", "Normalized light yields", 1000, 3000);
    lyCanvas->Divide(2, 8);
    lyCanvas->SetGrid();
    lyCanvas->SetTitle("Normalized light yields");

    TH1F *ly_plots[n_ch];

    for (int channel_counter = 0; channel_counter < 8; channel_counter++)
    {
        lyCanvas->cd(channel_counter + 1);

        ly_plots[channel_counter] = new TH1F(Form("ly_%d", channel_counter), Form("channel %d", channel_counter), 100, 0, 0.5);
        for (int i = 0; i<n_events; i++) ly_plots[channel_counter]->Fill(ly_event_channel[channel_counter][i] / ly_event_res[channel_counter]);
        ly_plots[channel_counter]->Draw();
        
    }
    TPaveText ly_title(0.1, 0.9, 0.9, 0.9, Form("Normalized light yields"));
    ly_title.SetTextSize(.7);
    ly_title.SetBorderSize(0);
    ly_title.Draw("same");
    lyCanvas->Print(ly_plot_filename.c_str());

    /*
    _______________________________GRAPHIC EXPORT new light yield plot______________________________________
    */
    //the canvas in which all light yield plots are saved
    TCanvas *new_lyCanvas = new TCanvas("newLY", "New Normalized light yields", 1000, 3000);
    new_lyCanvas->Divide(2, 8);
    new_lyCanvas->SetGrid();
    new_lyCanvas->SetTitle("New Normalized light yields");

    TH1F *new_ly_plots[n_ch];

    for (int channel_counter = 0; channel_counter < 8; channel_counter++)
    {
        new_lyCanvas->cd(channel_counter + 1);

        new_ly_plots[channel_counter] = new TH1F(Form("new_ly_%d", channel_counter), Form("channel %d", channel_counter), 100, 0, 0.0005);
        for (int i = 0; i<n_events; i++) new_ly_plots[channel_counter]->Fill(ly_event_channel[channel_counter][i] / ly_event_res[channel_counter] / ly_channel_av[channel_counter]);
        new_ly_plots[channel_counter]->Draw();
        
    }
    TPaveLabel new_ly_title(0.1, 0.9, 0.9, 0.9, Form("Normalized light yields"));
    new_ly_title.SetTextSize(.7);
    new_ly_title.SetBorderSize(0);
    new_ly_title.Draw("same");
    new_lyCanvas->Print(new_ly_plot_filename.c_str());


}