// reads finished C-0 run file and gives light yield distribution for SiPMs and in the x-y-plane

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

#include "/mnt/d/Work_SHK_Bachelor/RootReader/src/meanAngleFuncs.h"


using namespace std;

int main(int argc, char *argv[]) {

    //defining the directories
    string runName = "22_cosmics_pcbj_everyothergroup_vb42_PSmid";
    string in_file_path = "/mnt/d/Work_SHK_Bachelor/RootReader/runs/";
    string in_filename = in_file_path + "Sicherung_0/" + runName + ".root"; //check the naming of the files -- really important!!!

    string couplingCorrection_path = "/mnt/d/Work_SHK_Bachelor/RootReader/RootAnalysis/couplingCorrection/";

    string out_filename = couplingCorrection_path + "lyc0.txt";

    //opening the old root file
    cout << "Opening old root file at " << in_filename << endl;
    TFile oldfile(in_filename.c_str());

    if (oldfile.IsZombie())
    {
        cout << "Problem with file " << in_filename << "; check if file path is correct!" << endl;
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

    cout << "n_events = " << n_events << endl;

    float phi_list[n_ch] = {0, 315, 270, 225, 45, 90, 135, 180}; //angles corresponding to the channels in degree

    float ly_event_channel[n_ch][n_events]; //lightyield of each event, so the integral value of one event that is saved in the root file
    float ly_channel_av[n_ch];
    
    //reading the individual integral values from the root files into the ly_event list

    for (int channel_counter = 0; channel_counter < n_ch; channel_counter++)
    {
        for (int event_counter = 0; event_counter < n_events; event_counter++)
        {
            //for light yield
            TString integral_name;
            integral_name.Form("integral_hist_%d", channel_counter);
            oldtree->SetBranchAddress(integral_name, &ly_event_channel[channel_counter][event_counter]);
            oldtree->GetEntry(event_counter);

            ly_channel_av[channel_counter] += ly_event_channel[channel_counter][event_counter];
        }
    }

    for (int channel = 0; channel < 8; channel++) ly_channel_av[channel] /= n_events; //dividing by n_events to get the average light yield

    cout << "Finished opening and reading." << endl;

    cout << "\nOpening " << out_filename << endl;

    ofstream ly_file (out_filename.c_str());
    for (int channel = 0; channel < 8; channel++)   ly_file << showpoint << setprecision(4) << ly_channel_av[channel] << endl;

    ly_file.close();
    cout << "Created .txt file." << endl;

}