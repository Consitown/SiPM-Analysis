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
#include <TError.h> // root verbosity level
#include <TApplication.h>

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
#include <dirent.h>
#include <TH1F.h>
#include <THStack.h>
#include <algorithm>
#include <vector>

using namespace std;

void strip_ext(char *fname)
{
    char *end = fname + strlen(fname);

    while (end > fname && *end != '.')
    {
        --end;
    }

    if (end > fname)
    {
        *end = '\0';
    }
}

vector<string> split(const string &str, const string &delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos)
            pos = str.length();
        string token = str.substr(prev, pos - prev);
        if (!token.empty())
            tokens.push_back(token);
        prev = pos + delim.length();
    } while (pos < str.length() && prev < str.length());
    return tokens;
}
void extractIntegerWords(string str)
{
    stringstream ss;

    /* Storing the whole string into string stream */
    ss << str;

    /* Running loop till the end of the stream */
    string temp;
    int found;
    while (!ss.eof())
    {

        /* extracting word by word from stream */
        ss >> temp;

        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            cout << found << " ";

        /* To save from space at the end of string */
        temp = "";
    }
}

int main(int argc, char *argv[])
{

	TApplication *myapp = new TApplication("myapp", 0, 0);

    struct dirent *entry = nullptr;
    DIR *dp = nullptr;
    THStack *hs = new THStack("hs", "");
    TCanvas *effCanvas = new TCanvas("effCanvas", "Sum Histogram", 1000, 1000);
    TLegend *h_leg = new TLegend(0.45, 0.65, 0.90, 0.90);
    h_leg->SetTextSize(0.03);
    std::vector<TH1F *> histos;
    dp = opendir("../rootfiles/");
    int colorCounter = 1;

    if (dp != nullptr)
    {
        while ((entry = readdir(dp)))
        {
            if (entry->d_type == DT_DIR)
                continue;
            string fileLocation = "../rootfiles/" + string(entry->d_name);
            string run_name = string(entry->d_name).substr(0, string(entry->d_name).size() - 5);
            cout << "Doing:  " << entry->d_name << " runName: " << run_name << endl;
            vector<string> runInformations = split(run_name.c_str(), "_");

            float energy = 1.4;
            if (runInformations[2].find("5") != std::string::npos)
            {
                energy = 5.2;
            }
            else if (runInformations[2].find("2") != std::string::npos)
            {
                energy = 2.6;
            }
            int position = -1;
            string ps = split(runInformations[3], ".")[0];
            string posraw = ps.substr(8, ps.length());
            //cout<< "Doing:  "<<posraw<<endl;
            position = stoi(posraw);
            cout << "pos: " << position << "  energy: " << energy << endl;
            //wls_e_5.2GeV_position0

            TApplication *ROOTapp;

            TFile *file = new TFile(fileLocation.c_str());

            TTree *treeEnergy;
            file->GetObject("Energy deposition", treeEnergy);
            Int_t volume;
            Double_t Edep;
            treeEnergy->SetBranchAddress("volume", &volume); //1 scintillator, 2 walls
            treeEnergy->SetBranchAddress("Edep", &Edep);

            //effCanvas->SetGrid();

            int n_bins = 50;
            // int n_bins = (int)(Xmax - Xmin)*1;
            TH1F *hist = new TH1F("allHist", "",n_bins, -100, 1200);

            hist->GetXaxis()->SetTitle("pulse-height(charge) [mv #times ns]");
            // h_vec[i]->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}]");
            hist->GetXaxis()->SetTitleOffset(1.3);
            treeEnergy->Draw("Edep>>allHist", "volume==1", "");
            hist->SetFillColorAlpha(colorCounter, 0.4);
            hist->SetLineColorAlpha(1, 0);
            hist->SetLineColor(colorCounter);
            h_leg->AddEntry(hist, Form("%1.1f GeV  mean: %1.1f MeV", energy, hist->GetMean()), "f");
            histos.push_back(hist);

            colorCounter++;

            cout << "" << endl;
        }
    }
    std::sort(histos.begin(), histos.end(),
              [](TH1F *a, TH1F *b) { return a->GetMean() > b->GetMean(); });
    for (auto h : histos)
        hs->Add(h);

    //  hs->Add(hist);

    hs->Draw("nostack");
    TAxis *yaxisP = hs->GetYaxis();
    TAxis *xaxisP = hs->GetXaxis();
    yaxisP->SetLabelSize(0.03);
    yaxisP->SetTitle("counts");
    yaxisP->SetTitleSize(0.03);
    xaxisP->SetLabelSize(0.03);
    xaxisP->SetTitle("energy deposition [MeV]");
    xaxisP->SetTitleSize(0.03);

    h_leg->Draw();

    effCanvas->SaveAs("out.pdf");

    closedir(dp);
    //	myapp->Run();
    return 0;
}