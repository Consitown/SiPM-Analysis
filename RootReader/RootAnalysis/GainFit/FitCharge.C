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
#include <TGraphErrors.h>
#include <TError.h> // root verbosity level
#include <TVirtualFitter.h>

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
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;

// const int n_peaks = 5;
int g_argc;
char **g_argv;

// define a fit function, with single parameter corresponding to ptp distance
Double_t fitf(Double_t *x, Double_t *p)
{
  //0 - N0
  //1 - mu for poison
  //2 - muXT for crosstalk probability
  //3,4 -sogma0, sigma1
  //5,6 - G,B
  double sum = 0;
  for (int k = 0; k <= 20; k++)
  {
    Double_t sigma0 = p[3];
    Double_t sigma1 = p[4];
    Double_t sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
    Double_t mu = p[1];
    Double_t muXT = p[2];
    Double_t G = p[5];
    Double_t B = p[6];

    sum = sum + p[0] * mu * TMath::Power((mu + k * muXT), k - 1) * TMath::Exp(-(mu + k * muXT)) / TMath::Factorial(k) * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
  }
  return sum;
}

// fit function, sum of n_peaks Gaussians
Double_t alt_f(Double_t *x, Double_t *p)
{
  int n_peaks = atoi(g_argv[6]);
  Double_t gaus_single[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    gaus_single[i] = p[0 + 3 * i] * TMath::Exp(-TMath::Power(p[1 + 3 * i] - x[0], 2) / (2 * TMath::Power(p[2 + 3 * i], 2)));
  }

  double gaus_comb = 0;
  for (int i = 0; i < n_peaks; ++i)
  {
    gaus_comb = gaus_comb + gaus_single[i];
  }
  return gaus_comb;
}

/********************
__ FIT ROUTINE______
********************/

bool light = true;
int main(int argc, char **argv)
{
  ::g_argc = argc;
  ::g_argv = argv;
  // style options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatW(0.1);
  gStyle->SetStatH(0.1); // stats box size
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);

  /***** 
  __ INITIALIZE ___________________________
  *****/

  // parse arguments
  string run_name = (string)argv[1];
  int run_nr = atoi(argv[2]);
  string sipm_id = (string)argv[3];
  string wom_id = (string)argv[4];
  string sw_id = (string)argv[5];
  int n_peaks = atoi(argv[6]);
  string fileLocation = string(argv[7]);

  // __SETTINGS_____

  n_peaks = 6;
  int fitLimitEnd = 800;
  int fitLimitStart = -20; //40
  int xmin = -50;
  int xmax = 800;
  int nBins = (xmax - xmin) * 1;
  float range = 15;

  string in_dir = fileLocation;                                             // path to current root file
                                                                            // path to current root file
  string out_dir = Form("./calib_histograms/charge/%s/", run_name.c_str()); // to export plots and values
  string calib_factor_out_list = Form("calib_factor_run%d", run_nr);        // filename
  string BL_offset_out_list = Form("BL_offset_run%d", run_nr);              // filename values list

  if (!fs::is_directory(out_dir) || !fs::exists(out_dir))
  {                                // Check if src folder exists
    fs::create_directory(out_dir); // create src folder
  }

  // store results in .txt file
  string calib_list_filename = out_dir + calib_factor_out_list + ".txt";
  string BL_offset_list_filename = out_dir + BL_offset_out_list + ".txt";
  FILE *calib_factor_list;
  FILE *BL_offset_factor_list;
  calib_factor_list = fopen(calib_list_filename.c_str(), "w");
  BL_offset_factor_list = fopen(BL_offset_list_filename.c_str(), "w");

  // print date to file
  time_t now;
  time(&now);
  fprintf(calib_factor_list, "\ncharge calibration - %s\n", ctime(&now));
  fprintf(BL_offset_factor_list, "\ncharge calibration - %s\n", ctime(&now));

  // store results in vector
  vector<double> calib_factor_vec(8), calib_factor_err_vec(8), calib_factor_rel_err_vec(8);
  vector<double> calib_fit_rchi2_vec(8);
  vector<double> BL_offset_vec(8), BL_offset_err_vec(8);

  // fit range continuous fits
  float l_range = 20;
  float u_range = 25;

  // histogram settings
  // int nBins = 250;

  // y-range of calib factor plot
  // double l_values = 28, u_values =48;
  double l_values = 5, u_values = 9;
  // double l_values = 46, u_values =65;

  bool print_verb = 1;
  bool print_results_vec = 0;
  bool print_cf_array = 0;
  bool print_cf_err_array = 0;

  TCanvas *masterCanvas = new TCanvas("MC", "Overview", 1000, 1000);
  masterCanvas->Divide(3, 3);

  /*******************
  __ ANALYSIS ______
  ********************/

  printf("\nANALYZE CHARGE: run%d, SiPM: %s, %s, %s\n", run_nr, sipm_id.c_str(), wom_id.c_str(), sw_id.c_str());

  // loop over selected channels
  vector<int> channel_vec = {1, 2, 3, 4, 5, 6, 7, 8};
  // vector<int> channel_vec = {1};
  int n_channel = channel_vec.size();
  for (int i = 0; i < n_channel; ++i)
  {
    int channel = i;

    // __FIT RANGE__
    // read off spectrum, ranges for individual peaks
    std::vector<double> ranges;

    //Run 7 -> the one I used for Master Thesis: 7_calib_vb58_tune8700_pcbd
    switch (channel)
    {
    case 0:

      ranges = {25, 75, 120, 160, 200, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 226;
      fitLimitStart = 31; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 2.7;
      range = 10;
     
      //GOOD RESULTS

     /* minus
     
      ranges = {25, 75, 120, 160, 200, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 226;
      fitLimitStart = 30; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 3.302;
      range = 10;


    plus
      ranges = {25, 75, 120, 160, 200, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 226;
      fitLimitStart = 32; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 2.200; //2.2, 3.5
      range = 10;
*/





      break;
    case 1:

      ranges = {25, 70, 110, 155, 195, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 225;
      fitLimitStart = 33; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 3.3;
      range = 10;
      break;
    case 2:

      ranges = {25, 70, 110, 155, 195, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 251;
      fitLimitStart = 26; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 3.0;
      range = 10;

      break;

    case 3:

      ranges = {25, 75, 120, 160, 200, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 240;
      fitLimitStart = 35; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) *5.4;
      range = 10;
      

      break;
    case 4:
  

      ranges = {25, 75, 120, 160, 200, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 223;
      fitLimitStart = 33; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 2.6;
      range = 6;
   
      break;
    case 5:
    
      ranges = {25, 70, 110, 155, 195, 250, 300};
      n_peaks = 5;
      fitLimitEnd = 247;
      fitLimitStart = 40; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 3.2;
      range = 10;
   
      break;
    case 6:
     
      ranges = {-20, 20, 70, 120, 160, 200, 240, 280};
      n_peaks = 5;
      fitLimitEnd = 253;
      fitLimitStart = 41; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) * 3.49;
      range = 10;
  
      break;
    case 7:


      ranges = {-20, 20, 70, 120, 160, 200, 240, 280};
      n_peaks = 5;
      fitLimitEnd = 254;
      fitLimitStart = 36; //40
      xmin = -50;
      xmax = 400;
      nBins = (xmax - xmin) *3.7;
      range = 10;

      break;
    }
xmax =- 400;
    // open tree
    cout << "Doing: " << in_dir << endl;
    TFile *file = new TFile(in_dir.c_str());

    TTree *tree;
    file->GetObject("T", tree);

    // to show generalized poisson fit
    TCanvas *C1 = new TCanvas("C1", "GP fit", 1000, 800);
    TString h_title;
    if (sw_id == "none")
    {
      // h_title = Form("pulse-height spectrum - Hamamatsu S13360-%s, %s, ch%d", sipm_id.c_str(), wom_id.c_str(), channel);
    }
    else
    {
      // h_title = Form("pulse-height spectrum - Hamamatsu S13360-%s, %s, %s, ch%d", sipm_id.c_str(), wom_id.c_str(), sw_id.c_str(), channel);
    }

    TH1F *h = new TH1F("h", h_title, nBins, xmin, xmax);
    h->SetLineColorAlpha(kBlack, 0.7);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kBlack, 0.6);
    h->SetMarkerSize(0);
    TString cut("");
    //tree->Draw(Form("IntegralErrorP[%d]>>h", channel), cut);
    tree->Draw(Form("Integral[%d]>>h", channel), cut);
    //tree->Draw(Form("Integral[%d]>>h", channel), cut);

    h->Smooth(0);
    // to show alternative fit
    int hi = 800;
    if (light)
    {
      hi = 500;
    }

    TCanvas *C2 = new TCanvas("C2", "alt fit", 1000, hi);

    TH1F *h2 = new TH1F("h2", "", nBins, xmin, xmax);
    h2 = (TH1F *)h->Clone();
    h2->GetXaxis()->SetRangeUser(xmin, xmax);

    //C2->SetLogy();
    //C1->SetLogy();

    /***** 
    __ PRE FIT - Single Gauss fits___________________________
    *****/
    //TVirtualFitter::Fitter(h)->SetMaxIterations(100000);
    // TVirtualFitter::Fitter(h)->SetPrecision(1e-15);

    // get starting value for gauss mean parameter
    float pos_peak[n_peaks];
    for (int i = 0; i < n_peaks; ++i)
    {
      h->GetXaxis()->SetRangeUser(ranges[i], ranges[i + 1]);
      pos_peak[i] = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
    }
    h->GetXaxis()->UnZoom();

    // fit
    TF1 *peak_single[n_peaks];
    Double_t par_single[3 * n_peaks];
    for (int i = 0; i < n_peaks; ++i)
    {
      // run 47
      // float range = 12; // symmetric fit range for individual peaks
      peak_single[i] = new TF1("peak", "gaus", pos_peak[i] - range, pos_peak[i] + range);
      peak_single[i]->SetParameter(2, 10);
      peak_single[i]->SetParameter(1, pos_peak[i]);
      if (!light)
        h->Fit("peak", "RQ+");
      else
      {
        peak_single[i]->SetBit(TF1::kNotDraw);

        h->Fit("peak", "RQ");
      }
      peak_single[i]->GetParameters(&par_single[0 + 3 * i]);
    }

    /***** 
    __ GENERALIZED POISSON FIT * GAUSS FIT ___________________________
    *****/

    TF1 *f = new TF1("fitf", fitf, fitLimitStart, fitLimitEnd, 7);
    h->GetYaxis()->SetTitle("counts");
    h->GetXaxis()->SetTitle("#Lambda [mV #times ns]");
    if (light)
    {
      h->GetXaxis()->SetLabelSize(0.04);
      h->GetXaxis()->SetTitleSize(0.04);
      h->GetYaxis()->SetLabelSize(0.04);
      h->GetYaxis()->SetTitleSize(0.04);
    }

    f->SetLineColor(kGreen);
    f->SetLineStyle(1);
    f->SetNpx(1000);

    f->SetParName(0, "N0");
    f->SetParameter(0, h->Integral("width"));
    f->SetParName(1, "#mu");
    f->SetParameter(1, 3.5);
    f->SetParName(2, "#mu_{XT}");
    f->SetParameter(2, 0.2);
    f->SetParLimits(2, 0, 1);
    f->SetParName(3, "#sigma_{0 p.e.}");
    f->SetParameter(3, par_single[2]);
    f->SetParLimits(3, 0, 100);
    f->SetParName(4, "#sigma_{1 p.e.}");
    f->SetParameter(4, par_single[5]);
    f->SetParLimits(4, 0, 100);
    f->SetParName(5, "Gain");
    f->SetParameter(5, par_single[4] - par_single[1]);
    f->SetParName(6, "Base line");
    f->SetParameter(6, par_single[1]);

    if (!light)
      h->Fit("fitf", "RMQ");
    else
      h->Fit("fitf", "RMQ+");
    // C1->cd();
    // h->Draw("sameFUNC");

    // store fit results
    double calib_factor = f->GetParameter(5);
    double calib_factor_err = f->GetParError(5);
    double GP_chi2_ndof = f->GetChisquare() / f->GetNDF();
    double baseline = f->GetParameter(6);
    double baseline_err = f->GetParError(6);
    double norm = f->GetParameter(0);
    double norm_err = f->GetParError(0);
    double mu = f->GetParameter(1);
    double mu_err = f->GetParError(1);
    double mu_xt = f->GetParameter(2);
    double mu_xt_err = f->GetParError(2);
    double sig0 = f->GetParameter(3);
    double sig0_err = f->GetParError(3);
    double sig1 = f->GetParameter(4);
    double sig1_err = f->GetParError(4);

    TF1 *fCopy = new TF1("copy", fitf, -20, 1000, 7);
    fCopy->SetLineColor(kBlue);
    fCopy->SetLineStyle(2);
    fCopy->SetNpx(1000);

    fCopy->SetParName(0, "N0");
    fCopy->SetParameter(0, f->GetParameter(0));
    fCopy->SetParName(1, "#mu");
    fCopy->SetParameter(1, f->GetParameter(1));
    fCopy->SetParName(2, "#mu_{XT}");
    fCopy->SetParameter(2, f->GetParameter(2));
    fCopy->SetParName(3, "#sigma_{0 p.e.}");
    fCopy->SetParameter(3, f->GetParameter(3));
    fCopy->SetParName(4, "#sigma_{1 p.e.}");
    fCopy->SetParameter(4, f->GetParameter(4));
    fCopy->SetParName(5, "Gain");
    fCopy->SetParameter(5, f->GetParameter(5));
    fCopy->SetParName(6, "Base line");
    fCopy->SetParameter(6, f->GetParameter(6));

    /***** 
    __ ALTERNATIVE MULTI-GAUSS FIT ___________________________
    *****/

    C2->cd();
    // C2->SetTopMargin(0.90);
    C2->SetBottomMargin(0.12);
    C2->SetLeftMargin(0.125);
    C2->SetRightMargin(0.05);

    TF1 *alt = new TF1("alt", alt_f, pos_peak[0] - l_range, pos_peak[n_peaks - 1] + u_range, 3 * n_peaks);
    alt->SetNpx(1000);
    alt->SetLineColor(kRed);
    alt->SetLineStyle(5);
    alt->SetParameters(&par_single[0]);
    for (int j = 0; j < n_peaks; ++j)
    {
      alt->SetParameter(1 + 3 * j, baseline + j * calib_factor);
      alt->SetParameter(2 + 3 * j, sqrt(sig0 * sig0 + j * sig1 * sig1));
    }

    if (!light)
      h2->Fit("alt", "RQM");

    h2->GetYaxis()->SetTitle("counts");
    h2->GetXaxis()->SetTitle("#Lambda [mV #times ns]");
    if (light)
    {
      h2->GetXaxis()->SetLabelSize(0.04);
      h2->GetXaxis()->SetTitleSize(0.04);
      h2->GetYaxis()->SetLabelSize(0.04);
      h2->GetYaxis()->SetTitleSize(0.04);
    }

    h2->SetFillColorAlpha(kBlack, 0.1);
    h2->Draw("HISTE");

    masterCanvas->cd(i + 1);
    h2->Draw("HISTE");
    C2->cd();

    Double_t par_alt[3 + n_peaks];
    Double_t par_alt_err[3 + n_peaks];
    alt->GetParameters(&par_alt[0]);

    // draw individual gaussians
    TF1 *f_s[n_peaks];
    for (int i = 0; i < n_peaks; ++i)
    {
      f_s[i] = new TF1("s_peak", "gaus", 0, 1);
      f_s[i]->SetRange(par_alt[1 + 3 * i] - 100, par_alt[1 + 3 * i] + 100);
      f_s[i]->SetParameters(&par_alt[0 + 3 * i]);
      f_s[i]->SetLineColorAlpha(kRed, 0.4);
      f_s[i]->SetLineStyle(2);
      f_s[i]->SetNpx(1000);
      if (!light)
        f_s[i]->Draw("same");
    }
    if (!light)
    {
      fCopy->Draw("same");
      alt->Draw("same"); // draw fit above individual fit graphs
    }
    h->Draw("FUNCsame"); // also draw GP fit

    masterCanvas->cd(i + 1);
    if (!light)
      alt->Draw("same"); // draw fit above individual fit graphs
    h->Draw("FUNCsame"); // also draw GP fit

    C2->cd();

    double gauss_chi2_ndof = alt->GetChisquare() / alt->GetNDF();

    // custom legend
    TLegend *h_leg = new TLegend(0.65, 0.60, 0.99, 0.9);
    h_leg->SetTextSize(0.03);
    h_leg->AddEntry(h2, Form("#bf{pulse-charge spectrum}"), "lpef");
    if (!light)
      h_leg->AddEntry((TObject *)0, Form("entries: %1.0f", h2->GetEntries()), "");
    if (!light)
      h_leg->AddEntry((TObject *)0, Form("integration window 25 ns"), "");
    if (!light)
      h_leg->AddEntry(f, Form("Generalized Poisson #times Gaussian fit"), "l");
    h_leg->AddEntry((TObject *)0, Form("#chi^{2}_{red} = %1.3f", GP_chi2_ndof), "");
    // if (!light)
    h_leg->AddEntry((TObject *)0, Form("N = %1.f #pm %1.f", norm, norm_err), "");
    // if (!light)
    h_leg->AddEntry((TObject *)0, Form("#mu = %1.3f #pm %1.3f", mu, mu_err), "");
    // if (!light)
    h_leg->AddEntry((TObject *)0, Form("#lambda = %1.3f #pm %1.3f", mu_xt, mu_xt_err), "");
    if (!light)
      h_leg->AddEntry((TObject *)0, Form("#sigma_{0} = %1.3f #pm %1.3f", sig0, sig0_err), "");
    if (!light)
      h_leg->AddEntry((TObject *)0, Form("#sigma_{1} = %1.3f #pm %1.3f", sig1, sig1_err), "");
    h_leg->AddEntry((TObject *)0, Form("B = %1.3f #pm %1.3f", baseline, baseline_err), "");
        h_leg->AddEntry((TObject *)0, Form("#Lambda_{pixel} = %1.3f #pm %1.3f", calib_factor, calib_factor_err), "");

    if (!light)
    {
      h_leg->AddEntry(alt, Form("multi-Gaussian fit"), "l");
      h_leg->AddEntry((TObject *)0, Form("#chi^{2}/ndf = %1.2f", gauss_chi2_ndof), "");
      h_leg->AddEntry(f_s[1], Form("multi-Gaussian fit, individual peaks"), "l");
    }
    h_leg->Draw();
    masterCanvas->cd(i + 1);
    h_leg->Draw();
    C2->cd();
    double pos_alt[n_peaks], u_pos_alt[n_peaks];
    for (int i = 0; i < n_peaks; ++i)
    {
      pos_alt[i] = alt->GetParameter(1 + 3 * i);
      u_pos_alt[i] = alt->GetParError(1 + 3 * i);
    }

    // peak-to-peak distance
    double diff[n_peaks - 2], u_diff[n_peaks - 2], wght[n_peaks - 2];
    double w_mean, u_w_mean, sum1, sum2;
    for (int i = 0; i < n_peaks - 2; ++i)
    {
      diff[i] = (pos_alt[i + 2] - pos_alt[i]) / 2;

      u_diff[i] = 1. / 2 * sqrt(u_pos_alt[i + 2] * u_pos_alt[i + 2] + u_pos_alt[i] * u_pos_alt[i]);
      wght[i] = 1. / (u_diff[i] * u_diff[i]);
      sum1 += diff[i] * wght[i];
      sum2 += wght[i];
    }
    // weighted mean, neglecting correlation
    w_mean = sum1 / sum2;
    u_w_mean = sqrt(1. / sum2);

    /***** 
    __ SYSTEMATICS ___________________________
    *****/

    double var_sum = 0, calib_factor_sys, calib_factor_err_t;

    for (int i = 0; i < n_peaks - 2; ++i)
    {
      var_sum = TMath::Power(diff[i] - calib_factor, 2);
    }
    calib_factor_sys = sqrt(var_sum / (n_peaks - 2));
    calib_factor_err_t = sqrt(TMath::Power(calib_factor_sys, 2) + TMath::Power(calib_factor_err, 2));

    double rel_err = calib_factor_err_t / calib_factor;

    /***** 
    __ COMPARE MULTI-GAUSS <-> GEN.POISS*GAUSS ___________________________
    *****/

    TCanvas *C3 = new TCanvas("C3", "Compare Fit Results", 800, 500);
    C3->SetGrid();
    gPad->SetGridx();
    gPad->SetGridy();

    // construct x-coordinates
    double x_values[n_peaks - 1];
    double x_values_err[n_peaks - 1];
    for (int i = 0; i < n_peaks - 1; ++i)
    {
      x_values[i] = i + 1;
      x_values_err[i] = 0;
    }

    // graph: peak-to-peak distances
    TGraphErrors *gr_alt = new TGraphErrors(n_peaks - 2, x_values, diff, x_values_err, u_diff);
    if (sw_id == "none")
    {
      //gr_alt->SetTitle(Form("pulse-height spectrum: peak-to-peak distance - Hamamatsu S13360-%s, %s, ch%d", sipm_id.c_str(), wom_id.c_str(), channel));
    }
    else
    {
      //  gr_alt->SetTitle(Form("pulse-height spectrum: peak-to-peak distance - Hamamatsu S13360-%s, %s, %s, ch%d", sipm_id.c_str(), wom_id.c_str(), sw_id.c_str(), channel));
    }
    gr_alt->GetYaxis()->SetRangeUser(calib_factor - l_values, calib_factor + u_values);
    gr_alt->GetXaxis()->SetNdivisions(n_peaks - 2);
    gr_alt->SetName("gr_alt");
    gr_alt->SetMarkerColor(4);
    gr_alt->SetMarkerSize(0.5);
    gr_alt->SetMarkerStyle(21);
    gr_alt->GetXaxis()->SetTitle("peak number");
    gr_alt->GetYaxis()->SetTitle("peak-to-peak distance #Delta_{ptp} [mV #times ns]");
    gr_alt->Draw("AP");

    // construct data objects to draw
    double x_points[2];
    x_points[0] = 1;
    x_points[1] = n_peaks - 2;
    double x_points_err[2];
    x_points_err[0] = 0;
    x_points_err[1] = 0;
    double y_points[2];
    y_points[0] = calib_factor;
    y_points[1] = calib_factor;
    double y_points_err1[2];
    y_points_err1[0] = calib_factor_err_t;
    y_points_err1[1] = calib_factor_err_t;
    double y_points_err2[2];
    y_points_err2[0] = calib_factor_err;
    y_points_err2[1] = calib_factor_err;

    // graph: 1-sigma band: statistical error
    TGraphErrors *gr_comb = new TGraphErrors(2, x_points, y_points, x_points_err, y_points_err1);
    gr_comb->SetName("gr_comb");
    gr_comb->SetLineColor(2);
    gr_comb->SetMarkerColor(2);
    gr_comb->SetMarkerStyle(11);
    gr_comb->SetFillColor(2);
    gr_comb->SetFillColorAlpha(2, 0.2);
    gr_comb->Draw("l3");

    // graph: 1-sigma band: total error
    TGraphErrors *gr_comb2 = new TGraphErrors(2, x_points, y_points, x_points_err, y_points_err2);
    gr_comb2->SetName("gr_comb2");
    gr_comb2->SetLineColor(2);
    gr_comb2->SetMarkerColor(2);
    gr_comb2->SetMarkerStyle(11);
    gr_comb2->SetFillColor(2);
    gr_comb2->SetFillColorAlpha(2, 0.45);
    gr_comb2->Draw("l3");

    // print legend
    TLegend *leg = new TLegend(0.33, 0.56, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry("gr_alt", "peak-to-peak distance #Delta_{ptp,k} of k^{th} peak", "ep");
    leg->AddEntry("gr_comb2", Form("gain G_{fit,charge} = %1.2f mV #times ns", calib_factor), "l");
    leg->AddEntry((TObject *)0, Form("1#sigma-uncertainty bands:"), "");
    leg->AddEntry("gr_comb2", Form("stat: #sigma_{G_{fit,charge},stat} = %1.2f mV #times ns", calib_factor_err), "lf");
    leg->AddEntry("gr_comb", Form("stat+sys: #sigma_{G_{fit,charge},total} = (#sigma^{2}_{G_{fit,charge},stat} + #sigma^{2}_{#Delta_{ptp}})^{1/2} = %1.2f mV #times ns ", calib_factor_err_t), "lf");
    leg->AddEntry((TObject *)0, Form("rel. uncertainty #sigma_{G_{fit,charge},total} /  G_{fit,charge} : %1.2f %%", rel_err * 100), "");

    leg->Draw();

    /***** 
    __ Export Result ___________________________
    *****/

    if (print_verb)
    {
      // print results
      printf("ch%d : 1pe dist: %1.2f ± %1.2f ± %1.2f | red. chi2 = %1.2f | %1.2f %%\n", channel, calib_factor, calib_factor_err, calib_factor_sys, GP_chi2_ndof, rel_err * 100);
      printf("multi-gauss peak position: ");
      for (int j = 0; j < n_peaks; ++j)
      {
        printf("%1.1f  ", pos_alt[j]);
      }
      printf("\n");
    }
    if (print_cf_array)
    {
      printf("%f, ", calib_factor);
    }
    if (print_cf_err_array)
    {
      printf("%f, ", calib_factor_err);
    }

    // store in vector
    calib_factor_vec[i] = calib_factor;
    calib_factor_err_vec[i] = calib_factor_err;
    calib_factor_rel_err_vec[i] = rel_err * 100;
    calib_fit_rchi2_vec[i] = GP_chi2_ndof;
    BL_offset_vec[i] = baseline;
    BL_offset_err_vec[i] = baseline_err;

    // store in .txt file
    string pdf_filename2, pdf_filename3;
    if (sw_id == "none")
    {
      pdf_filename2 = out_dir + run_name + Form("_ch%d_fit.pdf", channel);
      pdf_filename3 = out_dir + run_name + Form("_ch%d_values.pdf", channel);
    }
    else
    {
      pdf_filename2 = out_dir + run_name + Form("_ch%d_%s_fit.pdf", channel, sw_id.c_str());
      pdf_filename3 = out_dir + run_name + Form("_ch%d_%s_values.pdf", channel, sw_id.c_str());
    }

    // C2->DrawClone("same");

    // export calibration factor
    fprintf(calib_factor_list, "ch%d : 1pe dist: %f ± %f ± %f ± %f  | red. chi2 = %f | %1.2f %%\n", channel, calib_factor, calib_factor_err, calib_factor_sys, calib_factor_err_t, GP_chi2_ndof, rel_err * 100);

    // export baseline offset
    fprintf(BL_offset_factor_list, "ch%d : pedestal peak position: %f ± %f | pedestal sigma = %f ± %f \n", channel, baseline, baseline_err, sig0, sig0_err);

    gErrorIgnoreLevel = kError; // suppress root terminal output
    // C1->Print(pdf_filename1.c_str());
    C2->Print(pdf_filename2.c_str());
    if (!light)
      C3->Print(pdf_filename3.c_str());
    gErrorIgnoreLevel = kUnset;

    delete C1;
    delete C2;
    delete C3;
  }
  // end, channel loop

  if (print_results_vec)
  {
    printf("calib factor:\n");
    for (int i = 0; i < 8; ++i)
    {
      printf("%f,", calib_factor_vec[i]);
    }
    printf("\n");

    printf("calib factor, statistical error:\n");
    for (int i = 0; i < 8; ++i)
    {
      printf("%f,", calib_factor_err_vec[i]);
    }
    printf("\n");

    printf("calib factor, relative stat+sys error:\n");
    for (int i = 0; i < 8; ++i)
    {
      printf("%f,", calib_factor_rel_err_vec[i]);
    }
    printf("\n");

    printf("baseline offset:\n");
    for (int i = 0; i < 8; ++i)
    {
      printf("%f,", BL_offset_vec[i]);
    }
    printf("\n");

    printf("baseline offset, error:\n");
    for (int i = 0; i < 8; ++i)
    {
      printf("%f,", BL_offset_err_vec[i]);
    }
    printf("\n");
  }

  fprintf(calib_factor_list, "%s={", run_name.c_str());
  for (int i = 0; i < 8; ++i)
  {
    fprintf(calib_factor_list, "%f/%f,", calib_factor_vec[i], calib_factor_err_vec[i]);
  }
  fprintf(calib_factor_list, "}\n");

  if (!light)
  {
    fprintf(calib_factor_list, "calib factor, statistical error:\n");
    for (int i = 0; i < 8; ++i)
    {
      fprintf(calib_factor_list, "%f,", calib_factor_err_vec[i]);
    }
    fprintf(calib_factor_list, "\n");

    fprintf(calib_factor_list, "calib factor, relative stat+sys error:\n");
    for (int i = 0; i < 8; ++i)
    {
      fprintf(calib_factor_list, "%f,", calib_factor_rel_err_vec[i]);
    }
    fprintf(calib_factor_list, "\n");

    fprintf(calib_factor_list, "calib fit, reduced chi2:\n");
    for (int i = 0; i < 8; ++i)
    {
      fprintf(calib_factor_list, "%f,", calib_fit_rchi2_vec[i]);
    }
    fprintf(calib_factor_list, "\n");

    fprintf(BL_offset_factor_list, "baseline offset:\n");
    for (int i = 0; i < 8; ++i)
    {
      fprintf(BL_offset_factor_list, "%f,", BL_offset_vec[i]);
    }
    fprintf(BL_offset_factor_list, "\n");

    fprintf(BL_offset_factor_list, "baseline offset, error:\n");
    for (int i = 0; i < 8; ++i)
    {
      fprintf(BL_offset_factor_list, "%f,", BL_offset_err_vec[i]);
    }
    fprintf(BL_offset_factor_list, "\n");
    fclose(BL_offset_factor_list);
  }

  fclose(calib_factor_list);

  string parent_dir = "./calib_histograms/charge/";
  string target_dir = parent_dir + run_name + "/";
  string overview_filename = out_dir + "overview_" + run_name + ".pdf";

  masterCanvas->Print(overview_filename.c_str());
  return 0;
}
