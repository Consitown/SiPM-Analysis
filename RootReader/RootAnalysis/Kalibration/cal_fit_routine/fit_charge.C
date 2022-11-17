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

//C, C++
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

const int n_peaks = 6;

// define a fit function, with single parameter corresponding to ptp distance
//ALLES Afterpulsing Zeugs weggelassen
//Oritniert an Super Funktion
//GP Fit-> Gain PArameter in GP-FIT
//Man sieht jetzt zwischen keinen 


//Man nimmt wieder roten Punkt. Auch wieder 
//Values.pdf-> kein funktionaler Aufwärtstrend mehr wie in Amp. Hinten schlechter weil Statistik der höheren Fit schlechter-> mehr Events nehmen
Double_t fitf(Double_t *x,Double_t *p)
{
  //0 - N0
  //1 - mu for poison
  //2 - muXT for crosstalk probability
  //3,4 -sogma0, sigma1
  //5,6 - G,B
  double sum = 0;
  for(int k=0;k<=10;k++)
  {
  	Double_t sigma0 = p[3];
  	Double_t sigma1 = p[4];
  	Double_t sigmaK = sqrt(sigma0*sigma0+k*sigma1*sigma1);
  	Double_t mu = p[1];
  	Double_t muXT = p[2];
  	Double_t G= p[5];
  	Double_t B= p[6];
  	
  	sum=sum+p[0]*mu*TMath::Power((mu+k*muXT),k-1)*TMath::Exp(-(mu+k*muXT))/TMath::Factorial(k) * (1./sqrt(2.*TMath::Pi())/sigmaK)*TMath::Exp(-TMath::Power( ( (x[0]-k*G+B)/sqrt(2)/sigmaK ),2 ));
  }
  return sum;
}

// fit function, sum of n_peaks Gaussians 
//Continious Fit
Double_t alt_f(Double_t *x,Double_t *p)
{
  Double_t gaus_single[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    gaus_single[i] = p[0+3*i]*TMath::Exp(-TMath::Power( p[1+3*i]-x[0],2) / (2*TMath::Power(p[2+3*i],2)) );
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

int main(int argc, char *argv[])
{
  // style options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.2); gStyle->SetStatH(0.1); // stats box size
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);


  /***** 
  __ INITIALIZE ___________________________
  *****/

  string path = "./"+(string)argv[1]+".root";
  int channel = atoi(argv[2]);
  TFile* file = new TFile(path.c_str());
  TTree* tree;
  file->GetObject("T",tree);
  
  int nBins = 250;
  int xmin = -80;
  int xmax = 500;

  // to show generalized poisson fit
  TCanvas* C1 = new TCanvas("C1","GP fit",800,800);
  TH1F* h = new TH1F("h","PE spectrum - integration window 35 ns",nBins,xmin,xmax);

  TString cut("");
  tree->Draw(Form("Integral_inRange[%d]>>h",channel),cut);//qS2 integration window width 55ns

  // to show alternative fit
  TCanvas *C2 = new TCanvas("C2","alt fit",800,800);
  TH1F* h2 = new TH1F("h2","",nBins,xmin,xmax);
  h2 = (TH1F*)h->Clone();
  
  /***** 
  __ FIT RANGE___________________________
  *****/

  // read off spectrum, ranges for individual peaks
  // std::vector<int> ranges = {-20,45,95,160,220,270}; // run 51
  // std::vector<int> ranges = {-25,50,110,160,230,300}; // run 52
  // std::vector<int> ranges = {-25,30,65,110,145,190}; // run 53
  // std::vector<int> ranges = {-25,30,65,110,150,190}; // run 53(1)
  // std::vector<int> ranges = {-10,25,60,90,130,160}; // run 56
  // std::vector<int> ranges = {-20,45,110,160,220,280,340}; // run 65
  std::vector<int> ranges = {-20,45,110,160,220,280,340}; // run 65
  // std::vector<int> ranges = {-20,30,75,120,170,220,260}; // run 77,76
  // std::vector<int> ranges = {-20,30,75,120,160,205,245}; // run 81,79
  // std::vector<int> ranges = {-20,35,75,115,160,195,235}; // run 82(1)
  // std::vector<int> ranges = {-20,25,65,100,130,165,215}; // run 82
  // std::vector<int> ranges = {-20,40,90,130,180,220,260}; // run 84
  // std::vector<int> ranges = {-20,40,90,130,180,220,260}; // run 85
  // std::vector<int> ranges = {-20,38,82,132,172,224,265}; // run 85(1)
  // std::vector<int> ranges = {-20,30,80,120,150,180,260}; // run 87

  // upper and lower range
    //Für alternative fit, schritte die man macht
  float l_range = 25;
  float u_range = 25;

  /***** 
  __ PRE FIT - Single Gauss fits___________________________
  //Gibt Startparameter für komplizierten Fit
  *****/
  
  // get starting value for gauss mean parameter
  float pos_peak[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    h->GetXaxis()->SetRange( h->GetXaxis()->FindBin(ranges[i]) ,h->GetXaxis()->FindBin(ranges[i+1]));
    pos_peak[i] = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  } 
  h->GetXaxis()->SetRange(0,nBins);

  TF1 * peak_single[n_peaks];
  Double_t par_single[3*n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  { 
    float range = 15; // symmetric fit range for individual peaks
    peak_single[i] = new TF1("peak","gaus",pos_peak[i]-range,pos_peak[i]+range);
    h->Fit("peak","RQ+");
    peak_single[i]->GetParameters(&par_single[0+3*i]);
  }

  /***** 
  __ GENERALIZED POISSON FIT ___________________________
  *****/
  //Range für kompletten Fit
  TF1 *f = new TF1("fitf",fitf,-20,400,7);
  f->SetParName(0,"N0");f->SetParameter(0,117600);
  f->SetParName(1,"#mu");f->SetParameter(1,0.95);
  f->SetParName(2,"#mu_{XT}");f->SetParameter(2,0.2);
  f->SetParName(3,"#sigma_{0 p.e.}");f->SetParameter(3,par_single[2]);
  f->SetParName(4,"#sigma_{1 p.e.}");f->SetParameter(4,par_single[5]);
  f->SetParName(5,"Gain");f->SetParameter(5,35);
  f->SetParName(6,"Base line");f->SetParameter(6,5);
  //f1->SetParName(6,"#sigma 2 p.e.");f1->SetParameter(6,25*sqrt(2.));
  //f1->SetParName(7,"#sigma 3 p.e.");f1->SetParameter(7,25*sqrt(3.));
  //f1->SetParName(8,"#sigma 4 p.e.");f1->SetParameter(8,25*sqrt(4.));
  //f1->Draw("same");
  h->GetYaxis()->SetTitle("#Entries");
  h->GetXaxis()->SetTitle("integral in mV#timesns");
  f->SetLineColor(3);
  h->Fit("fitf","RMQ");


  double calib_factor = f->GetParameter(5);
  double calib_factor_err = f->GetParError(5);
  double GP_chi2_ndof = f->GetChisquare()/f->GetNDF();

  printf("ch%d : 1pe dist: %1.2f ± %1.2f | red. chi2 = %1.2f\n",channel,calib_factor,calib_factor_err,GP_chi2_ndof );

  /***** 
  __ ALTERNATIVE FIT ___________________________
  *****/
//ERster Peak und letzter Peak und gehe dann nach Links und rechts.
  // TF1 *alt = new TF1("alt","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",pos_peak[0]-l_range,pos_peak[n_peaks-1]+u_range);
  TF1 *alt = new TF1("alt",alt_f, pos_peak[0]-l_range, pos_peak[n_peaks-1]+u_range, 3*n_peaks);
  alt->SetParameters(&par_single[0]);
  // for (int i = 0; i < 3*n_peaks+3; ++i)
  // {
  //  alt->SetParameter(i,par_single[i]);
  // }
  h2->Fit("alt","RQM");

  alt->SetLineColor(3);

  h2->GetYaxis()->SetTitle("#Entries");
  h2->GetXaxis()->SetTitle("integral in (mV #times ns)");
  h2->Draw();

  Double_t par_alt[3+n_peaks];
  Double_t par_alt_err[3+n_peaks];
  // alt->GetParameters(&par_alt[0]);
  for (int i = 0; i < 3*n_peaks; ++i)
  {
    par_alt[i] = alt->GetParameter(i);
    par_alt_err[i] = alt->GetParError(i);
  }

  double pos_alt[n_peaks], u_pos_alt[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    pos_alt[i] = alt->GetParameter(1+3*i);
    u_pos_alt[i] = alt->GetParError(1+3*i);
  }
  // results
  // printf("ch%d: ",channel );
  // for (int i = 0; i < n_peaks; ++i)
  // {
  //   printf("%1.2f±%1.2f ",pos_alt[i],u_pos_alt[i] );
  // }
  // printf("\n");

  // peak-to-peak distance
  double diff[n_peaks-1], u_diff[n_peaks-1], wght[n_peaks-1];
  double w_mean, u_w_mean, sum1, sum2;
  for (int i = 0; i < n_peaks-1; ++i)
  {
    diff[i] = pos_alt[i+1]-pos_alt[i];
    u_diff[i] = sqrt( u_pos_alt[i+1]*u_pos_alt[i+1] + u_pos_alt[i]*u_pos_alt[i] );
    wght[i] = 1./( u_diff[i] * u_diff[i] );
    sum1 += diff[i]*wght[i];
    sum2 += wght[i];
  }
  // weighted mean, neglecting correlation
  w_mean = sum1/sum2;
  u_w_mean = sqrt (1./sum2);

  /***** 
  __ COMPARE ALT <-> GEN.POISS ___________________________
  *****/

  double x_values[n_peaks-1];
  double x_values_err[n_peaks-1];
  for (int i = 0; i < n_peaks-1; ++i)
  {
    x_values[i] = i+1;
    x_values_err[i] = 0;
  }

  TCanvas *C3 = new TCanvas("C3","Compare Fit Results",700,500);
  gPad->SetGridx(); gPad->SetGridy();
  TGraphErrors * gr_alt = new TGraphErrors(n_peaks-1,x_values,diff,x_values_err,u_diff);
  gr_alt->SetName("gr_alt");
  gr_alt->SetTitle("PE Spectrum: Peak-To-Peak Distance");
  gr_alt->SetMarkerColor(4);
  gr_alt->SetMarkerStyle(21);
  gr_alt->GetXaxis()->SetTitle("peak number");
  gr_alt->GetYaxis()->SetTitle("peak-to-peak distance in (mV #times ns)");
  gr_alt->Draw("AP");


  double x_pint= 2.5;
  TGraphErrors * gr_comb = new TGraphErrors(1,&x_pint,&calib_factor,&x_values_err[1],&calib_factor_err);
  gr_comb->SetName("gr_comb");
  gr_comb->SetMarkerColor(2);
  gr_comb->SetMarkerStyle(11);
  gr_comb->Draw("P");

  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  // leg->SetHeader("The Legend Title");
  // leg->AddEntry(h1,"Histogram filled with random numbers","f");
  // leg->AddEntry("gr_alt","Function abs(#frac{sin(x)}{x})","l");
  leg->AddEntry("gr_comb","#splitline{average distance}{single parameter gauss fit}","lep");
  leg->AddEntry("gr_alt","#splitline{individual peak distances}{reg. cont. gauss fit}","lep");
  leg->Draw();

  /***** 
  __ Export Result ___________________________
  *****/

  string target_path = "./calib_histograms/charge/";
  string pdf_filename1 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_GP_fit.pdf";
  string pdf_filename2 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_contG_fit.pdf";
  string pdf_filename3 = target_path+(string)argv[1]+"_"+(string)argv[2]+"_values.pdf";
  string list_filename = target_path+"calib_factor.txt";

  FILE * factor_list;
  factor_list = fopen(list_filename.c_str(),"a");
  fprintf(factor_list, "ch%d : 1pe dist: %f ± %f | red. chi2 = %f\n",channel,calib_factor,calib_factor_err,GP_chi2_ndof );
  fclose(factor_list);

  gErrorIgnoreLevel = kError; // suppress root terminal output 
    C1->Print(pdf_filename1.c_str());
    C2->Print(pdf_filename2.c_str());
    C3->Print(pdf_filename3.c_str());
  
  return 0;

}


