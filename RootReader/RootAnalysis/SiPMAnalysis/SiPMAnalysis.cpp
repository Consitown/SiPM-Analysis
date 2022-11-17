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
#include <Math/VavilovAccuratePdf.h>

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

/*******************
__ FUNCTIONS ______
*******************/

// Vavilov fit function
struct Vavilov_Func { 
   Vavilov_Func() {}

   double operator() (const double *x, const double *p) { 
      double kappa = p[0]; 
      double beta2 = p[1];
      return p[4]*( pdf.Pdf( (x[0]-p[2])/p[3], kappa,beta2) );
   }

   ROOT::Math::VavilovAccurate pdf; 
};
/**
 * Helper function to process the run Informations from the name
 * */
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

int main(int argc, char *argv[]){

	bool store_result = 1;

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	gStyle->SetLineScalePS(1); // export high resolution .pdf
	float alpha = 0.1;
	int lineW = 2;

	string runName = (string)argv[1];
	int run_nr = atoi(argv[2]);
	string inDir = string(argv[3]);


	//Open Root File
	TFile *file = new TFile(inDir.c_str());
	if (file->IsZombie())
	{
		cout << "PROBLEM with the initialization of the output ROOT ntuple "
			 << inDir << ": check that the path is correct!!!"
			 << endl;

		exit(-1);
	}
	TTree *tree;
	file->GetObject("T", tree);






	

	//Du willst an die Informationen aus dem RunName, wie macht man das?
	/*
	Ein Teil der Infos werden ja bereits in der Shell Datei geparsed und dann hier schön als Argument übergeben ("runName,run_nr, inDir")
	Das könntest du eigentlich mit allen Informationen (Position,Energie) machen. Hab Ich aber nicht, weil nie gebraucht.
	Stattdessen kannst du deswegen aus dem runName (sowas wie: 14_pos4_angle0_e14_ch32) auch hier in C++ die Infos holen.
	Beispiel Energie, Position:
	14_pos4_angle0_e14_ch32 -> Split by "_" -> 14,pos4,angle0,e14,ch32.root   -> Energie= Index 3, siehe hier: 



	ACHTUNG: Ich konnte das hier nicht testen, du musst also selbst testen ob das geht. Denk daran: Loggen=Debugging

	*/
	vector<string> runInformations = split(runName.c_str(), "_");
	cout << "Doing:  " << runName << endl; //Debug

	float energy = 1.4;
	if (runInformations[3].find("5") != std::string::npos)
	{
		energy = 5.2;
	}
	else if (runInformations[3].find("2") != std::string::npos)
	{
		energy = 2.6;
	}
	int pos = -1;

	string ps =runInformations[1];
	string posraw = ps.substr(2, ps.length()); //entferne "pos" aus "pos4"
	pos = stoi(posraw); //in int umwandeln
	cout << "pos: " << pos << "  energy: " << energy << endl; //Schon wieder ein DEBUG CALL -> wenn dieser hier nichtmehr in deiner console erscheint, dann wird wohl zwischen "Doing" und hier etwas falsches passiert sein




	//string runName = (string)argv[1];
	//int run_nr = atoi(argv[2]);
	//int pos = atoi(argv[3]);
	//string WC_v = (string)argv[4]; //unnötig
	//string particle = (string)argv[5]; //unnötig, nur elektronen
	//string target_wom = (string)argv[6]; //vermutlich auch unnötig, nur WOM D
	string target_wom="WOM_D";
	string particle="electron"; //wird nur ausgedruckt, komplett egal


	int scan_nr; //Was ist das? Es wird irgendwie als Parameter übergeben, wenn das irgendwie hardcoded ist, dann kannst du das auch hier hardcoden ODER den Namen der Datei ändern, z.B. 14_pos4_angle0_e14_ch32.root in 14_pos4_angle0_e14_ch32_sc5.root
				// und dann wie oben sc als "scan_nr" parsen. Alternativ hier ein switch statement rein und diese Nummer anhand der Runnummer anpassen.
	if (argc == 8){scan_nr = atoi(argv[7]);}
	else{scan_nr=-99;}
	

	// open tree
	string out_dir = "./ly_histograms/charge/"+target_wom+"/"; // dir. to export plots and txt file
	string out_pdf_name = Form("light_yield_%s_pos%d_%s_run%d.pdf",target_wom.c_str(),pos,particle.c_str(),run_nr) ;// histogram pdf filename
	string out_values_list_name = Form("mpv_%s.txt",target_wom.c_str()) ; // MPV list filename
	string out_values_err_list_name = Form("mpv_err_%s.txt",target_wom.c_str()) ; // MPV err list


	// SiPM calibration uncertainty
	vector<double> cal_err_WOM_A = {0.889009,0.309643,0.380993,0.462777,0.870155,0.152640,1.006724};	
	vector<double> cal_err_WOM_B = {0.710870,0.342268,1.733274,2.257234,0.363562,0.386964,0.822406,0.273501};
	vector<double> cal_err_WOM_C = {0.429751,0.498424,0.398480,0.280804,0.782851,0.644864,0.212400};
	vector<double> cal_err_WOM_D = {1.344291,0.904895,0.367969,0.437667,1.559801,1.500992,0.312632,0.802277}; 
	// JOSCHA: Are these the values I need to use? Since from my understanding, they are for the same SiPM array but for a different measurement.
	/**
	 * Spielt keine Rolle, diese Werte stimmen sowieso nicht.  Er benutzt sie außerdem nichtmal. Wenn du Visual Studio Code benutzt kannst du Strg+Linksclick (alternative einfach danach suchen) und siehst, dass das a
	 * alles nur überbleibsel sind. Die Hälfte des Skripts ist eigentlich Müll. Kannst du also vielleicht mal neu machen :)
	 * */



	vector<double> cal_err_vec;
	if (target_wom=="WOM_A"){cal_err_vec = cal_err_WOM_A;}
	else if(target_wom=="WOM_B"){cal_err_vec = cal_err_WOM_B;}
	else if(target_wom=="WOM_C"){cal_err_vec = cal_err_WOM_C;}
	else if(target_wom=="WOM_D"){cal_err_vec = cal_err_WOM_D;}

	string spec_charge,exp_file_spec;
	//spec_charge = "_int";
	spec_charge = "";
	//exp_file_spec = "charge"; wird nichtmehr gebraucht?!
	
	bool is_charge=1;

	// default ranges
	int upperA = 100, upperB = 100, upperA_sum = 450, upperB_sum = 450;

	// set x range depending on measurement position
	if (pos == 1)
	{
		upperA = 70; upperA_sum = 180 ;
		upperB = 70; upperB_sum = 50;
	}
	if (pos == 2)
	{
		upperA = 70; upperA_sum = 220 ;
		upperB = 30; upperB_sum = 50;
	}
	if (pos == 4)
	{
		upperA = 70; upperA_sum = 120 ;
		upperB = 80; upperB_sum = 200;
	}
	if (pos == 5)
	{
		upperA = 70; upperA_sum = 120 ;
		upperB = 80; upperB_sum = 200;
	}
	if (pos == 5 && particle == "electron")
	{
		upperA = 25; upperA_sum = 300 ;
		upperB = 40; upperB_sum = 250;
	}
	if (pos == 6)
	{
		upperA = 70; upperA_sum = 120;
		upperB = 80; upperB_sum = 200;
	}
	if (pos == 7)
	{
		upperA = 60; upperA_sum = 120 ;
		upperB = 150; upperB_sum = 450;
	}
	if (pos == 7 && particle == "electron")
	{
		upperA = 15; upperA_sum = 50 ;
		upperB = 150; upperB_sum = 800;
	}
	if (pos == 8)
	{
		upperA = 60; upperA_sum = 120;
		upperB = 150; upperB_sum = 450;
	}
	if (pos == 9)
	{
		upperA = 60; upperA_sum = 120 ;
		upperB = 150; upperB_sum = 450;
	}
	if (pos == 9 && particle == "electron")
	{
		upperA = 20; upperA_sum = 100 ;
		upperB = 150; upperB_sum = 800;
	}
	if (pos == 10)
	{
		upperA = 60; upperA_sum = 150 ;
		upperB = 60; upperB_sum = 150;
	}
	if (pos == 11)
	{
		upperA = 70; upperA_sum = 110 ;
		upperB = 80; upperB_sum = 200;
	}
	if (pos == 13 )
	{
		upperA = 60; upperA_sum = 100 ;
		upperB = 100; upperB_sum = 250;
	}
	if (pos == 14)
	{
		upperA = 70; upperA_sum = 120 ;
		upperB = 80; upperB_sum = 200;
	}
	if (pos == 15)
	{
		upperA = 50; upperA_sum = 120 ;
		upperB = 50; upperB_sum = 120;
	}
	if (pos == 16)
	{
		upperA = 50; upperA_sum = 100;
		upperB = 125; upperB_sum = 300;
	}
	if ( is_charge==1 && run_nr >=88 && run_nr <=103 )
	{
		upperA = 30; upperA_sum = 60 ;
		upperB = 500; upperB_sum = 3000;
	}
	if ( is_charge==0 && run_nr >=88 && run_nr <=103 )
	{
		upperA = 30; upperA_sum = 60 ;
		upperB = 200; upperB_sum = 2000;
	}
		if ( run_nr >=104 && run_nr <=107 )
	{
		upperA = 30; upperA_sum = 60 ;
		upperB = 200; upperB_sum = 1000;
	}


	/***** 
	__ PE DISTRIBUTION _______
	*****/

	int n_ch ; // number of SiPM channels + sum data
	if (target_wom=="WOM_A" || target_wom=="WOM_C")
	{n_ch = 8;}
	if (target_wom=="WOM_B" || target_wom=="WOM_D")
	{n_ch = 9;}

	TH1F * h_vec[n_ch];
	TLine * ln_vec[n_ch];
	TF1 * fit_vec_g[n_ch], * fit_vec_l[n_ch], * fit_vec_v[n_ch];;
	double amp_max[n_ch], peak[n_ch], peak_err[n_ch], rchi2_g[n_ch];
	double peak_v[n_ch], peak_mean_v[n_ch], rchi2_v[n_ch];
	int h_entries[n_ch];

	TCanvas *C1;
	C1 = new TCanvas("ly_dist","",600,500);
	C1->Divide(3,3);

	//___ LOOP OVER CHANNELS & SUM ___
	for (int i = 0; i < n_ch; ++i)
	{
		TString h_name1, h_title1, draw_cmnd1, cut_cmnd1;


		/***** 
		__ SHOW Light Yield ______________________________
		*****/

		int Xmin = -5, Xmax;
		double lan_range_lo = 0.5, lan_range_hi;
		float frac_l=0.5,frac_u=0.4;




		/**
		 * Also:
		 * Die Kommandos haben die Form:
		 * draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i,i);
		 * spec_ wird irgendwie "specification_" bedeutet haben und wurde zur unterscheidung zwischen Amplitude und Ladung benutzt wurden
		 * ABER: Wie du hoffentlich noch aus meiner Arbeit weißt (und der von Julian) macht die Amplitude keinen Sinn (Probleme bei der Kalibration). Ergo wurden nur Integrale (=Charge) genutzt.
		 * Öffne Ubuntu, tippe: root, danach bist du im ROOT Shell, dann gibst du ein "new TBrowser". Es öffnet sich ein Fenster (wenn nicht, starte XMing), dort suchst du eine ROOT datei von mir/dir (eine neue halt).
		 * Doppelclick und du siehst einen Tree mit dem Namen "T". Dort drauflicken und du siehst die Branches. 
		 * 
		 * Dieser Befehl:
		 *  draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i,i);
		 * Erzeugt einen String von der Form: Integral_int[0]>>h0    -> im TBrowser den du gerade hoffentlich geöffnet hast, siehst du aber, dass es keinen Branch mit dem Namen: Integral_int gibt. Es gibt nur einen mit dem
		 * Namen: Integral. D.h. diese Befehle machen keinen Sinn. Ferner macht das schon von dem Namen keinen Sinn. Integral_int(egral)? 
		 * 
		 * Also Fazit: Es muss heißen: draw_cmnd1.Form("Integral[%d]>>h%d",i,i);  oder wahlweise alles so lassen und spec_charge=""; (habe Ich oben mal so umgestellt)
		 * 
		 * */


			/**
			 * JOSCHA: Judging by your (Jans) script, I don't get why it uses i+7. I tried just using i, which yields different but equally terrible results.
			 * 
			 * Das kommt daher, dass man in einem Branch die Daten aus mehreren WOMS speichert. Die Channel 0-7 sind dann z.B. für die Integrale von WOM A,C und die von 7-15 für WOM D,B
			 * Das ist aber alles nichtmehr so. Integral[channel] für WOM D=> channel 0-7. Die Ordnung wird im read.C angegeben (Zeile 670, WOM D,C,A,B -> 3,2,0,1, dort siehst du, dass i<8 wird WOM_ID=3 also WOM D zugewiesen)
			 * 
			 * Du musst das hier mal alles anpassen und die Logik überprüfen.
			 * 
			 * 
			 * */




		if (i<n_ch-1) // SiPM channels
		{	
			C1->cd(i+1);
			h_name1.Form("h%d",i);
			if (target_wom=="WOM_A")
			{
				h_title1.Form("WOM-A, run%d, pos%d, %s, ch%d",run_nr,pos,particle.c_str(),i);
				draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i,i); //JOSCHA: These are where it used to say "chPE%s[%d]>>h$d". With is_charge set to 0 (line 60) it reads chPE, if set to 1 it reads chPE_int.
				cut_cmnd1.Form("Integral%s[%d]>0",spec_charge.c_str(),i);
				Xmax = upperA;
				lan_range_hi = Xmax-50;
			}
			else if(target_wom=="WOM_C")
			{
				h_title1.Form("WOM-C, run%d, pos%d, %s, ch%d",run_nr,pos,particle.c_str(),i);
				draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i,i);
				cut_cmnd1.Form("Integral%s[%d]>0",spec_charge.c_str(),i);
				Xmax = upperA;
				lan_range_hi = Xmax-50; 
			}
			else if(target_wom=="WOM_B")
			{
				h_title1.Form("WOM-B, run%d, pos%d, %s, ch%d",run_nr,pos,particle.c_str(),i);
				draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i+7,i);
				cut_cmnd1.Form("Integral%s[%d]>0",spec_charge.c_str(),i+7);

				Xmax = upperB;
				lan_range_hi = Xmax-50; 
			}
			else if(target_wom=="WOM_D")
			{
				h_title1.Form("WOM-D, run%d, pos%d, %s, ch%d",run_nr,pos,particle.c_str(),i);
				draw_cmnd1.Form("Integral%s[%d]>>h%d",spec_charge.c_str(),i+7,i); //JOSCHA: Judging by your (Jans) script, I don't get why it uses i+7. I tried just using i, which yields different but equally terrible results.
				cut_cmnd1.Form("Integral%s[%d]>0",spec_charge.c_str(),i);
				Xmax = upperB;
				lan_range_hi = Xmax-50; 
			}			
		}
		else if (i == n_ch-1) // sum; JOSCHA: This part doesn't seem to do what it should either, but I guess I can just ignore it for now since I don't need the sum at the moment.
		{ //Vermutlich wollten sie hier irgendwas mit WOM-Summen machen. Die haben aber mittlerweile auch einen anderen Namen, read.C Zeile 1000.
			C1->cd(9);
			h_name1.Form("h%d",i);
			if (target_wom=="WOM_A")
			{
				h_title1.Form("WOM-A, run%d, pos%d, %s, sum",run_nr,pos,particle.c_str());
				draw_cmnd1.Form("PE_WOM1%s>>h%d",spec_charge.c_str(),i);
				cut_cmnd1.Form("PE_WOM1%s>0",spec_charge.c_str());
				Xmax = upperA_sum; 
				lan_range_hi = Xmax-75;
			}
			else if(target_wom=="WOM_C")
			{
				h_title1.Form("WOM-C, run%d, pos%d, %s, sum",run_nr,pos,particle.c_str());
				draw_cmnd1.Form("PE_WOM1%s>>h%d",spec_charge.c_str(),i);
				cut_cmnd1.Form("PE_WOM1%s>0",spec_charge.c_str());
				Xmax = upperA_sum; 
				lan_range_hi = Xmax-75;
			}
			else if(target_wom=="WOM_B")
			{
				h_title1.Form("WOM-B, run%d, pos%d, %s, sum",run_nr,pos,particle.c_str());
				draw_cmnd1.Form("PE_WOM2%s>>h%d",spec_charge.c_str(),i);
				cut_cmnd1.Form("PE_WOM2%s>0",spec_charge.c_str());
				Xmax = upperB_sum;
				lan_range_hi = Xmax-75; 
			}
			else if(target_wom=="WOM_D")
			{
				h_title1.Form("WOM-D, run%d, pos%d, %s, sum",run_nr,pos,particle.c_str());
				draw_cmnd1.Form("PE_WOM2%s>>h%d",spec_charge.c_str(),i);
				cut_cmnd1.Form("PE_WOM2%s>0",spec_charge.c_str());
				Xmax = upperB_sum; 
				lan_range_hi = Xmax-75;

			}			
		}
		
		// canvas style
		// gPad->SetLogy();
		// gPad->SetGridx(); gPad->SetGridy();
		gPad->SetRightMargin(0.00);
		gPad->SetLeftMargin(.12);
		
		// draw histogram
		// h_vec[i] = new TH1F(h_name1,h_title1,(Xmax - Xmin)*1,Xmin-0.5,Xmax-0.5);
		h_vec[i] = new TH1F(h_name1,h_title1,(Xmax - Xmin)*2,Xmin-0.75,Xmax-0.75);
		h_vec[i]->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}]");
		h_vec[i]->GetXaxis()->SetTitleOffset(1.3);
		h_vec[i]->GetYaxis()->SetTitle("entries");
		h_vec[i]->SetLineColorAlpha(kBlack,0.7);
		h_vec[i]->SetFillColorAlpha(kBlack,alpha);
		h_vec[i]->SetMarkerStyle(8);
		h_vec[i]->SetMarkerSize(0.2);
		h_vec[i]->SetMarkerColorAlpha(kBlack,0.6);
		tree->Draw(draw_cmnd1,"","HISTE");

		/***** 
		__ Maximum LY in range_________________
		*****/
		double h_mean =  h_vec[i]->GetMean();
    	double h_max_y = h_vec[i]->GetMaximum();

    	// get histogram maximum value & position
		h_vec[i]->GetXaxis()->SetRange(0,Xmax-30);
		amp_max[i] = h_vec[i]->GetXaxis()->GetBinCenter(h_vec[i]->GetMaximumBin());
		float max = h_vec[i]->GetMaximum();
		int lower_bin = h_vec[i]->GetMaximumBin();
		int upper_bin = h_vec[i]->GetMaximumBin();
		h_vec[i]->GetXaxis()->UnZoom();
		// printf("ch%d max amp: %f\n",i,amp_max[i]);		

		ln_vec[i] = new TLine(amp_max[i],0.1,amp_max[i],h_max_y);
		ln_vec[i]->SetLineColor(8);
		ln_vec[i]->SetLineWidth(1);
		ln_vec[i]->SetLineStyle(5);
		ln_vec[i]->Draw("same");

		/***** 
		__ FIT PEAK: GAUSS ______________________________
		*****/

		// fit range ± 0.5x maximum
		// limit lower range to +0.5 Npe
		while (h_vec[i]->GetBinContent(lower_bin) > max*frac_l && (h_vec[i]->GetBinLowEdge(lower_bin))>0.5 ) {lower_bin--;}
		while (h_vec[i]->GetBinContent(upper_bin) > max*frac_u) {upper_bin++;}

		float r1,r2;
		r1 = h_vec[i]->GetBinLowEdge(lower_bin);
		r2 = h_vec[i]->GetBinCenter(upper_bin)+0.5;

		// fit_vec_g[i] = new TF1("g_fit","gaus",0.5,30);
		fit_vec_g[i] = new TF1("g_fit","gaus",r1,r2);
		fit_vec_g[i]->SetParameter(1,10);
		fit_vec_g[i]->SetLineColor(4);
		fit_vec_g[i]->SetLineStyle(1);
		// fit_vec_g[i]->SetLineStyle(5);
		fit_vec_g[i]->SetNpx(1000); // draw function with high resolution 
		h_vec[i]->Fit("g_fit","RQM");
		h_vec[i]->Draw("sameFUNC");

		Double_t par_g[3]; // to store fit results
		fit_vec_g[i]->GetParameters(&par_g[0]);

		// results
		peak[i] = fit_vec_g[i]->GetParameter(1);
		peak_err[i] = fit_vec_g[i]->GetParError(1);
		
		// // take calibration uncertainty into account
		// --->> NO! would only introduce an error bar that is same for all points 
		// 		 wrt a given channels...
		// double abs_cal_err;
		// if (i<n_ch-1){	
		// 	abs_cal_err = cal_err_vec[i]/100*peak[i];
		// 	peak_err[i] = sqrt( peak_err[i]*peak_err[i] + abs_cal_err*abs_cal_err);
		// }

		double sigma_g = fit_vec_g[i]->GetParameter(2);
		double sigma_g_err = fit_vec_g[i]->GetParError(2);
		double rchi2_g = fit_vec_g[i]->GetChisquare()/fit_vec_g[i]->GetNDF();

		// draw line for MPV
		TLine * ln_gauss_mpv = new TLine(peak[i],0,peak[i],h_max_y);
		ln_gauss_mpv->SetLineColor(4);
		ln_gauss_mpv->SetLineStyle(3);
		ln_gauss_mpv->Draw("same");

		// custom histogram legend with fit results
		TLegend *h_lan_leg = new TLegend(0.53,0.50,1.0,0.9);
		h_lan_leg->SetTextSize(0.04);
		h_lan_leg->AddEntry(h_vec[i],Form("#bf{data}"),"elpf");
		h_lan_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i]->GetEntries()),"");
		h_lan_leg->AddEntry(ln_vec[i],Form("distribution max.: %1.1f N_{pe} ",amp_max[i]),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{L} = %1.2f #pm %1.2f ",lan_sigma,lan_sigma_err),"");
		h_lan_leg->AddEntry(fit_vec_g[i],Form("Gaussian fit"),"l");
		h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.1f",rchi2_g),"");
		h_lan_leg->AddEntry(ln_gauss_mpv,Form("MPV = %1.2f #pm %1.2f ",peak[i],peak_err[i]),"l");
		h_lan_leg->AddEntry((TObject*)0,Form("#sigma = %1.2f #pm %1.2f",sigma_g,sigma_g_err),"");
		h_lan_leg->Draw();


		/***** 
		__ FIT PEAK: LANDAU ______________________________
		*****/

		fit_vec_l[i] = new TF1("fit_l","landau",0.5,lan_range_hi);
		fit_vec_l[i]->SetLineColor(2);
		fit_vec_l[i]->SetNpx(1000); // draw function with high resolution 
		// fit_vec_l[i]->SetParameters(&par_g[0]);

		h_vec[i]->Fit("fit_l","RQMN");
		// h_vec[i]->Draw("sameFUNC");

		 // store fit results
		double mpshift  = -0.22278298;
		double lan_sigma = fit_vec_l[i]->GetParameter(2);
		double lan_sigma_err = fit_vec_l[i]->GetParError(2);
		float lan_mpv = fit_vec_l[i]->GetParameter(1);
		float lan_mpv_err = fit_vec_l[i]->GetParError(1);
		float lan_rchi2 = fit_vec_l[i]->GetChisquare()/fit_vec_l[i]->GetNDF();
		
		double lan_mpv_c = lan_mpv + mpshift*lan_sigma ;
		double lan_mpv_c_err = sqrt(lan_mpv_err*lan_mpv_err + (mpshift*lan_sigma_err)*(mpshift*lan_sigma_err)) ;
		
		

		// // draw line for MPV
		// TLine * ln_lan_mpv = new TLine(lan_mpv_c,0,lan_mpv_c,h_max_y);
		// ln_lan_mpv->SetLineColor(2);
		// ln_lan_mpv->SetLineStyle(2);
		// ln_lan_mpv->Draw("same");

		// // custom histogram legend with fit results
		// TLegend *h_lan_leg = new TLegend(0.53,0.30,1.0,0.9);
		// h_lan_leg->SetTextSize(0.04);
		// h_lan_leg->AddEntry(h_vec[i],Form("#bf{data}"),"elpf");
		// h_lan_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i]->GetEntries()),"");
		// h_lan_leg->AddEntry(ln_vec[i],Form("dist. max. = %1.2f ",amp_max[i]),"l");
		// h_lan_leg->AddEntry(fit_vec_l[i],Form("Landau fit"),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}_{L}/ndf = %1.1f",lan_rchi2),"");
		// h_lan_leg->AddEntry(ln_lan_mpv,Form("MPV_{L} = %1.2f #pm %1.2f ",lan_mpv_c,lan_mpv_c_err),"l");
		// // h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{L} = %1.2f #pm %1.2f ",lan_sigma,lan_sigma_err),"");
		// h_lan_leg->AddEntry(fit_vec_g[i],Form("Gauss fit"),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#chi^{2}_{G}/ndf = %1.1f",rchi2_g),"");
		// h_lan_leg->AddEntry(ln_gauss_mpv,Form("MPV_{G} = %1.2f #pm %1.2f ",peak[i],peak_err[i]),"l");
		// h_lan_leg->AddEntry((TObject*)0,Form("#sigma_{G} = %1.2f #pm %1.2f ",sigma_g,sigma_g_err),"");
		// h_lan_leg->Draw();

		/***** 
		__ FIT PEAK: VAVILOV ______________________________
		*****/

		Vavilov_Func * func = new Vavilov_Func();
		fit_vec_v[i] = new TF1("f1",func, 0.5,lan_range_hi,5,"Vavilov_Func");
		fit_vec_v[i]->SetParNames("kappa","beta2","mean","sigma","Amp"); 
		fit_vec_v[i]->SetParameters(0.01,0.03,fit_vec_l[i]->GetParameter(1),fit_vec_l[i]->GetParameter(2),h_vec[i]->GetEntries());
		fit_vec_v[i]->SetParLimits(0,0.01,10); // Vavilov model valid only in this kappa range
		// fit_vec_v[i]->SetParLimits(1,0,1);
		fit_vec_v[i]->FixParameter(1,1); // fix beta2 to "light speed", no dependence this parameter observed
		fit_vec_v[i]->SetNpx(1000); // draw function with high resolution 
		fit_vec_v[i]->SetLineColor(kGreen);

		h_vec[i]->Fit("f1","RQMN");

		// h_vec[i]->Draw("sameFUNC");

	    // store fit results
	    double vav_mean = fit_vec_v[i]->GetParameter(2);
	    double vav_max = fit_vec_v[i]->GetMaximumX(vav_mean-5,vav_mean);
	    double vav_rchi2 =  fit_vec_v[i]->GetChisquare()/fit_vec_v[i]->GetNDF();

	    peak_v[i] = vav_max;
		peak_mean_v[i] = vav_mean;
		rchi2_v[i] = vav_rchi2;
		h_entries[i] = h_vec[i]->GetEntries();

		// // draw lines for MPV and mean
		// TLine * ln_vav_mean = new TLine(vav_mean,0,vav_mean,h_max_y);
		// ln_vav_mean->SetLineColor(8);
		// ln_vav_mean->SetLineStyle(2);
		// ln_vav_mean->Draw("same");

		// TLine * ln_vav_max = new TLine(vav_max,0,vav_max,h_max_y);
		// ln_vav_max->SetLineColor(9);
		// ln_vav_max->SetLineStyle(9);
		// ln_vav_max->Draw("same");

		// h_vec[i]->Draw("sameFUNC");

		// // custom histogram legend with fit results
		// TLegend *h_vav_leg = new TLegend(0.50,0.60,0.9,0.9);
		// h_vav_leg->AddEntry(h_vec[i],Form("#bf{data}"),"elpf");
		// h_vav_leg->AddEntry((TObject*)0,Form("entries = %1.f",h_vec[i]->GetEntries()),"");
		// h_vav_leg->AddEntry(fit_vec_v[i],Form("Vavilov fit"),"l");
		// h_vav_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.1f",vav_rchi2),"");
		// h_vav_leg->AddEntry((TObject*)0,Form("#kappa = %1.3f #pm %1.3f",fit_vec_v[i]->GetParameter(0),fit_vec_v[i]->GetParError(0)),"");
		// h_vav_leg->AddEntry((TObject*)0,Form("#beta^{2} = %1.f",fit_vec_v[i]->GetParameter(1)),"");
		// h_vav_leg->AddEntry(ln_vav_mean,Form("mean = %1.2f #pm %1.2f",vav_mean,fit_vec_v[i]->GetParError(2)),"l");
		// h_vav_leg->AddEntry((TObject*)0,Form("#sigma = %1.3f #pm %1.3f",fit_vec_v[i]->GetParameter(3),fit_vec_v[i]->GetParError(3)),"");
		// h_vav_leg->AddEntry((TObject*)0,Form("amplitude = %1.1f #pm %1.1f",fit_vec_v[i]->GetParameter(4),fit_vec_v[i]->GetParError(4)),"");
		// h_vav_leg->AddEntry(ln_vav_max,Form("MPV = %1.2f",vav_max),"l");

		// // h_vav_leg->Draw();



		// print results
		if (i<n_ch-1)
		{printf("ch%d\n",i );}
		else if(i==n_ch-1)
		{printf("sum\n");}
		
		printf("Gauss Fit:   mean: %1.2f ± %1.2f | red. chi2 = %1.2f \n",peak[i],peak_err[i],rchi2_g);
		printf("Landau Fit:   MPV: %1.2f ± %1.2f | red. chi2 = %1.2f \n",lan_mpv_c,lan_mpv_err,lan_rchi2);
		printf("Vavilov Fit:  MPV: %1.2f        | red. chi2 = %1.2f \n",peak_v[i],rchi2_v[i]);
    	// printf("Vavilov Fit: MPV: %1.2f         | red. chi2 = %1.2f \n",vav_max,vav_rchi2);

	} // ---> end loop over channels + sum

	/***** 
	__ MPV SUM ______________________________
	*****/

	double peak_sum = 0, peak_err_sum = 0;	
	if (target_wom=="WOM_A" || target_wom=="WOM_C"){
		for (int i = 0; i < 7; ++i){
			peak_sum += peak[i];
			peak_err_sum += peak_err[i]*peak_err[i];
		}

	}
	if (target_wom=="WOM_B" || target_wom=="WOM_D"){
		for (int i = 0; i < 8; ++i){
			peak_sum += peak[i];
			peak_err_sum += peak_err[i]*peak_err[i];
		}
	}
	peak_err_sum = sqrt(peak_err_sum);


	/***** 
	__ EXPORT ______________________________
	*****/

	printf("%s %d %d \n",particle.c_str(),run_nr,pos);

	if (store_result)
	{
		FILE *mpv_list, *mpv_err_list;
		mpv_list = fopen(Form("%s%s",out_dir.c_str(),out_values_list_name.c_str()),"a");
		mpv_err_list =fopen(Form("%s%s",out_dir.c_str(),out_values_err_list_name.c_str()),"a");

		//   // print date to file
		// time_t now;
		// time(&now); 
		// fprintf(mpv_list,"\nlight yield analysis - %s\n",ctime(&now));
		// fprintf(mpv_err_list,"\nlight yield analysis - %s\n",ctime(&now));

		fprintf(mpv_list, "%s ",target_wom.c_str() );
		fprintf(mpv_err_list, "%s ",target_wom.c_str() );
		if (target_wom=="WOM_A" || target_wom=="WOM_C"){
			for (int i = 0; i < 7; ++i){fprintf(mpv_list, "%f ",peak[i] );}
			fprintf(mpv_list, "0.000000 %f %f ",peak[7],peak_sum);
			for (int i = 0; i < 7; ++i){fprintf(mpv_err_list, "%f ",peak_err[i] );}
			fprintf(mpv_err_list, "0.000000 %f %f ",peak_err[7],peak_err_sum);
		}
		if (target_wom=="WOM_B" || target_wom=="WOM_D"){
			for (int i = 0; i < 8; ++i){fprintf(mpv_list, "%f ",peak[i] );}
			fprintf(mpv_list, "%f %f ",peak[8], peak_sum);
			for (int i = 0; i < 8; ++i){fprintf(mpv_err_list, "%f ",peak_err[i] );}
			fprintf(mpv_err_list, "%f %f ",peak_err[8], peak_err_sum);
		}
		fprintf(mpv_list, "%s %d %d %d\n",particle.c_str(),run_nr,pos,scan_nr);
		fprintf(mpv_err_list, "%s %d %d %d\n",particle.c_str(),run_nr,pos,scan_nr);


		gErrorIgnoreLevel = kError; // suppress root terminal output 
		C1->SaveAs(Form("%s%s",out_dir.c_str(),out_pdf_name.c_str()));
		gErrorIgnoreLevel = kUnset;

	}

	//if (is_interactive){ROOTapp->Run();}
	

  	return 0;
}