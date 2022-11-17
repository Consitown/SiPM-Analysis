void fit()  {
      TFile* f = new TFile("out.root","READ");
      TH1F* his;

      f->GetObject("integral_histo",his);
      his->Draw();
      
      gPad->SetLogy();

      auto two_gauss = new TF1("three gaussians", "gaus(0)+gaus(3)+gaus(6)", -15, 80); two_gauss->SetTitle("Sum of three gauss");
	two_gauss->SetParameters(1e4, 0., 0.21, 2e2, 20, 3, 10, 40, 10);
      //two_gauss->FixParameter(4,20);
      //two_gauss->FixParameter(7,40);
      TFitResultPtr fresults = his->Fit(two_gauss, "RSL");
      two_gauss->Draw("same");

      auto pedestal = new TF1("pedestal", "gaus", -15, 80); pedestal->SetTitle("pedestal");
	pedestal->SetParameters(fresults->Parameter(0), fresults->Parameter(1), fresults->Parameter(2));
	pedestal->SetLineColor(3);
	pedestal->Draw("same");

	auto pepeak = new TF1("pepeak", "gaus", -15, 80); pepeak->SetTitle("pepeak");
	pepeak->SetParameters(fresults->Parameter(3), fresults->Parameter(4), fresults->Parameter(5));
	pepeak->SetLineColor(4);
	pepeak->Draw("same");

	auto secpepeak = new TF1("secpepeak", "gaus", -15, 80); secpepeak->SetTitle("secpepeak");
	secpepeak->SetParameters(fresults->Parameter(6), fresults->Parameter(7), fresults->Parameter(8));
	secpepeak->SetLineColor(5);
	secpepeak->Draw("same");
}