#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <iostream>

using namespace std;
TH1F *sumHistList[32]; //Declare it here to not get lost in the pointer parameter passing
int numberOfEntries;

int numberOfSumHistograms = 9;
void removeSubstring(char *s, const char *toremove)
{
	while (s = strstr(s, toremove))
		memmove(s, s + strlen(toremove), 1 + strlen(s + strlen(toremove)));
}

void addToSumHistogram(int counter, string token)
{
	string filePath = token.substr(0, token.size() - 2);
	TFile *file = new TFile(filePath.c_str());

	if (counter == 0)
	{
		//Get Number of Sum Histogramms
		TTree *tree;
		file->GetObject("T", tree);
		tree->SetBranchAddress("nCh", &numberOfSumHistograms);
		tree->GetEntry(1);

		for (int i = 0; i < numberOfSumHistograms; ++i)
		{
			string name = Form("hChSum_%d", (i));
			file->GetObject(name.c_str(), sumHistList[i]);
		}
	}
	else
	{
		for (int i = 0; i < numberOfSumHistograms; ++i)
		{
			TH1F *tempHist;
			string name = Form("hChSum_%d", (i));
			file->GetObject(name.c_str(), tempHist);
			 sumHistList[i]->Add(tempHist);
			 numberOfEntries += 1;
			//  cout << "number of entries " << numberOfEntries << endl;
		}
	}
}

int main(int argc, char const *argv[])
{
	//Argument: ./runs/Fast//22_muon6_pos4/0/out.root/T||./runs/Fast//22_muon6_pos4/1/out.root/T||./runs/Fast//22_muon6_pos4/2/out.root/T||
	//std::cout<<argv[0]<<"  ::::::::::: "<<argv[1]<<"::::::::"<<argc <<std::endl;
	TString inFileList;
	TString inDataFolder;
	TString outFile;
	TChain *chain = new TChain("T");
	std::string s = argv[1];
	TString outputFolder = argv[2];
	TString runName = argv[3];
	std::string delimiter = "||";

	std::cout << "MERGER STARTED--------------------------------" << std::endl;

	size_t pos = 0;
	std::string token;

	int counter = 0;
	while ((pos = s.find(delimiter)) != std::string::npos)
	{
		token = s.substr(0, pos);

		TString tok = token;

		if (tok.Length() > 0)
		{
			std::cout << "MERGER ADDED: " << token << std::endl;
			chain->Add(tok);

			addToSumHistogram(counter, token);

			counter++;
		}

		s.erase(0, pos + delimiter.length());
	}
	TString tok2 = s;
	std::cout << "MERGER ADDED END: " << s << std::endl;
	chain->Add(tok2);

	addToSumHistogram(counter, s);

	TString outputFile = outputFolder + "/" + runName + ".root";
	chain->Merge(outputFile);

	TFile *file = TFile::Open(outputFile, "UPDATE");
	
	
		for (int i = 0; i < numberOfSumHistograms; ++i)
		{
			if(sumHistList[i]->GetEntries()>0)
			sumHistList[i]->Write(Form("hChSum_%d",(i)), TObject::kSingleKey);
		}
	
	//file->ls();

	/* TList *l = new TList();
   TH1F *h1 = new TH1F("h1","h1",100.,0.,1.);
   TH1F *h2 = new TH1F("h2","h2",100.,0.,1.);
   TH1F *h3 = new TH1F("h3","h3",100.,0.,1.);
   l->Add(h1);
   l->Add(h2);
   l->Add(h3);
   TFile *f = new TFile("histlist.root","RECREATE");
   l->Write("histlist", TObject::kSingleKey);
   f->ls();
*/

	std::cout << "MERGER FINISHED-> Number of Events in the final File: " << chain->GetEntries() << std::endl;
	std::cout << "MERGER FINISHED-> Number of Events in the final SumHistogram[0]: " << (sumHistList[0]->GetEntries()) / 1024 << std::endl;
	// std::cout << "MERGER FINISHED-> Number of Events in the final SumHistogram[0] with counter " << numberOfEntries << std::endl;

	return 0;
}