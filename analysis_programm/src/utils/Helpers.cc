#include "Helpers.h"

/// @brief Helper. Creates a list of .bin data files in data folder to be read in
/// @param dirname Directory
/// @param ext File extension
/// @return String of line separated file names
string Helpers::ListFiles(const char* dirname, const char* ext) {
	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	
	TIter next(files);
	TObjString* objString;
	while ((objString = (TObjString*)next())) {
		const char* fileName = objString->GetString().Data();
		// Check if the filename or extension contains ".bin"
		if (strstr(fileName, ext) != nullptr) ss << fileName << "\n";
	}
	if (ss.str().empty()) throw runtime_error("Error: No .bin files found in " + static_cast<string>(dirname) + "\n");
	return ss.str();
}

/// @brief Print progress bar for a loop in steps of 10 percent
/// @param index Current loop index
/// @param length Length of loop
void Helpers::PrintProgressBar(int index, int length) {
	if ((index + 1) % (length / 10) == 0) {
		float progress = 10 * (index + 1) / length;
		cout << "\r[" << string(progress * 5, '#')
			<< std::string(50 - progress * 5, ' ')
			<< "] " << int(progress * 10) << "%";
		cout.flush();
		if (index + 1 == length) cout << endl << endl;
	}
}

/// @brief Translate a random number into a useful root color https://root.cern.ch/doc/master/classTColor.html
/// @param i Index of your plotting loop that is to be translated into a useful ROOT color index
/// @return ROOT color index
int Helpers::rcolor(unsigned int i) {
	const int nclrs = 17;
	int rclrs[nclrs] = { 1, 2, 3, 4, 6, 8, 9, 13, 20, 28, 30, 34, 38, 40, 31, 46, 49 };
	return rclrs[i - static_cast<int>(floor(i / nclrs)) * nclrs];
}

/// @brief Set consistent x-axis and y-axis range for all TH1 histograms on a canvas. 
/// 
/// Will only change the axes where ```min != max```.
/// 
/// @param c Canvas Canvas
/// @param x_range_min X-axis minimum
/// @param x_range_max X-axis maximum
/// @param y_range_min Y-axis minimum
/// @param y_range_max Y-axis maximum
void Helpers::SetRangeCanvas(TCanvas*& c, double x_range_min, double x_range_max, double y_range_min, double y_range_max) {

	// Lambda function to set the axis ranges
	auto setAxisRanges = [&](TH1* his) {
		if (x_range_min != x_range_max) his->GetXaxis()->SetRangeUser(x_range_min, x_range_max);
		if (y_range_min != y_range_max) his->GetYaxis()->SetRangeUser(y_range_min, y_range_max);
		};

	// Get the list of pads on the canvas and loop over pads
	TList* pads = c->GetListOfPrimitives();
	TIter nextPad(pads);
	TObject* object;
	while ((object = nextPad())) {
		if (object->InheritsFrom(TPad::Class())) {
			TPad* pad = static_cast<TPad*>(object);
			// Set the axis ranges for all plots on the current pad
			pad->cd();
			TList* primitives = pad->GetListOfPrimitives();
			TIter nextPrimitive(primitives);
			TObject* primitive;
			while ((primitive = nextPrimitive())) {
				if (primitive->InheritsFrom(TH1::Class())) {
					setAxisRanges(static_cast<TH1*>(primitive));
				}
			}
		}
		else if (object->InheritsFrom(TH1::Class())) {
			setAxisRanges(static_cast<TH1*>(object));
		}
	}
	c->Modified();
	c->Update();
}

/// @brief Helper to split canvas according to the number of channels to be plotted
/// @param c Canvas to be split
void Helpers::SplitCanvas(TCanvas*& c, vector<int> active_channels, vector<int> plot_active_channels) {
	// cross check if user input exists in data
	vector<int> rmv;
	for (int i = 0; i < static_cast<int>(plot_active_channels.size()); i++) {
		if (find(active_channels.begin(), active_channels.end(), plot_active_channels[i]) == active_channels.end()) {
			cout << "\n\n\n ------------ WARNING ------------\n";
			cout << "YOUR SELECTED CHANNEL " << plot_active_channels[i] << " DOES NOT EXIST IN DATA\n";
			cout << "PLEASE CHANGE plot_active_channels\n\n\n";
			rmv.push_back(plot_active_channels[i]);
		}
	}

	for (int i = 0; i < static_cast<int>(rmv.size()); i++) {
		auto it = find(plot_active_channels.begin(), plot_active_channels.end(), rmv[i]);
		if (it != plot_active_channels.end()) plot_active_channels.erase(it);
	}

	if (plot_active_channels.empty()) {
		c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	}
	else if (static_cast<int>(plot_active_channels.size()) > 1) {
		c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
	}
}

/// @brief Get array of y values for a histogram
/// @param his TH1F histogram
/// @return Y values of waveform
double* Helpers::gety(TH1F* his) {
	double* yvals = new double[his->GetNbinsX()];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i + 1);
	}
	return yvals;
}

/// @brief Get truncated array of y values for a certain waveform
/// @param his Waveform histogram
/// @param start_at Truncate from index
/// @param end_at Truncate to index
/// @return Truncated Y values of waveform
double* Helpers::gety(TH1F* his, int start_at, int end_at) {
	if (start_at < 0 || start_at >= his->GetNbinsX() || end_at >= his->GetNbinsX() || end_at - start_at < 1) {
		cout << "\nError: Helpers::gety out of range" << endl;
		return 0;
	}
	const int n_bins_new = end_at - start_at;
	double* yvals = new double[n_bins_new];
	for (int i = start_at; i < end_at; i++) {
		yvals[i - start_at] = his->GetBinContent(i + 1);
	}
	return yvals;
}