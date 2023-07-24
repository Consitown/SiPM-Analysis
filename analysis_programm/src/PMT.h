/// Functions specific for PMTs
#ifndef _PMT
#define _PMT

#include "ReadRun.h"

class PMT : public virtual ReadRun {
private:
	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrumPMT_cnt;
	/// @brief Index for multiple executions of the same plotting function
	int PrintChargeSpectrumPMTthreshold_cnt;

public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read 
	PMT(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) {
		PrintChargeSpectrumPMT_cnt = 0;
		PrintChargeSpectrumPMTthreshold_cnt = 0;
	}
	
	void PrintChargeSpectrumPMT(float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	/// @brief Starting values of the fit parameters for PrintChargeSpectrumPMT()
	vector<float> PrintChargeSpectrumPMT_pars;
	void PrintChargeSpectrumPMTthreshold(float = 0, float = 0, float = 0, float = 300, int = 750, double = 4, bool = false);
};
#endif