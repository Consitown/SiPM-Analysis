/// Functions specific for the cosmics box, a liquid scintillator prototype read out with a SiPM array coupled to a wavelength-shifting optical module WOM 
#ifndef _CosmicsBox
#define _CosmicsBox

#include "ReadRun.h"

class CosmicsBox : public ReadRun {
public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read 
	CosmicsBox(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) {}
	
	void Print_Phi_ew(vector<int>, vector<float>, vector<int>, float = 1, float = 1, float = 100, float = 140, int = 400, bool = true, bool = false);
};
#endif