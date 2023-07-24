/// Class contains method for ploting the Fourier transforms of waveforms. For analysis of noise/signal frequencies, potential filters, or input data for machine learning.
#ifndef _FFT_WF
#define _FFT_WF

#include "ReadRun.h"

class FFT_WF : public ReadRun {
public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read 
	FFT_WF(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) {}
	
	// print FFT
	void PrintFFTWF(int = 1, float = 0., float = 0., int = 1);
};
#endif