/// Class containing experimental functions
#ifndef _Experimental
#define _Experimental

#include "ReadRun.h"
#include <TRandom3.h>

class Experimental : public virtual ReadRun {
public:
	/// @brief Initializer will call initializer of ReadRun class
	/// @param no_of_bin_files_to_read 
	Experimental(int no_of_bin_files_to_read) : ReadRun(no_of_bin_files_to_read) { }
	
	void RebinAll(int = 2, float = 0, unsigned long = 0);
	void DerivativeAll();
};
#endif