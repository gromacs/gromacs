#ifndef _dens_filter_h
#define _dens_filter_h
	
#include "types/simple.h"

gmx_bool convolve1D(real* in, real* out, int dataSize, real* kernel, int kernelSize);
void gausskernel(real *out, int size, real var);

#endif
