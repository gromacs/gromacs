#ifndef _dens_filter_h
#define _dens_filter_h
	
#include "types/simple.h"

gmx_bool convolution(int dataSize, real* in, int kernelSize, real* kernel);
gmx_bool periodic_convolution(int dsize, real *in, int ksize, real* kernel);
void gausskernel(real *out, int size, real var);

#endif
