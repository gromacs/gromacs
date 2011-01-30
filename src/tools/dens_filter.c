/* dens_filter.c
 * Routines for Filters and convolutions
 */

#include <math.h>
#include "typedefs.h"
#include "dens_filter.h"
#include "smalloc.h"

#ifdef GMX_DOUBLE
#define EXP(x) (exp(x))
#else
#define EXP(x) (expf(x))
#endif



gmx_bool convolve1D(real* in, real* out, int dataSize, real* kernel, int kernelSize)
{
    int i, j, k;

    // check validity of params
    if(!in || !out || !kernel) return FALSE;
    if(dataSize <=0 || kernelSize <= 0) return FALSE;

    // start convolution from out[kernelSize-1] to out[dataSize-1] (last)
    for(i = kernelSize-1; i < dataSize; ++i)
    {
        out[i] = 0;                             // init to 0 before accumulate

        for(j = i, k = 0; k < kernelSize; --j, ++k)
            out[i] += in[j] * kernel[k];
    }

    // convolution from out[0] to out[kernelSize-2]
    for(i = 0; i < kernelSize - 1; ++i)
    {
        out[i] = 0;                             // init to 0 before sum

        for(j = i, k = 0; j >= 0; --j, ++k)
            out[i] += in[j] * kernel[k];
    }

    return TRUE;
}

/* returns discrete gaussian kernel of size n in *out, where n=2k+1=3,5,7,9,11 and k=1,2,3 is the order
 * NO checks are performed 
 */
void gausskernel(real *out, int n, real var){
	int i,j=0,k;
	real arg,tot=0;
	k=n/2;
	
	for(i=-k; i<=k; i++){
		arg=(i*i)/(2*var);
		tot+=out[j++]=EXP(-arg);
	}
	for(i=0;i<j;i++){
		out[i]/=tot;
	}
}


