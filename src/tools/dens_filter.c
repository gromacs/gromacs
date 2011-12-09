/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: densorder.c,v 0.9
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/* dens_filter.c
 * Routines for Filters and convolutions
 */

#include <math.h>
#include "typedefs.h"
#include "dens_filter.h"
#include "smalloc.h"
#include "vec.h"

#ifdef GMX_DOUBLE
#define EXP(x) (exp(x))
#else
#define EXP(x) (expf(x))
#endif

/* Modulo function assuming y>0 but with arbitrary INTEGER x */
static int MOD(int x, int y){
    if (x<0) x=x+y;
    return (mod(x,y));
}


gmx_bool convolution(int dataSize,real *x, int kernelSize, real* kernel)
{
    int i, j, k;
    real *out;
    snew(out,dataSize);
    /* check validity of params */
    if(!x || !kernel) return FALSE;
    if(dataSize <=0 || kernelSize <= 0) return FALSE;

    /* start convolution from out[kernelSize-1] to out[dataSize-1] (last) */
    for(i = kernelSize-1; i < dataSize; ++i)
    {
        for(j = i, k = 0; k < kernelSize; --j, ++k)
            out[i] += x[j] * kernel[k];
    }

    /* convolution from out[0] to out[kernelSize-2] */
    for(i = 0; i < kernelSize - 1; ++i)
    {
        for(j = i, k = 0; j >= 0; --j, ++k)
            out[i] += x[j] * kernel[k];
    }

    for(i=0;i<dataSize;i++) x[i]=out[i];
    sfree(out);	
    return TRUE;
}

/* Assuming kernel is shorter than x */

gmx_bool periodic_convolution(int datasize, real *x, int kernelsize, 
                              real *kernel)
{
    int i,j,k,nj;
    real *filtered;

    if (!x || !kernel) return FALSE;
    if (kernelsize<=0|| datasize<=0|| kernelsize > datasize) return FALSE;

    snew(filtered,datasize);
    
    for(i=0;(i<datasize); i++){
        for(j=0; (j<kernelsize); j++){
            filtered[i] += kernel[j]*x[MOD((i-j),datasize)];
        }
    }
    for(i=0; i<datasize; i++) x[i]=filtered[i];
    sfree(filtered);

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


