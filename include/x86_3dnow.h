/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _x86_3dnow_h
#define _x86_3dnow_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if (defined USE_X86_SSE_AND_3DNOW && !defined DOUBLE)

void check3dnow();
void vecinvsqrt_3dnow(float in[],float out[],int n);
void vecrecip_3dnow(float in[],float out[],int n);


void inl0100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const int type[],int ntype,const float nbfp[],
		   float Vnb[]);
void inl0110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const int type[],int ntype,const float nbfp[],
		   float Vnb[], int nsatoms[]);
void inl0300_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const int type[],int ntype,const float nbfp[],
		   float Vnb[],float tabscale,const float VFtab[]);
void inl0310_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const int type[],int ntype,const float nbfp[],
		   float Vnb[],float tabscale,const float VFtab[], int nsatoms[]);
void inl1000_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   const float shiftvec[],float fshift[],const int gid[],const float pos[],
                   float faction[],const float charge[],float facel,float Vc[]);
void inl1010_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   const float shiftvec[],float fshift[],const int gid[],const float pos[],
                   float faction[],const float charge[],float facel, float Vc[],
		   int nsatoms[]);
void inl1020_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   const float shiftvec[],float fshift[],const int gid[],const float pos[],
                   float faction[],const float charge[],float facel,float Vc[]);
void inl1030_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   const float shiftvec[],float fshift[],const int gid[],const float pos[],
                   float faction[],const float charge[],float facel,float Vc[]);
void inl1100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void inl1110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   int nsatoms[]);
void inl1120_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void inl1130_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void inl3000_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void inl3010_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[],const int nsatoms[]);
void inl3020_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void inl3030_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void inl3100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void inl3110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[],const int nsatoms[]);
void inl3120_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void inl3130_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void inl3300_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);
void inl3310_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[],const int nsatoms[]);
void inl3320_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);
void inl3330_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           const float shiftvec[],float fshift[],const int gid[],const float pos[],
	           float faction[],const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);


void mcinl0100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const int type[],int ntype,const float nbfp[],
		   float Vnb[]);
void mcinl0110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const int type[],int ntype,const float nbfp[],
		   float Vnb[],const int nsatoms[]);
void mcinl0300_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const int type[],int ntype,const float nbfp[],
		   float Vnb[],float tabscale,const float VFtab[]);
void mcinl0310_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const int type[],int ntype,const float nbfp[],
		   float Vnb[],float tabscale,const float VFtab[],const int nsatoms[]);
void mcinl1000_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   float shiftvec[],int gid[],const float pos[],
                   const float charge[],float facel,float Vc[]);
void mcinl1010_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   float shiftvec[],int gid[],const float pos[],
                   const float charge[],float facel, float Vc[],
		   const int nsatoms[]);
void mcinl1020_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   float shiftvec[],int gid[],const float pos[],
                   const float charge[],float facel,float Vc[]);
void mcinl1030_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
                   float shiftvec[],int gid[],const float pos[],
                   const float charge[],float facel,float Vc[]);
void mcinl1100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void mcinl1110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   const int nsatoms[]);
void mcinl1120_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void mcinl1130_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[]);
void mcinl3000_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void mcinl3010_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[],const int nsatoms[]);
void mcinl3020_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void mcinl3030_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   float tabscale,const float VFtab[]);
void mcinl3100_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void mcinl3110_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[],const int nsatoms[]);
void mcinl3120_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void mcinl3130_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale, const float VFtab[]);
void mcinl3300_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);
void mcinl3310_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[],const int nsatoms[]);
void mcinl3320_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);
void mcinl3330_3dnow(int nri,const int iinr[],const int jindex[],const int jjnr[],const int shift[],
	           float shiftvec[],int gid[],const float pos[],
	           const float charge[],float facel,float Vc[],
		   const int type[],int ntype,const float nbfp[],float Vnb[],
		   float tabscale,const float VFtab[]);

  
#endif
#endif

