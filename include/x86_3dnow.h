/*
 * $Id$
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
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _x86_3dnow_h
#define _x86_3dnow_h

static char *SRCID_x86_3dnow_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if (defined USE_X86_ASM && !defined DOUBLE)

void check3dnow();
void vecinvsqrt_3dnow(float in[],float out[],int n);
void vecrecip_3dnow(float in[],float out[],int n);


void inl0100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[]);
void inl0110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[], int nsatoms[]);
void inl0300_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[],float tabscale,float VFtab[]);
void inl0310_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[],float tabscale,float VFtab[], int nsatoms[]);
void inl1000_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1010_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel, float Vc[],
		   int nsatoms[]);
void inl1020_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1030_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl1110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   int nsatoms[]);
void inl1120_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl1130_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl3000_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3010_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[], int nsatoms[]);
void inl3020_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3030_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[], int nsatoms[]);
void inl3120_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3130_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3300_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);
void inl3310_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[], int nsatoms[]);
void inl3320_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);
void inl3330_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);

 
#endif
#endif

