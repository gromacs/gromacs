/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _pp2shift_h
#define _pp2shift_h

#include "typedefs.h"
	
/* must correspond with 'leg' g_chi.c:727 */
enum { edPhi=0, edPsi, edOmega, edChi1, edChi2, edChi3, edChi4, edChi5, edChi6, edMax };

enum { edPrintST=0,edPrintRO } ; 

#define NHISTO 360
#define NONCHI 3
#define MAXCHI edMax-NONCHI
#define NROT 4  /* number of rotamers: 1=g(-), 2=t, 3=g(+), 0=other */ 

typedef struct {
  int minO,minC,H,N,C,O,Cn[MAXCHI+3];
} t_dihatms; /* Cn[0]=N, Cn[1]=Ca, Cn[2]=Cb etc. */

typedef struct {
  char name[12];
  int  resnr;
  int  index;       /* Index for amino acids (histograms) */
  int  j0[edMax];   /* Index in dih array (phi angle is first...) */
  t_dihatms  atm;
  int  b[edMax];
  int  ntr[edMax];
  real S2[edMax];
  real rot_occ[edMax][NROT];

} t_dlist;

extern void do_pp2shifts(FILE *fp,int nframes,
			 int nlist,t_dlist dlist[],real **dih);

extern bool has_dihedral(int Dih,t_dlist *dl);

extern t_dlist *mk_dlist(FILE *log, 
			 t_atoms *atoms, int *nlist,
			 bool bPhi, bool bPsi, bool bChi, int maxchi,
			 int r0,int naa,char **aa);
			 
extern void pr_dlist(FILE *fp,int nl,t_dlist dl[],real dt,  int printtype,
bool bPhi, bool bPsi,bool bChi,bool bOmega, int maxchi);

extern int pr_trans(FILE *fp,int nl,t_dlist dl[],real dt,int Xi);

extern void mk_chi_lookup (int **lookup, int maxchi, real **dih, 
			   int nlist, t_dlist dlist[]) ; 

extern void get_chi_product_traj (real **dih,int nframes,int nangles, 
			   int nlist,int maxchi, t_dlist dlist[], real time[], 
			   int **lookup,int *xity,bool bRb,bool bNormalize,
			   real core_frac); 

extern void print_one (char *base,char *name,char *title, char *ylabel,
		      int nf,real time[],real data[]); 

#endif




