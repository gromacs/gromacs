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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#define _correct_h

#include <stdio.h>
#include "typedefs.h"
#include "network.h"

#ifdef debug_gmx
#undef debug_gmx
#endif

#ifdef HAVE_MPI
#define debug_gmx() do { FILE *fp=debug ? debug : (stdlog ? stdlog : stderr);\
fprintf(fp,"NODEID=%d, %s  %d\n",gmx_cpu_id(),__FILE__,__LINE__); fflush(fp); } while (0)
#else
#define debug_gmx()
#endif

enum { edcNONBOND, edcBOND, edcDISRE, edcNR }; 

extern char *edc_names[edcNR+1];

typedef struct {
  int  ai,aj,cons_type;
  real lb,ub,wi;
} t_dist;

typedef struct {
  int ai,aj,ak,al,am;
} t_quadruple;

/* Garbage bin with flags and data for the correct routine */
typedef struct {
  int         maxnit,nbcheck,nstprint,nstranlist,ngrow;
  bool        bExplicit,bChiral,bPep,bDump,bLowerOnly,bRanlistFirst;
  bool        bCubic,bBox,bCenter;
  real        lodev;
  int         maxdist,ndist,npep,nimp;
  int         *ip,*tag;
  t_dist      *d;
  t_quadruple *pepbond,*imp;
  real        *omega,*idih,*weight;
  bool        *bViol;
} t_correct;

extern t_correct *init_corr(int maxnit,int nstprint,int nbcheck,int nstranlist,
			    int ngrow,bool bExplicit,bool bChiral,bool bPep,
			    bool bDump,real lowdev,bool bLowerOnly,
			    bool bRanlistFirst,bool bCubic,bool bBox,bool bCenter);
/* Initiate the data structure and set some of the parameters */

extern void init_corr2(t_correct *c,int natom);
/* Make tags (indices) for optimization of shaking and 
 * initiate a number of other arrays.
 */

extern int shake_coords(FILE *log,bool bVerbose,int nstruct,
			int natom,rvec xref[],rvec x[],int *seed,
			matrix box,t_correct *c,int *niter);
/* Do the actual shaking. Returns the final number of violations */

extern int quick_check(FILE *log,int natom,rvec x[],matrix box,t_correct *c);
/* Check a structure once for the number of violations (return value) */

extern real *read_weights(char *fn,int natom);
/* Read the weights from the occupancy field in the pdb file */  

extern void define_peptide_bonds(FILE *log,t_atoms *atoms,t_correct *c);
/* Fill the peptide bond structure */

extern void define_impropers(FILE *log,t_atoms *atoms,t_correct *c);
/* Fill the improper structure */

extern void read_dist(FILE *log,char *fn,int natom,t_correct *c);
/* Read distances from a dist.dat file produced by cdist */

extern void pr_conv_stat(FILE *fp,int ntry,int nconv,double tnit);
/* Statistics */

extern void pr_corr(FILE *log,t_correct *c);
/* Print parameters from the corr structure */

extern void pr_distances(FILE *fp,t_correct *c);
/* Print the distances in a cdist compatible file format */

extern void center_in_box(int natom,rvec xin[],matrix box,rvec xout[]);
/* Center the coordinates in the box. xin and xout may be the same array */

extern void check_dist(FILE *log,t_correct *c);
/* Check internal consistency of the distances */

extern void check_final(FILE *log,t_correct *c,rvec x[]);
/* Check structure for distances */

extern void rand_coords(int natom,rvec x[],rvec xref[],real weight[],
			bool bCenter,rvec xcenter[],rvec box,int *seed);
/* Generate random coordinates */

extern void rand_box(bool bUserBox,
		     matrix box,rvec boxsize,int nres,bool bCubic,int *seed);
/* Generate random box */

extern void disco_slave(t_commrec *cr,FILE *log);
/* Slave process for parallel disco */

extern void disco_master(t_commrec *cr,FILE *log,char *outfn,char *keepfn,t_correct *c,
			 bool bVerbose,t_atoms *atoms,
			 rvec xref[],rvec xcenter[],
			 int nstruct,int *seed,
			 bool bFit,int nfit,atom_id fit_ind[],
			 bool bPrintViol,char *violfn,rvec boxsize);
/* Master process for parallel disco */



