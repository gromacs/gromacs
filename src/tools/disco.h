/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_disco_h = "$Id$";

#define _correct_h

#include <stdio.h>
#include "typedefs.h"

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

extern bool shake_coords(FILE *log,bool bVerbose,int nstruct,
			 int natom,rvec xref[],rvec x[],int *seed,
			 matrix box,t_correct *c,int *niter);
/* Do the actual shaking. Return TRUE when converged */

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



