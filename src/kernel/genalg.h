/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_genalg_h = "$Id$";

#ifndef _genalg_h
#define _genalg_h
typedef struct {
  int  np;          /* Number of points           */
  int  atype;       /* Atom type                  */
  int  ptype;       /* Parameter type             */
  real rmin,rmax;   /* Minimum and maximum value  */
  real dr;
  real rval;        /* Current value              */
} t_range;

typedef struct {
  int  NP,D;
  int  strategy;
  int  seed;
  int  ipop,gen;        /* Current population member and generation */
  int  imin;            /* Member with lowest energy */
  real CR,FF;
  real **pold,**pnew;   /* Old and new populations */
  real *best,*bestit,*cost,*tmp,*rmsf,*energy;
} t_genalg;

enum { eseSIGMA, eseEPSILON, eseBHAMA, eseBHAMB, eseBHAMC, eseNR };

extern real value_rand(t_range *r,int *seed);

extern t_genalg *init_ga(char *infile,int D,t_range range[]);

extern void update_ga(FILE *fpout_ptr,t_range range[],t_genalg *ga);

extern bool print_ga(FILE *fp,t_genalg *ga,real rmsf,real energy,t_range range[],
		     real tol);

extern real cost(real rmsf,real energy);

#endif
