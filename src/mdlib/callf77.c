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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */

/* Fortran versions of shake and settle */

#ifdef USE_FORTRAN
extern void F77_FUNC(forsettle,FORSETTLE)(int *nshake,int owptr[],real
					  b4[],real after[],real
					  *dOH,real *dHH,real *mO,real
					  *mH,int *error);
extern void F77_FUNC(forshake,FORSHAKE)(atom_id iatom[],int *ncon,
					int *nit, int *maxnit,
					real dist2[],real xp[],
					real rij[],real m2[],
					real *omega,
					real invmass[],real tt[],
					real lambda[],int *error);
extern void F77_FUNC(forlincs,FORLINCS)(real *x,real *xp,int *nc,
					int *bla1,int *bla2,int *blnr,
					int *blbnb,real *bllen,
					real *blc,real *blcc,real *blm,
					int *nit,int *nrec,real *invmass,
					real *r,real *temp1,real *temp2,
					real *temp3,real *wangle,
					int *warn,real *lambda);
#endif 

void fsettle(int *nshake,int owptr[],real b4[],real after[],real *dOH,real *dHH,real *mO,real *mH,int *error)
{
#ifdef USE_FORTRAN
  F77_FUNC(forsettle,FORSETTLE) (nshake,owptr,b4,after,dOH,dHH,mO,mH,error);
#else
  gmx_fatal(FARGS,"fsettle called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

void fshake(atom_id iatom[],int *ncon,int *nit,int *maxnit,real dist2[],real xp[],real rij[],real m2[],real *omega,real invmass[],real tt[],real lambda[],int *error)
{
#ifdef USE_FORTRAN
  F77_FUNC(forshake,FORSHAKE)(iatom,ncon,nit,maxnit,dist2,xp,rij,m2,omega,
			      invmass,tt,lambda,error);
#else
  gmx_fatal(FARGS,"fshake called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

/* LINCS */

void flincs(real *x,real *xp,int *nc,int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,real *blc,real *blcc,real *blm,int *nit,int *nrec,real *invmass,real *r,real *temp1,real *temp2,real *temp3,real *wangle,int *warn,real *lambda)
{
#ifdef USE_FORTRAN
  F77_FUNC(forlincs,FORLINCS)(x,xp,nc,bla1,bla2,blnr,blbnb,bllen,blc,blcc,
  	blm,nit,nrec,invmass,r,temp1,temp2,temp3,wangle,warn,lambda);
#else
  gmx_fatal(FARGS,"flincs called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}
