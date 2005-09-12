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

#ifndef _pme_h
#define _pme_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "gmxcomplex.h"
#include "fftgrid.h"

typedef real *splinevec[DIM];

extern real do_pme(FILE *log,       bool bVerbose,
		   t_inputrec *ir,
		   rvec x[],        rvec f[],
		   real chargeA[],  real chargeB[],
		   matrix box,      t_commrec *cr,
		   t_nsborder *nsb, t_nrnb *nrnb,
		   matrix lrvir,    real ewaldcoeff,
		   bool bFreeEnergy,
		   real lambda,     real *dvdlambda,
		   bool bGatherOnly);
    
/* Do a PME calculation for the long range electrostatics. 
 * If bGatherOnly is set, the energy from the last computation will be used, and 
 * the forces will be interpolated at the new positions. No new solving is done then.
 */

extern void sum_qgrid(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,
		      int pme_order,bool bForward);

extern void sum_qgrid_dd(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,
                         int pme_order, bool bForward);

/* the sum_qgrid routine with domain decomposition (dd) yields better 
 * performance in parallel runs
 */	       

extern t_fftgrid *init_pme(FILE *log,t_commrec *cr,
			   int nkx,int nky,int nkz,int pme_order,int homenr,
			   bool bFreeEnergy,bool bOptFFT,int ewald_geometry);

/* Routine for spreading something on a grid. Can be misused for non-PME
 * related things. init_pme must be called before this guy.
 */
typedef struct {
  int atom,index;
} t_pmeatom;

extern void spread_on_grid(FILE *logfile,   t_commrec *cr,
			   t_fftgrid *grid, int homenr,
			   int pme_order,   rvec x[],
			   real charge[],   matrix box,      
			   bool bGatherOnly,bool bHaveSplines);
		
/*
 * the following three routines are for PME/PP node splitting:
 */

extern void send_coordinates(t_nsborder *nsb, rvec x[], 
                             real chargeA[],  real chargeB[], 
			     bool bLastTime, bool bFreeEnergy,
			     t_commrec *cr);

/* send the particle coordinates (and charges) from the PP to the PME nodes,
 * so that they can do a PME step */

extern void receive_lrforces(t_commrec *cr, t_nsborder *nsb, 
                      rvec f[]     , matrix vir, 
		      real *energy,  int step);

/* PP nodes receive the long range forces from the PME nodes */


extern void do_pmeonly(FILE *logfile, 
                       t_inputrec *ir,    t_commrec *cr,
	  	       matrix box, 
		       t_nsborder *nsb,   t_nrnb *mynrnb,
		                          real ewaldcoeff,
		       real lambda,
		       bool bGatherOnly   );

/* Called on the nodes that do PME exclusively (as slaves) 
 */
		   
#endif
