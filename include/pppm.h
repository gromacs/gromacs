/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

#ifndef _pppm_h
#define _pppm_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "gmxcomplex.h"
#include "fftgrid.h"

extern void init_pppm(FILE *log,t_commrec *cr,t_nsborder *nsb,
		      bool bVerbose,bool bOld,
		      matrix box,char *ghatfn,t_inputrec *ir);
/* Setup stuff for PPPM. 
 * Either reads a ghat function from file (when the file exists)
 * or generate a ghat function from scratch.
 */

extern real do_pppm(FILE *log,       bool bVerbose,
		    rvec x[],        rvec f[],
		    real charge[],   rvec box,
		    real phi[],      t_commrec *cr,
		    t_nsborder *nsb, t_nrnb *nrnb,
		    int pme_order);
/* Do a PPPM calculation for the long range electrostatics. */
 
extern real do_opt_pppm(FILE *log,       bool bVerbose,
			t_inputrec *ir,  int natoms,
			rvec x[],        rvec f[],
			real charge[],   rvec box,
			real phi[],      t_commrec *cr,
			t_nrnb *nrnb,    rvec beta,
			t_fftgrid *grid, bool bOld);
/* Do a PPPM setup (generate grid etc.) and a calculation as well 
 * the grid should be initiated beforehand.
 */

extern void calc_invh(rvec box,int nx,int ny,int nz,rvec invh);
		    
extern void spread_q(FILE *log,bool bVerbose,
		     int start,int nr,
		     rvec x[],real charge[],rvec box,
		     t_fftgrid *grid,t_nrnb *nrnb);
 
#endif


