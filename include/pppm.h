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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _pppm_h
#define _pppm_h

static char *SRCID_pppm_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "complex.h"
#include "fftgrid.h"

extern void init_pppm(FILE *log,t_commrec *cr,t_nsborder *nsb,
		      bool bVerbose,bool bOld,
		      rvec box,char *ghatfn,t_inputrec *ir);
/* Setup stuff for PPPM. 
 * Either reads a ghat function from file (when the file exists)
 * or generate a ghat function from scratch.
 */

extern real do_pppm(FILE *log,       bool bVerbose,
		    rvec x[],        rvec f[],
		    real charge[],   rvec box,
		    real phi[],      t_commrec *cr,
		    t_nsborder *nsb, t_nrnb *nrnb);
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


