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

#ifndef _pme_h
#define _pme_h

static char *SRCID_pme_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "complex.h"
#include "fftgrid.h"

typedef real *splinevec[DIM];

extern real do_pme(FILE *log,       bool bVerbose,
		   t_inputrec *ir,
		   rvec x[],        rvec f[],
		   real charge[],   matrix box,
		   t_commrec *cr,
		   t_nsborder *nsb, t_nrnb *nrnb,
		   matrix lrvir,real ewaldcoeff,
		   bool bGatherOnly);
    
/* Do a PME calculation for the long range electrostatics. 
 * If bGatherOnly is set, the energy from the last computation will be used, and 
 * the forces will be interpolated at the new positions. No new solving is done then.
 */

extern void sum_qgrid(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,bool bForward);

extern void init_pme(FILE *log,t_commrec *cr,
		     int nkx,int nky,int nkz,int pme_order,int homenr,
		     bool bOptFFT);

/* Routine for spreading something on a grid. Can be misused for non-PME
 * related things. init_pme must be called before this guy.
 */
extern t_fftgrid *spread_on_grid(FILE *logfile,   int homenr,
				 int pme_order,   rvec x[],
				 real charge[],   matrix box,
				 bool bGatherOnly);

#endif
