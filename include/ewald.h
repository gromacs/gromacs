/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _ewald_h
#define _ewald_h

static char *SRCID_ewald_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "complex.h"
#include "fftgrid.h"

extern real calc_ewaldcoeff(real rc,real dtol);
/* Determines the Ewald parameter, both for Ewald and PME */

extern real do_ewald(FILE *log,       bool bVerbose,
		     t_inputrec *ir,
		    rvec x[],        rvec f[],
		    real charge[],   rvec box,
		    real phi[],      t_commrec *cr,
		     t_nsborder *nsb, t_nrnb *nrnb,
		     matrix lrvir, real ewaldcoeff);
/* Do an Ewald calculation for the long range electrostatics. */
 
#endif


