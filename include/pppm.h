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

#ifndef _pppm_h
#define _pppm_h

static char *SRCID_pppm_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "complex.h"
#include "fftgrid.h"

extern real do_pppm(FILE *log,       bool bVerbose,
		    bool bGenerGhat, char *ghatfn,
		    t_inputrec *ir,  int natoms,
		    rvec x[],        rvec f[],
		    real charge[],   rvec box,
		    real phi[],      t_commrec *cr,
		    t_nrnb *nrnb,    bool bNew);
/* Do a PPPM calculation for the long range electrostatics.
 */
 
extern real do_opt_pppm(FILE *log,       bool bVerbose,
			t_inputrec *ir,  int natoms,
			rvec x[],        rvec f[],
			real charge[],   rvec box,
			real phi[],      t_commrec *cr,
			t_nrnb *nrnb,    rvec beta,
			t_fftgrid *grid, bool bNew);
/* Do a PPPM setup (generate grid etc.) and a calculation as well 
 * the grid should be initiated beforehand.
 */

		    
extern real do_ewald(FILE *log,       t_inputrec *ir,
		     int natoms,      rvec x[],rvec f[],
		     real charge[],   rvec box,
		     real phi[],      t_commrec *cr,
		     bool bNew);
/* Do an Ewald summation on a fixed grid as given in inputrec.
 * The spread function is David's function, rather than a gaussian.
 */
 
extern real gather_f(FILE *log,bool bVerbose,
		     int natoms,rvec x[],rvec f[],real charge[],rvec box,
		     real pot[],t_fftgrid *grid,rvec beta,t_nrnb *nrnb);
/* Gather the forces and potential from a grid */
 
#endif


