/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
#ifndef _pppmcore_h
#define _pppmcore_h

#include <stdio.h>
#include "typedefs.h"
#include "complex.h"

extern real do_pppm(FILE *log,       bool bVerbose,
		    bool bGenerGhat, char *ghatfn,
		    t_inputrec *ir,  int natoms,
		    rvec x[],        rvec f[],
		    real charge[],   rvec box,
		    real phi[],      t_commrec *cr,
		    t_nrnb *nrnb);
/* Do a PPPM calculation for the long range electrostatics.
 */
		    
extern real do_ewald(FILE *log,       t_inputrec *ir,
		     int natoms,      rvec x[],rvec f[],
		     real charge[],   rvec box,
		     real phi[],      t_commrec *cr);
/* Do an Ewald summation on a fixed grid as given in inputrec.
 * The spread function is David's function, rather than a gaussian.
 */
 
#endif


