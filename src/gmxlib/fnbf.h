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
 
#ifndef _fnbf_h
#define _fnbf_h

#include "typedefs.h"

extern void do_fnbf(FILE *log,int ftype,t_forcerec *fr,
		    rvec x[],rvec f[],t_mdatoms *md,
		    real egnb[],real egcoul[],rvec box_size,
		    t_nrnb *nrnb,real lambda,real *dvdlambda);

extern void fdo_flr(FILE *log,int nri,atom_id i_atoms[],int shift,
		    int njcg,atom_id jcg[],atom_id index[],atom_id acg[],
		    rvec x[],real egcoul[],
		    t_mdatoms *md,int ngener,t_forcerec *fr);

#endif
