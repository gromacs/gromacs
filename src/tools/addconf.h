/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_addconf_h = "$Id$";

#include "typedefs.h"

extern void add_conf(t_atoms *atoms_1,rvec *x_1,rvec *v_1,real *r_1,
		     int NTB,matrix box_1,
		     t_atoms *atoms_2,rvec *x_2,rvec *v_2,real *r_2,
		     bool bVerbose,int maxmol);
/* Add two conformations together */
