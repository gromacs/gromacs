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
static char *SRCID_addconf_h = "$Id$";

#include "typedefs.h"

extern void add_conf(t_atoms *atoms, rvec **x, real **r, bool bSrenew,
		     int NTB,matrix box,
		     t_atoms *atoms_solvt, rvec *x_solvt, real *r_solvt, 
		     bool bVerbose);
/* Add two conformations together, without generating overlap */

extern void orient_mol(t_atoms *atoms,char *indexnm,rvec x[],rvec *v);
/* Orient a molecule along its principal component axes. 
 * indexnm may be the name of an index file or null.
 */
