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

#ifndef _genhydro_h
#define _genhydro_h

static char *SRCID_genhydro_h = "$Id$";

#include "pdbio.h"
#include "hackblock.h"

extern int add_h(t_atoms **pdbaptr, rvec *xptr[], int nah, t_hackblock ah[],
		 t_hackblock *ntdb, t_hackblock *ctdb, 
		 int nterpairs, int *rN, int *rC);
/* Generate hydrogen atoms and N and C terminal patches.
 * ntdb and ctdb may be NULL, no replacement will be done then.
 * int nterpairs is the number of termini pairs in the molecule
 * rN is the residue number of the N-terminus,
 * rC is the residue number of the C-terminus
 * return the New total number of atoms 
 */
 
extern void protonate(t_atoms **atoms,rvec **x);
/* Protonate protein molecule */

extern void deprotonate(t_atoms *atoms,rvec *x);
/* Deprotonate any molecule: all atoms whose name begins with H will be 
 * removed 
 */

#endif

