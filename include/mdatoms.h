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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _mdatoms_h
#define _mdatoms_h

static char *SRCID_mdatoms_h = "$Id$";

#include "typedefs.h"

extern t_mdatoms *atoms2md(t_atoms *atoms,bool bPert);
/* This routine copies the atoms->atom struct into a t_mdatoms struct
 * and then frees the atoms->atom struct.
 */

extern void md2atoms(t_mdatoms *md,t_atoms *atoms);
/* And vice versa */
#endif
