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

#ifndef _mdatoms_h
#define _mdatoms_h

static char *SRCID_mdatoms_h = "$Id$";

#include "typedefs.h"

extern t_mdatoms *atoms2md(FILE *fp,t_atoms *atoms,ivec nFreeze[],
			   bool bLD,bool bPert,bool bFree);
/* This routine copies the atoms->atom struct into a t_mdatoms struct
 * and then frees the atoms->atom struct if bFree is set.
 */

extern void md2atoms(t_mdatoms *md,t_atoms *atoms,bool bFree);
/* And vice versa */
#endif
