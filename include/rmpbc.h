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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _rmpbc_h
#define _rmpbc_h

static char *SRCID_rmpbc_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) rmpbc.h 1.6 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
	
extern void rm_pbc(t_idef *idef,int natoms,matrix box,rvec x[],rvec x_s[]);
/* Remove periodic boundary conditions.
 * natoms is the size of x and x_s and can be smaller than the number 
 * of atoms in idef, but should only contain whole molecules
 */

#endif	/* _rmpbc_h */
