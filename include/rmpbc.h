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
 * GROtesk MACabre and Sinister
 */

#ifndef _rmpbc_h
#define _rmpbc_h

static char *SRCID_rmpbc_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) rmpbc.h 1.6 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
	
extern void rm_pbc(t_idef *idef,int natoms,matrix box,rvec x[],rvec x_s[]);
/* Remove periodic boundary conditions */

#endif	/* _rmpbc_h */
