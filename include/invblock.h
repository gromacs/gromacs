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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef _invblock_h
#define _invblock_h

static char *SRCID_invblock_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) invblock.h 1.8 2/2/97"
#endif /* HAVE_IDENT */
#include <typedefs.h>

extern atom_id *make_invblock(t_block *block,int nr);
/* Inverse the block structure. nr is the maximum entry in the inversed
 * array, and therefore the dimension of the returned array
 */

#endif	/* _invblock_h */
