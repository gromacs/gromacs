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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _pgutil_h
#define _pgutil_h

static char *SRCID_pgutil_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pgutil.h 1.14 9/30/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"

extern int search_atom(char *type,int start,
		       int natoms,t_atom at[],char **anm[]);
/* Search an atom in array of pointers to strings, starting from start
 * if type starts with '-' then searches backwards from start.
 */

extern void set_at(t_atom *at,real m,real q,int type,int resnr);

#endif	/* _pgutil_h */
