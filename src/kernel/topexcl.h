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
 * S  C  A  M  O  R  G
 */

#ifndef _topexcl_h
#define _topexcl_h

static char *SRCID_topexcl_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) topexcl.h 1.11 11/23/92"
#endif /* HAVE_IDENT */

#include "topio.h"

typedef struct {
  int nr;		/* nr atoms (0 <= i < nr) (atoms->nr)	      	*/
  int nrex;		/* with nrex lists of neighbours		*/
			/* respectively containing zeroth, first	*/
			/* second etc. neigbours (0 <= nre < nrex)	*/
  int **nrexcl;		/* with (0 <= nrx < nrexcl[i][nre]) neigbours 	*/ 
			/* per list stored in one 2d array of lists	*/ 
  int ***a;		/* like this: a[i][nre][nrx]			*/
} t_nextnb;

extern void init_nnb(t_nextnb *nnb, int nr, int nrex);
/* Initiate the arrays for nnb (see above) */

extern void done_nnb(t_nextnb *nnb);
/* Cleanup the nnb struct */

extern void print_nnb(t_nextnb *nnb, char *s);
/* Print the nnb struct */

extern void gen_nnb(t_nextnb *nnb,t_params plist[]);
/* Generate a t_nextnb structure from bond information. 
 * With the structure you can either generate exclusions
 * or generate angles and dihedrals. The structure must be
 * initiated using init_nnb.
 */

extern void generate_excl (int nrexcl, int nratoms,
			   t_params plist[],t_block *excl);
/* Generate an exclusion block from bonds and constraints in
 * plist.
 */
#endif	/* _topexcl_h */
