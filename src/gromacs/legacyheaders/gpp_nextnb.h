/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _gpp_nextnb_h
#define _gpp_nextnb_h

#include "grompp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int nr;		/* nr atoms (0 <= i < nr) (atoms->nr)	      	*/
  int nrex;		/* with nrex lists of neighbours		*/
			/* respectively containing zeroth, first	*/
			/* second etc. neigbours (0 <= nre < nrex)	*/
  int **nrexcl;		/* with (0 <= nrx < nrexcl[i][nre]) neigbours 	*/ 
			/* per list stored in one 2d array of lists	*/ 
  int ***a;		/* like this: a[i][nre][nrx]			*/
} t_nextnb;

void init_nnb(t_nextnb *nnb, int nr, int nrex);
/* Initiate the arrays for nnb (see above) */

void done_nnb(t_nextnb *nnb);
/* Cleanup the nnb struct */

#ifdef DEBUG_NNB
#define print_nnb(nnb, s) __print_nnb(nnb, s)
void print_nnb(t_nextnb *nnb, char *s);
/* Print the nnb struct */
#else
#define print_nnb(nnb, s)
#endif

void gen_nnb(t_nextnb *nnb,t_params plist[]);
/* Generate a t_nextnb structure from bond information. 
 * With the structure you can either generate exclusions
 * or generate angles and dihedrals. The structure must be
 * initiated using init_nnb.
 */

void nnb2excl (t_nextnb *nnb, t_blocka *excl);
/* generate exclusions from nnb */

void generate_excl (int nrexcl, int nratoms,
			   t_params plist[],t_nextnb *nnb,t_blocka *excl);
/* Generate an exclusion block from bonds and constraints in
 * plist.
 */

#ifdef __cplusplus
}
#endif

#endif	/* _gpp_nextnb_h */
