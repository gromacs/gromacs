/*
 * $Id$
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _typedefs_h
#define _typedefs_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

#define STRLEN 4096
#define NOTSET -12345

#include <sys/types.h>
#include <sysstuff.h>
#include <types/simple.h>
#include <types/enums.h>
#include <types/block.h>
#include <types/symtab.h>
#include <types/idef.h>
#include <types/atoms.h>
#include <types/trx.h>
#include <types/topology.h>
#include <types/energy.h>
#include <types/inputrec.h>
#include <types/ishift.h>
#include <types/graph.h>
#include <types/nrnb.h>
#include <types/nblist.h>
#include <types/nsgrid.h>
#include <types/commrec.h>
#include <types/forcerec.h>
#include <types/fcdata.h>
#include <types/mdatom.h>
#include <types/pbc.h>
#include <types/ifunc.h>
#include <types/filenm.h>
#include <types/group.h>
#include <types/state.h>
#include <types/matrix.h>
#include <types/edsams.h>

extern void set_over_alloc(bool set);
  /* Turns over allocation on or off, default is off */
  
extern int over_alloc(int n);
  /* Returns n when over allocation is off.
   * Returns 1.1*n + 10 when over allocation in on.
   * This is to avoid frequent reallocation
   * during domain decomposition in mdrun.
   */

/* Functions to initiate and delete structures *
 * These functions are defined in gmxlib/typedefs.c 
 */
extern void init_block(t_block *block);
extern void init_atom (t_atoms *at);
extern void init_top (t_topology *top);
extern void init_inputrec(t_inputrec *ir);
extern void init_gtc_state(t_state *state,int ngtc);
extern void init_state(t_state *state,int natoms,int ngtc);
extern void done_block(t_block *block);
extern void done_atom (t_atoms *at);
extern void done_top(t_topology *top);
extern void done_inputrec(t_inputrec *ir);
extern void done_state(t_state *state);

extern void stupid_fill(t_block *grp, int natom,bool bOneIndexGroup);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * If bOneIndexGroup, then all atoms are  lumped in one index group,
 * otherwise there is one atom per index entry
 */

extern void init_t_atoms(t_atoms *atoms, int natoms, bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */  

extern void free_t_atoms(t_atoms *atoms);
/* free all the arrays and set the nr and nres to 0 */

#ifdef CPLUSPLUS
}
#endif


#endif	/* _typedefs_h */
