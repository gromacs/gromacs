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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _typedefs_h
#define _typedefs_h


#define STRLEN 4096
#define NOTSET -12345
#define BIG_STRLEN 1048576

#include <sys/types.h>
#include "sysstuff.h"
#include "types/simple.h"
#include "types/enums.h"
#include "types/block.h"
#include "types/symtab.h"
#include "types/idef.h"
#include "types/atoms.h"
#include "types/trx.h"
#include "types/topology.h"
#include "types/energy.h"
#include "types/inputrec.h"
#include "types/ishift.h"
#include "types/graph.h"
#include "types/nrnb.h"
#include "types/nblist.h"
#include "types/nsgrid.h"
#include "types/commrec.h"
#include "types/forcerec.h"
#include "types/fcdata.h"
#include "types/mdatom.h"
#include "types/pbc.h"
#include "types/ifunc.h"
#include "types/filenm.h"
#include "types/group.h"
#include "types/state.h"
#include "types/shellfc.h"
#include "types/constr.h"
#include "types/matrix.h"
#include "types/oenv.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Memory (re)allocation can be VERY slow, especially with some
 * MPI libraries that replace the standard malloc and realloc calls.
 * To avoid slow memory allocation we use over_alloc to set the memory
 * allocation size for large data blocks. Since this scales the size
 * with a factor, we use log(n) realloc calls instead of n.
 * This can reduce allocation times from minutes to seconds.
 */
/* This factor leads to 4 realloc calls to double the array size */
#define OVER_ALLOC_FAC 1.19

void set_over_alloc_dd(gmx_bool set);
  /* Turns over allocation for variable size atoms/cg/top arrays on or off,
   * default is off.
   */
  
int over_alloc_dd(int n);
  /* Returns n when domain decomposition over allocation is off.
   * Returns OVER_ALLOC_FAC*n + 100 when over allocation in on.
   * This is to avoid frequent reallocation
   * during domain decomposition in mdrun.
   */

/* Over allocation for small data types: int, real etc. */
#define over_alloc_small(n) (OVER_ALLOC_FAC*(n) + 8000)

/* Over allocation for large data types: complex structs */
#define over_alloc_large(n) (OVER_ALLOC_FAC*(n) + 1000)

int gmx_large_int_to_int(gmx_large_int_t step,const char *warn);
/* Convert a gmx_large_int_t value to int.
 * If warn!=NULL a warning message will be written
 * to stderr when step does not fit in an int,
 * the first line is:
 * "WARNING during %s:", where warn is printed in %s.
 */

#define STEPSTRSIZE 22

char *gmx_step_str(gmx_large_int_t i,char *buf);
/* Prints a gmx_large_int_t value in buf and returns the pointer to buf.
 * buf should be large enough to contain i: STEPSTRSIZE (22) chars.
 * When multiple gmx_large_int_t values are printed in the same printf call,
 * be sure to call gmx_step_str with different buffers.
 */

/* Functions to initiate and delete structures *
 * These functions are defined in gmxlib/typedefs.c 
 */
void init_block(t_block *block);
void init_blocka(t_blocka *block);
void init_atom (t_atoms *at);
void init_mtop(gmx_mtop_t *mtop);
void init_top (t_topology *top);
void init_inputrec(t_inputrec *ir);
void init_energyhistory(energyhistory_t * enerhist);
void done_energyhistory(energyhistory_t * enerhist);
void init_gtc_state(t_state *state,int ngtc, int nnhpres, int nhchainlength);
void init_state(t_state *state,int natoms,int ngtc, int nnhpres, int nhchainlength);

void copy_blocka(const t_blocka *src,t_blocka *dest);

void done_block(t_block *block);
void done_blocka(t_blocka *block);
void done_atom (t_atoms *at);
void done_moltype(gmx_moltype_t *molt);
void done_molblock(gmx_molblock_t *molb);
void done_mtop(gmx_mtop_t *mtop,gmx_bool bDoneSymtab);
void done_top(t_topology *top);
void done_inputrec(t_inputrec *ir);
void done_state(t_state *state);

void set_box_rel(t_inputrec *ir,t_state *state);
/* Set state->box_rel used in mdrun to preserve the box shape */

void preserve_box_shape(t_inputrec *ir,matrix box_rel,matrix b);
/* Preserve the box shape, b can be box or boxv */

void stupid_fill_block(t_block *grp, int natom,gmx_bool bOneIndexGroup);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * If bOneIndexGroup, then all atoms are  lumped in one index group,
 * otherwise there is one atom per index entry
 */

void stupid_fill_blocka(t_blocka *grp, int natom);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * There is one atom per index entry
 */

void init_t_atoms(t_atoms *atoms, int natoms, gmx_bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */  

t_atoms *copy_t_atoms(t_atoms *src);
/* copy an atoms struct from src to a new one */
  
void add_t_atoms(t_atoms *atoms,int natom_extra,int nres_extra);
/* allocate extra space for more atoms and or residues */
 
void t_atoms_set_resinfo(t_atoms *atoms,int atom_ind,t_symtab *symtab,
				const char *resname,int resnr,unsigned char ic,
				int chainnum, char chainid);
/* Set the residue name, number, insertion code and chain identifier
 * of atom index atom_ind.
 */

void free_t_atoms(t_atoms *atoms,gmx_bool bFreeNames);
/* Free all the arrays and set the nr and nres to 0.
 * bFreeNames tells if to free the atom and residue name strings,
 * don't free them if they still need to be used in e.g. the topology struct.
 */

t_atoms *mtop2atoms(gmx_mtop_t *mtop);
/* generate a t_atoms struct for the system from gmx_mtop_t */ 

real max_cutoff(real cutoff1,real cutoff2);
/* Returns the maximum of the cut-off's, taking into account that 0=inf. */

#ifdef __cplusplus
}
#endif


#endif	/* _typedefs_h */
