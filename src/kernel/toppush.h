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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _toppush_h
#define _toppush_h

#include "typedefs.h"
#include "toputil.h"

typedef struct {
  int     nr;	/* The number of entries in the list 			*/
  int     nra2; /* The total number of entries in a			*/
  atom_id *nra;	/* The number of entries in each a array (dim nr) 	*/
  atom_id **a;	/* The atom numbers (dim nr) the length of each element	*/
		/* i is nra[i]						*/
} t_block2;

extern void generate_nbparams(int comb,int funct,t_params plist[],
			      t_atomtype *atype);
			      
extern void push_at (t_symtab *symtab,t_atomtype *at,t_bond_atomtype *bat,char *line,int nb_funct,
		     t_nbparam ***nbparam,t_nbparam ***pair);

extern void push_bt(directive d,t_params bt[], int nral, 
		    char ***typenames, int ntypes, char *line);

extern void push_dihedraltype(directive d,t_params bt[],
			      char ***typenames, int ntypes,char *line);

extern void push_nbt(directive d,t_nbparam **nbt,t_atomtype *atype,
		     char *plines,int nb_funct);

extern void push_atom(t_symtab   *symtab, 
		      t_block    *cgs,
		      t_atoms    *at,
		      t_atomtype *atype,
		      char       *line,
		      int        *lastcg);

extern void push_bondnow (t_params *bond, t_param *b);

extern void push_bond(directive d,t_params bondtype[],t_params bond[],
		      t_atoms *at,t_atomtype *atype,char *line,
		      bool bBonded,bool bGenPairs,
		      bool bZero,bool *bWarn_copy_A_B);

extern void push_mol(int nrmols,t_molinfo mols[],char *pline,
		     int *whichmol,int *nrcopies);

extern void push_molt(t_symtab *symtab,int *nmol,t_molinfo **mol,char *line);

extern void init_block2(t_block2 *b2, int natom);
	
extern void done_block2(t_block2 *b2);

extern void push_excl(char *line, t_block2 *b2);

extern void merge_excl(t_block *excl, t_block2 *b2);

extern void b_to_b2(t_block *b, t_block2 *b2);

extern void b2_to_b(t_block2 *b2, t_block *b);

#endif	/* _toppush_h */
