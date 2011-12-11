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

#ifndef _toppush_h
#define _toppush_h

#include "typedefs.h"
#include "toputil.h"
#include "gpp_atomtype.h"
#include "gpp_bond_atomtype.h"
#include "warninp.h"

typedef struct {
  int     nr;	/* The number of entries in the list 			*/
  int     nra2; /* The total number of entries in a			*/
  atom_id *nra;	/* The number of entries in each a array (dim nr) 	*/
  atom_id **a;	/* The atom numbers (dim nr) the length of each element	*/
		/* i is nra[i]						*/
} t_block2;

extern void generate_nbparams(int comb,int funct,t_params plist[],
			      gpp_atomtype_t atype,
			      warninp_t wi);
			      
extern void push_at (t_symtab *symtab,gpp_atomtype_t at,
		     t_bond_atomtype bat,char *line,int nb_funct,
		     t_nbparam ***nbparam,t_nbparam ***pair,
		     warninp_t wi);

extern void push_bt(directive d,t_params bt[], int nral,
		    gpp_atomtype_t at,t_bond_atomtype bat,char *line,
		    warninp_t wi);

extern void push_dihedraltype(directive d,t_params bt[],
			      t_bond_atomtype bat,char *line,
			      warninp_t wi);

extern void push_cmaptype(directive d,t_params bt[], int nral, gpp_atomtype_t at,
			  t_bond_atomtype bat,char *line,
			  warninp_t wi);

extern void push_nbt(directive d,t_nbparam **nbt,gpp_atomtype_t atype,
		     char *plines,int nb_funct,
		     warninp_t wi);

void
push_gb_params(gpp_atomtype_t atype,
	       char       *line,
	       warninp_t wi);

extern void push_atom(t_symtab   *symtab, 
		      t_block    *cgs,
		      t_atoms    *at,
		      gpp_atomtype_t atype,
		      char       *line,
		      int        *lastcg,
		      warninp_t  wi);

extern void push_bond(directive d,t_params bondtype[],t_params bond[],
		      t_atoms *at,gpp_atomtype_t atype,char *line,
		      gmx_bool bBonded,gmx_bool bGenPairs,real fudgeQQ,
		      gmx_bool bZero,gmx_bool *bWarn_copy_A_B,
		      warninp_t wi);

extern void push_cmap(directive d, t_params bondtype[], t_params bond[],
		      t_atoms *at, gpp_atomtype_t atype, char *line,
		      gmx_bool *bWarn_copy_A_B,
		      warninp_t wi);

extern void push_vsitesn(directive d,t_params bondtype[],t_params bond[],
			 t_atoms *at,gpp_atomtype_t atype,char *line,
			 warninp_t wi);

extern void push_mol(int nrmols,t_molinfo mols[],char *pline,
		     int *whichmol,int *nrcopies,
		     warninp_t wi);

extern void push_molt(t_symtab *symtab,int *nmol,t_molinfo **mol,char *line,
		      warninp_t wi);

extern void init_block2(t_block2 *b2, int natom);
	
extern void done_block2(t_block2 *b2);

extern void push_excl(char *line, t_block2 *b2);

extern void merge_excl(t_blocka *excl, t_block2 *b2);

extern void b_to_b2(t_blocka *b, t_block2 *b2);

extern void b2_to_b(t_block2 *b2, t_blocka *b);

extern int add_atomtype_decoupled(t_symtab *symtab,gpp_atomtype_t at,
				  t_nbparam ***nbparam,t_nbparam ***pair);
/* Add an atom type with all parameters set to zero (no interactions).
 * Returns the atom type number.
 */

extern void convert_moltype_couple(t_molinfo *mol,int atomtype_decouple,
				   real fudgeQQ,
				   int couple_lam0,int couple_lam1,
				   gmx_bool bCoupleIntra,
				   int nb_funct,t_params *nbp);
/* Setup mol such that the B-state has no interaction with the rest
 * of the system, but full interaction with itself.
 */

#endif	/* _toppush_h */
