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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef _toppush_h
#define _toppush_h

static char *SRCID_toppush_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) toppush.h 1.31 9/30/97"
#endif /* HAVE_IDENT */

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
			      
extern void push_at (t_symtab *symtab, t_atomtype *at, char *line,int nb_funct);

extern void push_bt(directive d,t_params bt[], int nral, 
		    t_atomtype *at, char *line);

extern void push_nbt(directive d,t_params nbt[],t_atomtype *atype,
		     char *plines,int nb_funct);

extern void push_atom(t_symtab   *symtab, 
		      t_block    *cgs,
		      t_atoms    *at,
		      t_atomtype *atype,
		      char       *line);

extern void push_bondnow (t_params *bond, t_param *b);

extern void push_bond(directive d,t_params bondtype[],t_params bond[],
		      t_atoms *at,char *line);

extern void push_mol(int nrmols,t_molinfo mols[],char *pline,
		     int *whichmol,int *nrcopies);

extern void push_molt(t_symtab *symtab,t_molinfo *newmol,char *line);

extern void init_block2(t_block2 *b2, int natom);
	
extern void done_block2(t_block2 *b2);

extern void push_excl(char *line, t_block2 *b2);

extern void merge_excl(t_block *excl, t_block2 *b2);

extern void b_to_b2(t_block *b, t_block2 *b2);

extern void b2_to_b(t_block2 *b2, t_block *b);

#endif	/* _toppush_h */
