/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */

#ifndef	_protein_h
#define	_protein_h

#ifdef HAVE_IDENT
#ident	"@(#) protein.h 1.30 9/30/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "grompp.h"
#include "topexcl.h"
#include "pdb2gmx.h"

typedef struct t_alist {
  t_atom         at;
  char           **name;
  rvec           x;
  struct t_alist *next;
} t_alist;

typedef struct {
  t_symtab   symtab;
  char       **header;
  int        nres;
  char       ***resname;
  int        natom;
  t_alist    *al;
} t_seq;

extern t_alist *init_al(t_atom *at,rvec x,char **name);

extern void set_at(t_atom *at,real m,real q,int type,int resnr);

extern int search_atom(char *type,int start,int natoms,char **atom[]);
/* Search an atom index starting from start. 
 * If the first character of type is a '-' searching will be done backward.
 * Returns the atom number on success, or -1 when failing.
 */

extern t_atomtype *read_atype();
/* Read an atomtype database from adb */

extern void gen_h(t_seq *seq,int nah,t_addh ah[]);
/* Add hydrogens to the sequence, using database ah */

#endif	/* _protein_h */
