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

#ifndef	_pdbio_h
#define	_pdbio_h

#ifdef HAVE_IDENT
#ident	"@(#) pdbio.h 1.12 7/28/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"
#include "symtab.h"

enum { epdbATOM, epdbHETATM, epdbNR };

typedef struct {
  int  pdbtp;
  int  atomnr;
  int  resnr;
  char atomnm[12];
  char resnm[12];
  char pdbresnr[12];
  char chain;
  rvec x;
  real bfac,dummy;
  real m,q;
  int  type;
} t_pdbatom;

extern void pdb_use_ter(bool bSet);
/* set read_pdbatoms to read upto 'TER' of 'ENDMDL' (default, bSet=FALSE) */

extern int read_pdbatoms(FILE *in,t_pdbatom **pdbaptr,matrix box,bool bFilterH);
/* Read a pdb file and extract ATOM and HETATM fields.
 * The pdbaptr will point to an array of pdb atoms on return.
 * Function returns number of atoms found.
 * atomnm and resnama still contain the spaces from the
 * pdb file. 
 * Read a box from the CRYST1 line, return 0 box when no CRYST1 is found.
 * It is possible to filter out hydrogen atoms (set bFilterH)
 */

extern void get_pdb_coordnum(char *infile,int *natoms);
/* Read a pdb file and count the ATOM and HETATM fields. */

extern bool is_hydrogen(char *nm);
/* Return whether atom nm is a hydrogen */

extern void pdba_trimnames(int natom,t_pdbatom pdba[]);
/* Remove leading and trailing spaces from all names */

extern void print_pdbatoms(FILE *out,int natom,t_pdbatom pdba[],matrix box);
/* Print the pdbatoms in proper format */

extern void renumber_pdb(int natom,t_pdbatom pdba[]);
/* Renumber residues starting from 0, and atoms starting from 0 */

extern int pdbasearch_atom(char *name,int resnr,int natom,t_pdbatom pdba[]);
/* Return the atom number of atom name in residue resnr, or -1 when not 
 * found.
 */

extern void pdb2atoms(int natom,t_pdbatom pdba[],t_atoms *atoms,rvec **x,
		      t_symtab *symtab);
/* Convert pdbatoms to atoms, allocates memory for x and components of
 * atoms. Assumes symtab is open and does not close symtab
 * Call renumber_pdb at start.
 */

extern t_pdbatom *atoms2pdba(t_atoms *atoms,rvec x[]);
/* Convert atoms to pdbatom */

#endif	/* _pdbio_h */
