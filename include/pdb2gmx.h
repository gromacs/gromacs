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

#ifndef	_pdb2gmx_h
#define	_pdb2gmx_h

#ifdef HAVE_IDENT
#ident	"@(#) pdb2gmx.h 1.7 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "pdbio.h"

/* BONDS */
typedef struct {
  char 	*ai;		/* Atom i (may start with '-' to indicat prev residue*/
  char	*aj;		/* is bonded to atom j				*/
  real c[MAXFORCEPARAM];
} t_rbond;

typedef struct {
  char    *resname;	/* Residue name					*/
  int      nb;		/* Number of bonds			        */
  t_rbond *rbond;	/* The atom names of the bonded atoms	        */
} t_resbond;

typedef struct {
  char *ai;
  char *aj;
  char *ak;
  real c[MAXFORCEPARAM];
} t_rang;

/* ANGLES */
typedef struct {
  char    *resname;	/* Residue name					*/
  int      na;		/* Number of angles			        */
  t_rang  *rang; 	/* The atom names of the atoms in the angles    */
} t_resang;

/* RESIDUES */
typedef struct {
  char   *resname;
  int    natom;
  t_atom *atom;
  char   ***atomname;
  int    *cgnr;
} t_restp;

/* IMPROPERS */
typedef struct {
  char *ai[MAXATOMLIST];
  real c[MAXFORCEPARAM];
} t_idih;

typedef struct {
  char   *resname;
  int    nidih;
  t_idih *idih;
} t_idihres;

/* Structures for the h-database */
typedef struct {
  int 		nh;		/* Number of h-atoms		*/
  int 		tp;		/* Type of attachment (1..5)	*/
  char 		*na[4];		/* Control atoms i,j,k,l	*/
} t_add_block;

typedef struct {
  char 		*resname;	/* Residue name			*/
  int		n_add;		/* Number of add blocks		*/
  t_add_block 	*ab;		/* Array of add blocks		*/
} t_addh;

/* Block to hack terminal residues */
typedef struct {
  char        *bname;		/* Name of replace block		*/
  int         nreplace;		/* Number of atoms to replace		*/
  char        **nm_repl;	/* Names of atoms to be replaced 	*/
  char	      **new_nm;		/* Names of atoms after replacement	*/
  t_atom      *repl_by;		/* Atom data of replacing atoms		*/
  int         nadd;		/* Number of atoms to add		*/
  t_add_block *ab;		/* Addition data for add blocks		*/
  t_atom      *adder;		/* Atom data for adding atoms		*/
  char        **add_nm;		/* Names for added atoms		*/
  int         nidih;		/* Number of impropers to be added	*/
  t_idih      *idih;		/* Improper list			*/
  int         ndel;		/* Number of atoms to delete		*/
  char        **nm_del;		/* Names of atoms to delete		*/
} t_terblock;

/*extern real min_dist(int ncys,real **d,int *ii,int *jj,bool bCheckID);*/
/* Find the lowest number in a distance matrix d */

extern real distance(rvec a,rvec b);
/* calcluate distance between atoms without taking PBC into account */

extern void set_histp(int natom,t_pdbatom pdba[],real angle,real distance);
/* calculate HIStidine protonation state */


#endif	/* _pdb2gmx_h */
