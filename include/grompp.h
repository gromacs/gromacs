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
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _grompp_h
#define _grompp_h

static char *SRCID_grompp_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) grompp.h 1.37 2/2/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "macros.h"

typedef struct {
  atom_id a[MAXATOMLIST];	/* The atom list (eg. bonds: particle	*/
				/* i = a[0] (AI), j = a[1] (AJ))	*/
  real 	  c[MAXFORCEPARAM];	/* Force parameters (eg. b0 = c[0])	*/
} t_param;

typedef struct {
  int		nr;		/* The number of bonds in this record 	*/
  t_param 	*param;		/* Array of parameters (dim: nr)	*/
} t_params;

typedef struct {
  char          **name;
  int		nrexcl;		/* Number of exclusions per atom	*/
  bool		excl_set;	/* Have exclusions been generated?	*/
  bool          bProcessed;     /* Has the mol been processed           */
  t_atoms       atoms;          /* Atoms                                */
  t_block       cgs;            /* Charge groups                        */
  t_block       mols;           /* Molecules                            */
  t_params      plist[F_NRE];   /* Parameters in old style              */
} t_molinfo;

typedef struct {
  int           nr;		/* The number of atomtypes		*/
  t_atom	*atom;		/* Array of atoms			*/
  char          ***atomname;	/* Names of the atomtypes		*/
  t_param	*nb;		/* Nonbonded force default params	*/
} t_atomtype;

typedef enum {
  d_defaults,
  d_atomtypes,
  d_bondtypes,
  d_pairtypes,
  d_angletypes,
  d_dihedraltypes,
  d_nonbond_params,
  d_blocktype,
  d_moleculetype,
  d_atoms,
  d_dum2,
  d_dum3,
  d_bonds,
  d_exclusions,
  d_pairs,
  d_angles,
  d_dihedrals,
  d_constraints,
  d_settles,
  d_system,
  d_molecules,
  d_position_restraints,
  d_distance_restraints,
  d_maxdir,
  d_invalid,
  d_none
} directive;

#endif	/* _grompp_h */
