/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

#ifndef _grompp_h
#define _grompp_h

static char *SRCID_grompp_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef HAVE_IDENT
#ident	"@(#) grompp.h 1.37 2/2/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "macros.h"

#define MAXSLEN 32

typedef struct {
  bool bSet;                    /* Has this combination been set        */
  real c[3];                    /* The non-bonded parameters            */
} t_nbparam;
/* The t_nbparam struct is used to temporary store the explicit
 * non-bonded parameter combinations, which will be copied to t_params.
 */

typedef struct {
  atom_id a[MAXATOMLIST];	/* The atom list (eg. bonds: particle	*/
				/* i = a[0] (AI), j = a[1] (AJ))	*/
  real 	  c[MAXFORCEPARAM];	/* Force parameters (eg. b0 = c[0])	*/
  char    s[MAXSLEN];           /* A string (instead of parameters),    *
				 * read from the .rtp file in pdb2gmx   */
} t_param;

typedef struct {
  int		nr;		/* The number of bonds in this record 	*/
  t_param 	*param;		/* Array of parameters (dim: nr)	*/
} t_params;

typedef struct {
  int           nr;             /* The number of exclusions             */
  atom_id       *e;             /* The excluded atoms                   */
} t_excls;

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
  int           *bondatomtype;  /* bond_atomtype for each atomtype      */
  t_param	*nb;		/* Nonbonded force default params	*/
} t_atomtype;

typedef struct {
  int           nr;             /* Number of bond_atomtypes             */
  char          ***atomname;    /* Names of the bond_atomtypes          */
} t_bond_atomtype;

typedef enum {
  d_defaults,
  d_atomtypes,
  d_bondtypes,
  d_constrainttypes,
  d_pairtypes,
  d_angletypes,
  d_dihedraltypes,
  d_nonbond_params,
  d_blocktype,
  d_moleculetype,
  d_atoms,
  d_dum2,
  d_dum3,
  d_dum4,
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
  d_angle_restraints,
  d_angle_restraints_z,
  d_distance_restraints,
  d_orientation_restraints,
  d_maxdir,
  d_invalid,
  d_none
} directive;

static char *ds[d_maxdir+1] = {
  "defaults",
  "atomtypes",
  "bondtypes",
  "constrainttypes",
  "pairtypes",
  "angletypes",
  "dihedraltypes",
  "nonbond_params",
  "blocktype",
  /* All the directives above can not appear after moleculetype */
  "moleculetype",
  "atoms",
  "dummies2",
  "dummies3",
  "dummies4",
  "bonds",
  "exclusions",
  "pairs",
  "angles",
  "dihedrals",
  "constraints",
  "settles",
  "system",
  "molecules",
  "position_restraints",
  "angle_restraints",
  "angle_restraints_z",
  "distance_restraints",
  "orientation_restraints",
  "invalid"
  };

extern void convert_harmonics(int nrmols,t_molinfo mols[],t_atomtype *atype);

#endif	/* _grompp_h */
