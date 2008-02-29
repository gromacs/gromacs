/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _grompp_h
#define _grompp_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "macros.h"

#define MAXSLEN 128

typedef struct {
  bool bSet;                    /* Has this combination been set        */
  real c[4];                    /* The non-bonded parameters            */
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
  t_param	*nb;		/* Nonbonded force default params	*/
  int           *bondatomtype;  /* The bond_atomtype for each atomtype  */
  real          *radius;        /* Radius for GBSA stuff                */
  real          *vol;           /* Effective volume for GBSA            */
  real          *surftens;      /* Surface tension with water, for GBSA */
  int           *atomnumber;    /* Atomic number, used for QM/MM        */
} t_atomtype;

typedef struct {
  int           nr;             /* Number of bond_atomtypes		*/
  char          ***atomname;    /* Names of the bond_atomtypes		*/
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
  d_moleculetype,
  d_atoms,
  d_vsites2,
  d_vsites3,
  d_vsites4,
  d_bonds,
  d_exclusions,
  d_pairs,
  d_angles,
  d_dihedrals,
  d_constraints,
  d_settles,
  d_polarization,
  d_water_polarization,
  d_thole_polarization,
  d_system,
  d_molecules,
  d_position_restraints,
  d_angle_restraints,
  d_angle_restraints_z,
  d_distance_restraints,
  d_orientation_restraints,
  d_dihedral_restraints,
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
  /* All the directives above can not appear after moleculetype */
  "moleculetype",
  "atoms",
  "virtual_sites2",
  "virtual_sites3",
  "virtual_sites4",
  "bonds",
  "exclusions",
  "pairs",
  "angles",
  "dihedrals",
  "constraints",
  "settles",
  "polarization",
  "water_polarization",
  "thole_polarization",
  "system",
  "molecules",
  "position_restraints",
  "angle_restraints",
  "angle_restraints_z",
  "distance_restraints",
  "orientation_restraints",
  "dihedral_restraints",
  "invalid"
  };

extern void convert_harmonics(int nrmols,t_molinfo mols[],t_atomtype *atype);

#endif	/* _grompp_h */








