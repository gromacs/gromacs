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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum {
  eptAtom, eptNucleus, eptShell, eptBond, eptVSite, eptNR
};
/* The particle type */
 
enum {
  egcTC,    egcENER,   egcACC, egcFREEZE, 
  egcUser1, egcUser2,  egcVCM, egcXTC,
  egcORFIT,
  egcNR 
};

typedef struct {
  real 		m,q;		/* Mass and charge			*/
  real 		mB,qB;		/* Mass and charge for Free Energy calc */
  unsigned short type;		/* Atom type				*/
  unsigned short typeB;		/* Atom type for Free Energy calc	*/
  int           ptype;		/* Particle type			*/
  int 		resnr;		/* Residue number			*/
  unsigned char grpnr[egcNR];   /* Group numbers			*/
  unsigned char chain;          /* chain identifier                     */
} t_atom;

typedef struct {
  int  type;                    /* PDB record name                      */
  int  atomnr;                  /* PDB atom number                      */
  char altloc;                  /* Alternate location indicator         */
  char pdbresnr[6];             /* PDB res number                       */
  real occup;                   /* Occupancy                            */
  real bfac;                    /* B-factor                             */
  bool bAnisotropic;            /* (an)isotropic switch                 */
  int  uij[6];                  /* Anisotropic B-factor                 */
} t_pdbinfo;

typedef struct {
  int  nr;			/* Number of different groups		*/
  int  *nm_ind;                 /* Index in the group names             */
} t_grps;

typedef struct {
  int           nr;             /* Nr of atoms                          */
  t_atom	*atom;		/* Array of atoms (dim: nr)		*/
				/* The following entries will not 	*/
				/* allways be used (nres==0)	 	*/
  char		***atomname;	/* Array of pointers to atom name	*/
				/* use: (*(atomname[i]))		*/
  char		***atomtype;	/* Array of pointers to atom types	*/
				/* use: (*(atomtype[i]))		*/
  char		***atomtypeB;	/* Array of pointers to B atom types	*/
				/* use: (*(atomtypeB[i]))		*/
  int		nres;		/* Nr of residue names			*/
  char		***resname; 	/* Array of pointers to residue names 	*/
				/* use: (*(resname[i]))	       	*/
  int           ngrpname;        /* Number of groupnames                 */
  char          ***grpname;	/* Names of the groups		        */
  t_block	excl;		/* Exclusions		       	*/
  t_grps        grps[egcNR];    /* Groups of things                     */
  t_pdbinfo     *pdbinfo;       /* PDB Information, such as aniso. Bfac */
} t_atoms;

typedef struct {
  int           nr;              /* number of atomtypes                     */
  real         *radius;         /* GBSA radius for each atomtype        */
  real         *vol;            /* GBSA efective volume for each atomtype   */
  real         *surftens;       /* implicit solvent surftens for each atomtype */
} t_atomtypes;


#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))
