/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_TOPOLOGY_ATOMS_H
#define GMX_TOPOLOGY_ATOMS_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_symtab;

enum {
    eptAtom, eptNucleus, eptShell, eptBond, eptVSite, eptNR
};
/* The particle type */

typedef struct t_atom
{
    real           m, q;        /* Mass and charge                      */
    real           mB, qB;      /* Mass and charge for Free Energy calc */
    unsigned short type;        /* Atom type                            */
    unsigned short typeB;       /* Atom type for Free Energy calc       */
    int            ptype;       /* Particle type                        */
    int            resind;      /* Index into resinfo (in t_atoms)      */
    int            atomnumber;  /* Atomic Number or NOTSET              */
    char           elem[4];     /* Element name                         */
} t_atom;

typedef struct t_resinfo
{
    char          **name;       /* Pointer to the residue name          */
    int             nr;         /* Residue number                       */
    unsigned char   ic;         /* Code for insertion of residues       */
    int             chainnum;   /* Iincremented at TER or new chain id  */
    char            chainid;    /* Chain identifier written/read to pdb */
    char          **rtp;        /* rtp building block name (optional)   */
} t_resinfo;

typedef struct t_pdbinfo
{
    int      type;              /* PDB record name                      */
    int      atomnr;            /* PDB atom number                      */
    char     altloc;            /* Alternate location indicator         */
    char     atomnm[6];         /* True atom name including leading spaces */
    real     occup;             /* Occupancy                            */
    real     bfac;              /* B-factor                             */
    gmx_bool bAnisotropic;      /* (an)isotropic switch                 */
    int      uij[6];            /* Anisotropic B-factor                 */
} t_pdbinfo;

typedef struct t_grps
{
    int   nr;                   /* Number of different groups           */
    int  *nm_ind;               /* Index in the group names             */
} t_grps;

typedef struct t_atoms
{
    int            nr;          /* Nr of atoms                          */
    t_atom        *atom;        /* Array of atoms (dim: nr)             */
                                /* The following entries will not       */
                                /* always be used (nres==0)             */
    char          ***atomname;  /* Array of pointers to atom name       */
                                /* use: (*(atomname[i]))                */
    char          ***atomtype;  /* Array of pointers to atom types      */
                                /* use: (*(atomtype[i]))                */
    char          ***atomtypeB; /* Array of pointers to B atom types    */
                                /* use: (*(atomtypeB[i]))               */
    int              nres;      /* The number of resinfo entries        */
    t_resinfo       *resinfo;   /* Array of residue names and numbers   */
    t_pdbinfo       *pdbinfo;   /* PDB Information, such as aniso. Bfac */
} t_atoms;

typedef struct t_atomtypes
{
    int           nr;           /* number of atomtypes                          */
    real         *radius;       /* GBSA radius for each atomtype                */
    real         *vol;          /* GBSA efective volume for each atomtype       */
    real         *surftens;     /* implicit solvent surftens for each atomtype  */
    real         *gb_radius;    /* GB radius for each atom type                 */
    real         *S_hct;        /* Overlap factors for HCT/OBC GB models        */
    int          *atomnumber;   /* Atomic number, used for QM/MM                */
} t_atomtypes;

#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))

void init_atom(t_atoms *at);
void init_atomtypes(t_atomtypes *at);
void done_atom(t_atoms *at);
void done_atomtypes(t_atomtypes *at);

void init_t_atoms(t_atoms *atoms, int natoms, gmx_bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */

t_atoms *copy_t_atoms(t_atoms *src);
/* copy an atoms struct from src to a new one */

void add_t_atoms(t_atoms *atoms, int natom_extra, int nres_extra);
/* allocate extra space for more atoms and or residues */

void t_atoms_set_resinfo(t_atoms *atoms, int atom_ind, struct t_symtab *symtab,
                         const char *resname, int resnr, unsigned char ic,
                         int chainnum, char chainid);
/* Set the residue name, number, insertion code and chain identifier
 * of atom index atom_ind.
 */

void free_t_atoms(t_atoms *atoms, gmx_bool bFreeNames);
/* Free all the arrays and set the nr and nres to 0.
 * bFreeNames tells if to free the atom and residue name strings,
 * don't free them if they still need to be used in e.g. the topology struct.
 */

#ifdef __cplusplus
}
#endif

#endif
