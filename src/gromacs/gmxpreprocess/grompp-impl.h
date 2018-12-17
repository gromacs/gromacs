/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GROMPP_IMPL_H
#define GMX_GMXPREPROCESS_GROMPP_IMPL_H

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_block;
struct t_blocka;

#define MAXSLEN 32

typedef struct {
    bool     bSet;          /* Has this combination been set        */
    real     c[4];          /* The non-bonded parameters            */
} t_nbparam;
/* The t_nbparam struct is used to temporary store the explicit
 * non-bonded parameter combinations, which will be copied to t_params.
 */

struct t_param {
    //! The atom list (eg. bonds: particle i = a[0] (ai), j = a[1] (aj))
    int        a[MAXATOMLIST];
    //! Force parameters (eg. b0 = c[0])
    real       c[MAXFORCEPARAM];
    //! A string (instead of parameters), read from the .rtp file in pdb2gmx
    char       s[MAXSLEN];
    const int &ai() const { return a[0]; }
    int       &ai() { return a[0]; }
    const int &aj() const { return a[1]; }
    int       &aj() { return a[1]; }
    const int &ak() const { return a[2]; }
    int       &ak() { return a[2]; }
    const int &al() const { return a[3]; }
    int       &al() { return a[3]; }
    const int &am() const { return a[4]; }
    int       &am() { return a[4]; }

    real      &c0() { return c[0]; }
    real      &c1() { return c[1]; }
    real      &c2() { return c[2]; }
};

struct t_params {        // NOLINT (clang-analyzer-optin.performance.Padding)
    //! The number of bonds in this record
    int          nr;
    //! The amount of elements in the array
    int          maxnr;
    //! Array of parameters (dim: nr or nr*nr)
    t_param     *param;

    /* CMAP tmp data, there are probably better places for this */
    //! Cmap grid spacing TODO remove here.
    int         grid_spacing;
    //! Number of cmap angles TODO remove here.
    int         nc;
    //! Temporary storage of the raw cmap grid data TODO remove here.
    real       *cmap;
    //! Number of allocated elements in cmap grid TODO remove here.
    int         ncmap;
    //! Store the five atomtypes followed by a number that identifies the type
    int        *cmap_types;
    //! Number of allocated elements in cmap_types
    int         nct;

};

typedef struct {
    int            nr;      /* The number of exclusions             */
    int           *e;       /* The excluded atoms                   */
} t_excls;

struct t_molinfo {
    //! Entry in symbol table for name.
    SymbolPtr         name;
    //! Number of exclusions per atom
    int               nrexcl = 0;
    //! Have exclusions been generated?
    bool              excl_set = false;
    //! Has the mol been processed
    bool              bProcessed = false;
    //! Molecule atoms
    t_atoms           atoms;
    //! Molecule charge groups
    t_block           cgs;
    //! Molecule molecules :)
    t_block           mols;
    //! Molecule exclusions
    t_blocka          excls;
    //! Molecule Parameters in old style
    t_params          plist[F_NRE];
};

typedef struct {
    char *name;
    int   nr;
} t_mols;

bool is_int(double x);
/* Returns TRUE when x is integer */

/* Must correspond to strings in topdirs.c */
typedef enum {
    d_defaults,
    d_atomtypes,
    d_bondtypes,
    d_constrainttypes,
    d_pairtypes,
    d_angletypes,
    d_dihedraltypes,
    d_nonbond_params,
    d_implicit_genborn_params,
    d_implicit_surface_params,
    d_cmaptypes,
    d_moleculetype,
    d_atoms,
    d_vsites2,
    d_vsites3,
    d_vsites4,
    d_vsitesn,
    d_bonds,
    d_exclusions,
    d_pairs,
    d_pairs_nb,
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
    d_cmap,
    d_intermolecular_interactions,
    d_maxdir,
    d_invalid,
    d_none
} directive;

#endif
