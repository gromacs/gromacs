/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef GMX_LEGACYHEADERS_TYPES_GENBORN_H
#define GMX_LEGACYHEADERS_TYPES_GENBORN_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int  nbonds;
    int  bond[10];
    real length[10];
} genborn_bonds_t;

typedef struct gbtmpnbls *gbtmpnbls_t;

/* Struct to hold all the information for GB */
typedef struct
{
    int nr;                   /* number of atoms, length of arrays below */
    int n12;                  /* number of 1-2 (bond) interactions       */
    int n13;                  /* number of 1-3 (angle) terms             */
    int n14;                  /* number of 1-4 (torsion) terms           */
    int nalloc;               /* Allocation of local arrays (with DD)    */


    /* Arrays below that end with _globalindex are used for setting up initial values of
     * all gb parameters and values. They all have length natoms, which for DD is the
     * global atom number.
     * Values are then taken from these arrays to local copies, that have names without
     * _globalindex, in the routine make_local_gb(), which is called once for single
     * node runs, and for DD at every call to dd_partition_system
     */

    real       *gpol;              /* Atomic polarisation energies */
    real       *gpol_globalindex;  /*  */
    real       *gpol_still_work;   /* Work array for Still model */
    real       *gpol_hct_work;     /* Work array for HCT/OBC models */
    real       *bRad;              /* Atomic Born radii */
    real       *vsolv;             /* Atomic solvation volumes */
    real       *vsolv_globalindex; /*  */
    real       *gb_radius;         /* Radius info, copied from atomtypes */
    real       *gb_radius_globalindex;

    int        *use;                /* Array that till if this atom does GB */
    int        *use_globalindex;    /* Global array for parallelization */

    real        es;                 /* Solvation energy and derivatives */
    real       *asurf;              /* Atomic surface area */
    rvec       *dasurf;             /* Surface area derivatives */
    real        as;                 /* Total surface area */

    real       *drobc;              /* Parameters for OBC chain rule calculation */
    real       *param;              /* Precomputed factor rai*atype->S_hct for HCT/OBC */
    real       *param_globalindex;  /*  */

    real       *log_table;          /* Table for logarithm lookup */

    real        obc_alpha;          /* OBC parameters */
    real        obc_beta;           /* OBC parameters */
    real        obc_gamma;          /* OBC parameters */
    real        gb_doffset;         /* Dielectric offset for Still/HCT/OBC */
    real        gb_epsilon_solvent; /*   */
    real        epsilon_r;          /* Used for inner dielectric */

    real        sa_surface_tension; /* Surface tension for non-polar solvation */

    real       *work;               /* Used for parallel summation and in the chain rule, length natoms         */
    real       *buf;                /* Used for parallel summation and in the chain rule, length natoms         */
    int        *count;              /* Used for setting up the special gb nblist, length natoms                 */
    gbtmpnbls_t nblist_work;        /* Used for setting up the special gb nblist, dim natoms*nblist_work_nalloc */
    int         nblist_work_nalloc; /* Length of second dimension of nblist_work                                */
}
gmx_genborn_t;

#ifdef __cplusplus
}
#endif
#endif
