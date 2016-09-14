/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_INTERACTION_CONST_H
#define GMX_MDTYPES_INTERACTION_CONST_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Used with force switching or a constant potential shift:
 * rsw       = max(r - r_switch, 0)
 * force/p   = r^-(p+1) + c2*rsw^2 + c3*rsw^3
 * potential = r^-p + c2/3*rsw^3 + c3/4*rsw^4 + cpot
 * With a constant potential shift c2 and c3 are both 0.
 */
struct shift_consts_t
{
    real c2;
    real c3;
    real cpot;
};

/* Used with potential switching:
 * rsw        = max(r - r_switch, 0)
 * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
 * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
 * force      = force*dsw - potential*sw
 * potential *= sw
 */
struct switch_consts_t
{
    real c3;
    real c4;
    real c5;
};

struct interaction_const_t
{
    int             cutoff_scheme;

    /* VdW */
    int                    vdwtype;
    int                    vdw_modifier;
    double                 reppow;
    real                   rvdw;
    real                   rvdw_switch;
    struct shift_consts_t  dispersion_shift;
    struct shift_consts_t  repulsion_shift;
    struct switch_consts_t vdw_switch;
    gmx_bool               useBuckingham;
    real                   buckinghamBMax;
    /* TODO: remove this variable, used for not modyfing the group kernels,
     * it is equal to -dispersion_shift->cpot
     */
    real sh_invrc6;

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;
    int  coulomb_modifier;

    /* Coulomb */
    real rcoulomb;
    real rcoulomb_switch;

    /* PME/Ewald */
    real ewaldcoeff_q;
    real ewaldcoeff_lj;
    int  ljpme_comb_rule; /* LJ combination rule for the LJ PME mesh part */
    real sh_ewald;        /* -sh_ewald is added to the direct space potential */
    real sh_lj_ewald;     /* sh_lj_ewald is added to the correction potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r;
    real epsfac;

    /* Constants for reaction-field or plain cut-off */
    real epsilon_rf;
    real k_rf;
    real c_rf;

    /* Force/energy interpolation tables, linear in force, quadratic in V */
    real  tabq_scale;
    int   tabq_size;
    /* Coulomb force table, size of array is tabq_size (when used) */
    real *tabq_coul_F;
    /* Coulomb energy table, size of array is tabq_size (when used) */
    real *tabq_coul_V;
    /* Coulomb force+energy table, size of array is tabq_size*4,
       entry quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
       this is used with single precision x86 SIMD for aligned loads */
    real *tabq_coul_FDV0;

    /* Vdw force table for LJ-PME, size of array is tabq_size (when used) */
    real *tabq_vdw_F;
    /* Vdw energy table for LJ-PME, size of array is tabq_size (when used) */
    real *tabq_vdw_V;
    /* Vdw force+energy table for LJ-PME, size of array is tabq_size*4, entry
       quadruplets are: F[i], F[i+1]-F[i], V[i], 0, this is used with
       single precision x86 SIMD for aligned loads */
    real *tabq_vdw_FDV0;
};

#ifdef __cplusplus
}
#endif

#endif
