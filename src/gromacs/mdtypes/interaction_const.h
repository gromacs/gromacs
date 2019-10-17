/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include <memory>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

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

/* Convenience type for vector with aligned memory */
template<typename T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;

/* Force/energy interpolation tables for Ewald long-range corrections
 *
 * Interpolation is linear for the force, quadratic for the potential.
 */
struct EwaldCorrectionTables
{
    // 1/table_spacing, units 1/nm
    real scale = 0;
    // Force table
    AlignedVector<real> tableF;
    // Energy table
    AlignedVector<real> tableV;
    // Coulomb force+energy table, size of array is tabq_size*4,
    // entry quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
    // this is used with 4-wide SIMD for aligned loads
    AlignedVector<real> tableFDV0;
};

/* The physical interaction parameters for non-bonded interaction calculations
 *
 * This struct contains copies of the physical interaction parameters
 * from the user input as well as processed values that are need in
 * non-bonded interaction kernels.
 *
 * The default constructor gives plain Coulomb and LJ interactions cut off
 * a 1 nm without potential shifting and a Coulomb pre-factor of 1.
 */
struct interaction_const_t
{
    int cutoff_scheme = ecutsVERLET;

    /* VdW */
    int                    vdwtype          = evdwCUT;
    int                    vdw_modifier     = eintmodNONE;
    double                 reppow           = 12;
    real                   rvdw             = 1;
    real                   rvdw_switch      = 0;
    struct shift_consts_t  dispersion_shift = { 0, 0, 0 };
    struct shift_consts_t  repulsion_shift  = { 0, 0, 0 };
    struct switch_consts_t vdw_switch       = { 0, 0, 0 };
    gmx_bool               useBuckingham    = false;
    real                   buckinghamBMax   = 0;

    /* type of electrostatics */
    int eeltype          = eelCUT;
    int coulomb_modifier = eintmodNONE;

    /* Coulomb */
    real rcoulomb        = 1;
    real rcoulomb_switch = 0;

    /* PME/Ewald */
    real ewaldcoeff_q    = 0;
    real ewaldcoeff_lj   = 0;
    int  ljpme_comb_rule = eljpmeGEOM; /* LJ combination rule for the LJ PME mesh part */
    real sh_ewald        = 0;          /* -sh_ewald is added to the direct space potential */
    real sh_lj_ewald     = 0;          /* sh_lj_ewald is added to the correction potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r = 1;
    real epsfac    = 1;

    /* Constants for reaction-field or plain cut-off */
    real epsilon_rf = 1;
    real k_rf       = 0;
    real c_rf       = 0;

    // Coulomb Ewald correction table
    std::unique_ptr<EwaldCorrectionTables> coulombEwaldTables;
    /* Note that a Van der Waals Ewald correction table
     * of type EwaldCorrectionTables can be added here if wanted.
     */
};

#endif
