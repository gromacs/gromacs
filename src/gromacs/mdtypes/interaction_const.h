/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MDTYPES_INTERACTION_CONST_H
#define GMX_MDTYPES_INTERACTION_CONST_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"

struct t_lambda;
struct t_inputrec;
struct gmx_mtop_t;

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
    /* This struct contains the soft-core parameters from t_lambda,
     * but processed for direct use in the kernels.
     */
    struct SoftCoreParameters
    {
        // Constructor
        SoftCoreParameters(const t_lambda& fepvals);

        // Alpha parameter for Van der Waals interactions
        real alphaVdw;
        // Alpha parameter for Coulomb interactions
        real alphaCoulomb;
        // Exponent for the dependence of the soft-core on lambda
        int lambdaPower;
        // Value for sigma^6 for LJ interaction with C6<=0 and/or C12<=0
        real sigma6WithInvalidSigma;
        // Minimum value for sigma^6, used when soft-core is applied to Coulomb interactions
        real sigma6Minimum;
        // soft-core function
        SoftcoreType softcoreType;
        // (gapsys sc) linearization point scaling for vdW interactions
        real gapsysScaleLinpointVdW;
        // (gapsys sc) linearization point scaling for Coulomb interactions
        real gapsysScaleLinpointCoul;
        // (gapsys sc) lower bound/replacement for c12/c6 in vdw interactions
        real gapsysSigma6VdW;
    };

    /* VdW */
    VanDerWaalsType        vdwtype          = VanDerWaalsType::Cut;
    InteractionModifiers   vdw_modifier     = InteractionModifiers::None;
    double                 reppow           = 12;
    real                   rvdw             = 1;
    real                   rvdw_switch      = 0;
    struct shift_consts_t  dispersion_shift = { 0, 0, 0 };
    struct shift_consts_t  repulsion_shift  = { 0, 0, 0 };
    struct switch_consts_t vdw_switch       = { 0, 0, 0 };
    bool                   useBuckingham    = false;
    real                   buckinghamBMax   = 0;

    /* type of electrostatics */
    CoulombInteractionType eeltype          = CoulombInteractionType::Cut;
    InteractionModifiers   coulomb_modifier = InteractionModifiers::None;

    /* Coulomb */
    real rcoulomb        = 1;
    real rcoulomb_switch = 0;

    /* PME/Ewald */
    real         ewaldcoeff_q  = 0;
    real         ewaldcoeff_lj = 0;
    LongRangeVdW ljpme_comb_rule = LongRangeVdW::Geom; /* LJ combination rule for the LJ PME mesh part */
    real         sh_ewald        = 0; /* -sh_ewald is added to the direct space potential */
    real         sh_lj_ewald = 0;     /* sh_lj_ewald is added to the correction potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r = 1;
    real epsfac    = 1;

    /* Constants for reaction-field or plain cut-off */
    //! Dielectric constant for reaction field beyond the cutoff distance
    real reactionFieldPermitivity = 1;
    //! Coefficient for reaction field; scales relation between epsilon_r and reactionFieldPermitivity
    real reactionFieldCoefficient = 0;
    //! Constant shift to reaction field Coulomb interaction to make potential an integral of force
    real reactionFieldShift = 0;

    // Coulomb Ewald correction table
    std::unique_ptr<EwaldCorrectionTables> coulombEwaldTables;
    // Van der Waals Ewald correction table
    std::unique_ptr<EwaldCorrectionTables> vdwEwaldTables;

    // Free-energy parameters, only present when free-energy calculations are requested
    std::unique_ptr<SoftCoreParameters> softCoreParameters;
};

/*! \brief Construct interaction constants
 *
 * This data is used (particularly) by search and force code for
 * short-range interactions. Many of these are constant for the whole
 * simulation; some are constant only after PME tuning completes.
 */
interaction_const_t init_interaction_const(FILE*             fp,
                                           const t_inputrec& ir,
                                           const gmx_mtop_t& mtop,
                                           bool              systemHasNetCharge);

#endif
