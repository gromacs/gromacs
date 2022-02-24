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
#ifndef GMX_GMXLIB_NONBONDED_SOFTCORE_H
#define GMX_GMXLIB_NONBONDED_SOFTCORE_H

#include "config.h"

#include "gromacs/math/functions.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"

/* linearized electrostatics */
template<bool computeForces, class RealType, class BoolType>
static inline void quadraticApproximationCoulomb(const RealType qq,
                                                 const RealType rInvQ,
                                                 const RealType r,
                                                 const real     lambdaFac,
                                                 const real     dLambdaFac,
                                                 RealType gmx_unused* force,
                                                 RealType*            potential,
                                                 RealType*            dvdl,
                                                 BoolType             dvdlMask)
{
    RealType constFac = qq * rInvQ;
    RealType linFac   = constFac * r * rInvQ;
    RealType quadrFac = linFac * r * rInvQ;

    /* Computing Coulomb force and potential energy */
    if constexpr (computeForces)
    {
        *force = -2 * quadrFac + 3 * linFac;
    }

    *potential = quadrFac - 3 * (linFac - constFac);

    RealType lambdaFacRevInv = gmx::maskzInv(1 - lambdaFac, dvdlMask);
    *dvdl = dLambdaFac * 0.5_real * (lambdaFac * lambdaFacRevInv) * (quadrFac - 2 * linFac + constFac);
}

/* reaction-field linearized electrostatics */
template<bool computeForces, class RealType, class BoolType>
static inline void reactionFieldQuadraticPotential(const RealType qq,
                                                   const real     facel,
                                                   const RealType r,
                                                   const real     rCutoff,
                                                   const real     lambdaFac,
                                                   const real     dLambdaFac,
                                                   const RealType alphaEff,
                                                   const real     krf,
                                                   const real     potentialShift,
                                                   RealType gmx_unused* force,
                                                   RealType*            potential,
                                                   RealType*            dvdl,
                                                   BoolType             mask)
{
    RealType one(1);
    RealType zero(0);

    /* check if we have to use the hardcore values */
    BoolType computeValues = mask && (lambdaFac < one && zero < alphaEff && facel != zero);
    if (gmx::anyTrue(computeValues))
    {
        RealType lambdaFacRev = gmx::selectByMask(1 - lambdaFac, computeValues);

        RealType rQ = gmx::cbrt(lambdaFacRev);
        rQ          = gmx::sqrt(rQ) * (1 + gmx::abs(qq / facel));
        rQ          = rQ * alphaEff;

        // ensure that the linearization point doesn't go beyond rCutoff
        BoolType beyondCutoff = rCutoff < rQ;
        BoolType withinCutoff = rQ <= rCutoff;
        if (gmx::anyTrue(beyondCutoff))
        {
            rQ = gmx::blend(rQ, rCutoff, beyondCutoff);
        }

        computeValues = computeValues && (r < rQ);
        if (gmx::anyTrue(computeValues))
        {
            RealType rInvQ = gmx::maskzInv(rQ, computeValues);

            // Compute quadratic force, potential and dvdl
            RealType forceQuad(0);
            RealType potentialQuad(0);
            RealType dvdlQuad(0);
            quadraticApproximationCoulomb<computeForces>(
                    qq, rInvQ, r, lambdaFac, dLambdaFac, &forceQuad, &potentialQuad, &dvdlQuad, computeValues);

            // rf modification
            forceQuad     = forceQuad - qq * 2 * krf * r * r;
            potentialQuad = potentialQuad + qq * (krf * r * r - potentialShift);

            // update
            if constexpr (computeForces)
            {
                *force = gmx::blend(*force, forceQuad, computeValues);
            }
            *potential = gmx::blend(*potential, potentialQuad, computeValues);
            *dvdl      = *dvdl + gmx::selectByMask(dvdlQuad, computeValues && withinCutoff);
        }
    }
}

/* ewald linearized electrostatics */
template<bool computeForces, class RealType, class BoolType>
static inline void ewaldQuadraticPotential(const RealType qq,
                                           const real     facel,
                                           const RealType r,
                                           const real     rCutoff,
                                           const real     lambdaFac,
                                           const real     dLambdaFac,
                                           const RealType alphaEff,
                                           const real     potentialShift,
                                           RealType gmx_unused* force,
                                           RealType*            potential,
                                           RealType*            dvdl,
                                           BoolType             mask)
{
    RealType one(1);
    RealType zero(0);

    /* check if we have to use the hardcore values */
    BoolType computeValues = mask && (lambdaFac < one && zero < alphaEff && facel != zero);
    if (gmx::anyTrue(computeValues))
    {
        RealType lambdaFacRev = gmx::selectByMask(1 - lambdaFac, computeValues);

        RealType rQ = gmx::cbrt(lambdaFacRev);
        rQ          = gmx::sqrt(rQ) * (1 + gmx::abs(qq / facel));
        rQ          = rQ * alphaEff;

        // ensure that the linearization point doesn't go beyond rCutoff
        BoolType beyondCutoff = rCutoff < rQ;
        BoolType withinCutoff = rQ <= rCutoff;
        if (gmx::anyTrue(beyondCutoff))
        {
            rQ = gmx::blend(rQ, rCutoff, beyondCutoff);
        }

        computeValues = computeValues && (r < rQ);
        if (gmx::anyTrue(computeValues))
        {
            RealType rInvQ = gmx::maskzInv(rQ, computeValues);

            // Compute quadratic force, potential and dvdl
            RealType forceQuad(0);
            RealType potentialQuad(0);
            RealType dvdlQuad(0);
            quadraticApproximationCoulomb<computeForces>(
                    qq, rInvQ, r, lambdaFac, dLambdaFac, &forceQuad, &potentialQuad, &dvdlQuad, computeValues);

            // ewald modification
            potentialQuad = potentialQuad - qq * potentialShift;

            // update
            if constexpr (computeForces)
            {
                *force = gmx::blend(*force, forceQuad, computeValues);
            }
            *potential = gmx::blend(*potential, potentialQuad, computeValues);
            *dvdl      = *dvdl + gmx::selectByMask(dvdlQuad, computeValues && withinCutoff);
        }
    }
}

/* cutoff LJ with quadratic appximation of lj-potential */
template<bool computeForces, class RealType, class BoolType>
static inline void lennardJonesQuadraticPotential(const RealType c6,
                                                  const RealType c12,
                                                  const RealType r,
                                                  const RealType rsq,
                                                  const real     lambdaFac,
                                                  const real     dLambdaFac,
                                                  const RealType sigma6,
                                                  const RealType alphaEff,
                                                  const real     repulsionShift,
                                                  const real     dispersionShift,
                                                  RealType gmx_unused* force,
                                                  RealType*            potential,
                                                  RealType*            dvdl,
                                                  BoolType             mask)
{
    constexpr real c_twentySixSeventh = 26.0_real / 7.0_real;
    constexpr real c_oneSixth         = 1.0_real / 6.0_real;
    constexpr real c_oneTwelth        = 1.0_real / 12.0_real;
    constexpr real c_half             = 1.0_real / 2.0_real;

    RealType one(1);
    RealType zero(0);

    /* check if we have to use the hardcore values */
    BoolType computeValues = mask && (lambdaFac < one && zero < alphaEff);
    if (gmx::anyTrue(computeValues))
    {
        RealType lambdaFacRev    = gmx::selectByMask(1 - lambdaFac, computeValues);
        RealType lambdaFacRevInv = gmx::maskzInv(1 - lambdaFac, computeValues);

        RealType rQ = gmx::cbrt(c_twentySixSeventh * sigma6 * lambdaFacRev);
        rQ          = gmx::sqrt(rQ);
        rQ          = rQ * alphaEff;

        computeValues = (computeValues && r < rQ);
        if (gmx::anyTrue(computeValues))
        {
            /* scaled values for c6 and c12 */
            RealType c6s, c12s;
            c6s  = c_oneSixth * c6;
            c12s = c_oneTwelth * c12;
            /* Temporary variables for inverted values */
            RealType rInvQ = gmx::maskzInv(rQ, computeValues);
            RealType rInv14C, rInv13C, rInv12C;
            RealType rInv8C, rInv7C, rInv6C;
            rInv6C  = rInvQ * rInvQ * rInvQ;
            rInv6C  = rInv6C * rInv6C;
            rInv7C  = rInv6C * rInvQ;
            rInv8C  = rInv7C * rInvQ;
            rInv14C = c12s * rInv7C * rInv7C * rsq;
            rInv13C = c12s * rInv7C * rInv6C * r;
            rInv12C = c12s * rInv6C * rInv6C;
            rInv8C  = rInv8C * c6s * rsq;
            rInv7C  = rInv7C * c6s * r;
            rInv6C  = rInv6C * c6s;

            /* Temporary variables for A and B */
            RealType quadrFac, linearFac, constFac;
            quadrFac  = 156 * rInv14C - 42 * rInv8C;
            linearFac = 168 * rInv13C - 48 * rInv7C;
            constFac  = 91 * rInv12C - 28 * rInv6C;

            /* Computing LJ force and potential energy */
            RealType gmx_unused forceQuad     = -quadrFac + linearFac;
            RealType            potentialQuad = c_half * quadrFac - linearFac + constFac;
            RealType            dvdlQuad      = dLambdaFac * 28 * (lambdaFac * lambdaFacRevInv)
                                * ((6.5_real * rInv14C - rInv8C) - (13 * rInv13C - 2 * rInv7C)
                                   + (6.5_real * rInv12C - rInv6C));

            potentialQuad = potentialQuad
                            + gmx::selectByMask(((c12s * repulsionShift) - (c6s * dispersionShift)),
                                                computeValues);
            if constexpr (computeForces)
            {
                *force = gmx::blend(*force, forceQuad, computeValues);
            }
            *potential = gmx::blend(*potential, potentialQuad, computeValues);
            *dvdl      = *dvdl + gmx::selectByMask(dvdlQuad, computeValues);
        }
    }
}

#endif
