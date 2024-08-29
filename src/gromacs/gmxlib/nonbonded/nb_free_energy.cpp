/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "nb_free_energy.h"

#include "config.h"

#include <cmath>
#include <cstdint>

#include <algorithm>
#include <memory>
#include <set>
#include <vector>

#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "nb_softcore.h"

//! Scalar (non-SIMD) data types.
struct ScalarDataTypes
{
    using RealType = real; //!< The data type to use as real.
    using IntType  = int;  //!< The data type to use as int.
    using BoolType = bool; //!< The data type to use as bool for real value comparison.
    static constexpr int simdRealWidth = 1; //!< The width of the RealType.
};

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS
//! SIMD data types.
struct SimdDataTypes
{
    using RealType = gmx::SimdReal;  //!< The data type to use as real.
    using IntType  = gmx::SimdInt32; //!< The data type to use as int.
    using BoolType = gmx::SimdBool;  //!< The data type to use as bool for real value comparison.
    static constexpr int simdRealWidth = GMX_SIMD_REAL_WIDTH; //!< The width of the RealType.
#    if GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE
    static constexpr int simdIntWidth = GMX_SIMD_DINT32_WIDTH; //!< The width of the IntType.
#    else
    static constexpr int simdIntWidth = GMX_SIMD_FINT32_WIDTH; //!< The width of the IntType.
#    endif
};
#endif

/*! \brief Lower limit for square interaction distances in nonbonded kernels.
 *
 * This is a mimimum on r^2 to avoid overflows when computing r^6.
 * This will only affect results for soft-cored interaction at distances smaller
 * than 1e-6 and will limit extremely high foreign energies for overlapping atoms.
 * Note that we could use a somewhat smaller minimum in double precision.
 * But because invsqrt in double precision can use single precision, this number
 * can not be much smaller, we use the same number for simplicity.
 */
constexpr real c_minDistanceSquared = 1.0e-12_real;

/*! \brief Higher limit for r^-6 used for Lennard-Jones interactions
 *
 * This is needed to avoid overflow of LJ energy and force terms for excluded
 * atoms and foreign energies of hard-core states of overlapping atoms.
 * Note that in single precision this value leaves room for C12 coefficients up to 3.4e8.
 */
constexpr real c_maxRInvSix = 1.0e15_real;

//! The free-energy kernel Lennard-Jones functional shape
enum class LJKernelType
{
    Cutoff,          //!< Plain cut-off with or without potential shift
    ForceSwitch,     //!< Force switch
    PotentialSwitch, //!< Potential Switch
    Ewald            //!< Ewald for dispersion
};

template<bool computeScalarForce, class RealType>
static inline void
pmeCoulombCorrectionVF(const RealType rSq, const real beta, RealType* pot, RealType gmx_unused* force)
{
    const RealType brsq = rSq * beta * beta;
    if constexpr (computeScalarForce)
    {
        *force = -brsq * beta * gmx::pmeForceCorrection(brsq);
    }
    *pot = beta * gmx::pmePotentialCorrection(brsq);
}

template<bool computeScalarForce, class RealType, class BoolType>
static inline void pmeLJCorrectionVF(const RealType rInv,
                                     const RealType rSq,
                                     const real     ewaldLJCoeffSq,
                                     const real     ewaldLJCoeffSixDivSix,
                                     RealType*      pot,
                                     RealType gmx_unused* force,
                                     const BoolType       mask,
                                     const BoolType       bIiEqJnr)
{
    // The expression 1 - exp(-x) * (1 + x + x^2/2) gives divergent
    // rounding errors when x goes towards zero. To avoid this, we approximate
    // this expression with its series expansion around x = 0:
    //   x^3 (1/6 - x/8 + x^2/20 - x^3/72 + ...)
    // The relative error of the full expression is x^-3 as x -> 0.
    // The relative error of the quadratic approximation is 1/8 x^3.
    // We set the switch point for coeffSqRSq such that the maximum error
    // of rounding in the full term and the approximation error are equal.
    // At this point both errors are sqrt(2)/4 sqrt(GMX_REAL_EPS).
    // This error is large relative to this force component, but is small
    // relative to the total forces on atoms.
    const real c_coeffSqRSqSwitch = std::pow(8.0_real * GMX_REAL_EPS, 1.0_real / 6.0_real);

    // We mask rInv to get zero force and potential for masked out pair interactions
    const RealType rInvSq  = rInv * rInv;
    const RealType rInvSix = rInvSq * rInvSq * rInvSq;
    // Mask rSq to avoid underflow in exp()
    const RealType coeffSqRSq       = ewaldLJCoeffSq * gmx::selectByMask(rSq, mask);
    const RealType expNegCoeffSqRSq = gmx::exp(-coeffSqRSq);
    const RealType poly             = 1.0_real + coeffSqRSq + 0.5_real * coeffSqRSq * coeffSqRSq;
    const RealType fullTerm         = rInvSix * (1.0_real - expNegCoeffSqRSq * poly);
    const RealType approximation =
            ewaldLJCoeffSixDivSix * (1.0_real + coeffSqRSq * (-0.75_real + 0.3_real * coeffSqRSq));
    const RealType term = gmx::blend(fullTerm, approximation, coeffSqRSq < c_coeffSqRSqSwitch);

    if constexpr (computeScalarForce)
    {
        *force = term - expNegCoeffSqRSq * ewaldLJCoeffSixDivSix;
        *force = *force * rInvSq;
    }
    // The self interaction is the limit for r -> 0 which we need to compute separately
    *pot = gmx::blend(term, 0.5_real * ewaldLJCoeffSixDivSix, bIiEqJnr);
}

//! Computes r^(1/6) and 1/r^(1/6)
template<class RealType>
static inline void sixthRoot(const RealType r, RealType* sixthRoot, RealType* invSixthRoot)
{
    RealType cbrtRes = gmx::cbrt(r);
    *invSixthRoot    = gmx::invsqrt(cbrtRes);
    *sixthRoot       = gmx::inv(*invSixthRoot);
}

template<class RealType>
static inline RealType calculateRinv6(const RealType rInvV)
{
    RealType rInv6 = rInvV * rInvV;
    return (rInv6 * rInv6 * rInv6);
}

template<class RealType>
static inline RealType calculateVdw6(const RealType c6, const RealType rInv6)
{
    return (c6 * rInv6);
}

template<class RealType>
static inline RealType calculateVdw12(const RealType c12, const RealType rInv6)
{
    return (c12 * rInv6 * rInv6);
}

/* reaction-field electrostatics */
template<class RealType>
static inline RealType reactionFieldScalarForce(const RealType qq,
                                                const RealType rInv,
                                                const RealType r,
                                                const real     reactionFieldCoefficient,
                                                const real     two)
{
    return (qq * (rInv - two * reactionFieldCoefficient * r * r));
}
template<class RealType>
static inline RealType reactionFieldPotential(const RealType qq,
                                              const RealType rInv,
                                              const RealType r,
                                              const real     reactionFieldCoefficient,
                                              const real     potentialShift)
{
    return (qq * (rInv + reactionFieldCoefficient * r * r - potentialShift));
}

/* Ewald electrostatics */
template<class RealType>
static inline RealType ewaldScalarForce(const RealType coulomb, const RealType rInv)
{
    return (coulomb * rInv);
}
template<class RealType>
static inline RealType ewaldPotential(const RealType coulomb, const RealType rInv, const real potentialShift)
{
    return (coulomb * (rInv - potentialShift));
}

/* cutoff LJ */
template<class RealType>
static inline RealType lennardJonesScalarForce(const RealType v6, const RealType v12)
{
    return (v12 - v6);
}
template<class RealType>
static inline RealType lennardJonesPotential(const RealType v6,
                                             const RealType v12,
                                             const RealType c6,
                                             const RealType c12,
                                             const real     repulsionShift,
                                             const real     dispersionShift,
                                             const real     oneSixth,
                                             const real     oneTwelfth)
{
    return ((v12 + c12 * repulsionShift) * oneTwelfth - (v6 + c6 * dispersionShift) * oneSixth);
}

/* Ewald LJ */
template<class RealType>
static inline RealType ewaldLennardJonesGridSubtract(const RealType c6grid,
                                                     const real     potentialShift,
                                                     const real     oneSixth)
{
    return (c6grid * potentialShift * oneSixth);
}

//! Generates intermediate quantities for force switch interactions
template<class RealType>
static inline void computeForceSwitchVariables(const RealType r,
                                               const real     rSwitch,
                                               RealType*      rSwitched,
                                               RealType*      rSwitchedSquared,
                                               RealType*      rSwitchedSquaredTimesR)
{
    *rSwitched        = gmx::max(r - rSwitch, RealType(0.0_real));
    *rSwitchedSquared = *rSwitched * *rSwitched;

    *rSwitchedSquaredTimesR = *rSwitchedSquared * r;
}

//! Return the force for force-switch interactions
template<class RealType, class BoolType>
static inline RealType forceSwitchScalarForceMod(const RealType fScalarInp,
                                                 const RealType rSwitched,
                                                 const RealType rSwitchedSquaredTimesR,
                                                 const real     c2,
                                                 const real     c3,
                                                 const BoolType mask)
{
    /* The mask should select on rV < rVdw */
    return (gmx::selectByMask(fScalarInp + (c3 * rSwitched + c2) * rSwitchedSquaredTimesR, mask));
}

//! Return the modification of the potential for force-switch interactions
template<class RealType, class BoolType>
static inline RealType forceSwitchPotentialMod(const RealType rSwitched,
                                               const RealType rSwitchedSquaredTimesR,
                                               const real     c3,
                                               const real     c4,
                                               const BoolType mask)
{
    /* The mask should select on rV < rVdw */
    return (gmx::selectByMask((c4 * rSwitched + c3) * (rSwitchedSquaredTimesR * rSwitched), mask));
}

/* LJ Potential switch */
template<class RealType, class BoolType>
static inline RealType potSwitchScalarForceMod(const RealType fScalarInp,
                                               const RealType potential,
                                               const RealType sw,
                                               const RealType r,
                                               const RealType dsw,
                                               const BoolType mask)
{
    /* The mask should select on rV < rVdw */
    return (gmx::selectByMask(fScalarInp * sw - r * potential * dsw, mask));
}
template<class RealType, class BoolType>
static inline RealType potSwitchPotentialMod(const RealType potentialInp, const RealType sw, const BoolType mask)
{
    /* The mask should select on rV < rVdw */
    return (gmx::selectByMask(potentialInp * sw, mask));
}

//! Templated free-energy non-bonded kernel
template<typename DataTypes, KernelSoftcoreType softcoreType, bool scLambdasOrAlphasDiffer, bool elecInteractionTypeIsEwald, LJKernelType ljKernelType, bool computeForces>
static void nb_free_energy_kernel(const t_nblist&                                  nlist,
                                  const gmx::ArrayRefWithPadding<const gmx::RVec>& coords,
                                  const int                                        ntype,
                                  const interaction_const_t&           interactionParameters,
                                  gmx::ArrayRef<const gmx::RVec>       shiftvec,
                                  gmx::ArrayRef<const real>            nbfp,
                                  gmx::ArrayRef<const real> gmx_unused nbfp_grid,
                                  gmx::ArrayRef<const real>            chargeA,
                                  gmx::ArrayRef<const real>            chargeB,
                                  gmx::ArrayRef<const int>             typeA,
                                  gmx::ArrayRef<const int>             typeB,
                                  int                                  flags,
                                  gmx::ArrayRef<const real>            lambda,
                                  t_nrnb* gmx_restrict                 nrnb,
                                  gmx::ArrayRefWithPadding<gmx::RVec>  threadForceBuffer,
                                  rvec gmx_unused*    threadForceShiftBuffer,
                                  gmx::ArrayRef<real> threadVCoul,
                                  gmx::ArrayRef<real> threadVVdw,
                                  gmx::ArrayRef<real> threadDvdl)
{
#define STATE_A 0
#define STATE_B 1
#define NSTATES 2

    using RealType = typename DataTypes::RealType;
    using IntType  = typename DataTypes::IntType;
    using BoolType = typename DataTypes::BoolType;

    // We need the scalar force to compute the Beutler soft-core contribution to dV/dlambda
    constexpr bool computeScalarForce = computeForces || softcoreType == KernelSoftcoreType::Beutler;

    constexpr real oneTwelfth = 1.0_real / 12.0_real;
    constexpr real oneSixth   = 1.0_real / 6.0_real;
    constexpr real zero       = 0.0_real;
    constexpr real half       = 0.5_real;
    constexpr real one        = 1.0_real;
    constexpr real two        = 2.0_real;
    constexpr real six        = 6.0_real;

    // Extract pair list data
    const int                nri    = nlist.nri;
    gmx::ArrayRef<const int> iinr   = nlist.iinr;
    gmx::ArrayRef<const int> jindex = nlist.jindex;
    gmx::ArrayRef<const int> jjnr   = nlist.jjnr;
    gmx::ArrayRef<const int> shift  = nlist.shift;
    gmx::ArrayRef<const int> gid    = nlist.gid;

    const real lambdaCoul = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
    const real lambdaVdw  = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];

    // Extract softcore parameters
    const auto&           scParams               = *interactionParameters.softCoreParameters;
    const real            lambdaPower            = scParams.lambdaPower;
    const real gmx_unused alphaCoulomb           = scParams.alphaCoulomb;
    const real gmx_unused alphaVdw               = scParams.alphaVdw;
    const real gmx_unused sigma6WithInvalidSigma = scParams.sigma6WithInvalidSigma;
    const real gmx_unused sigma6Minimum          = scParams.sigma6Minimum;

    const real gmx_unused gapsysScaleLinpointCoul = scParams.gapsysScaleLinpointCoul;
    const real gmx_unused gapsysScaleLinpointVdW  = scParams.gapsysScaleLinpointVdW;
    const real gmx_unused gapsysSigma6VdW         = scParams.gapsysSigma6VdW;

    const bool gmx_unused doShiftForces = ((flags & GMX_NONBONDED_DO_SHIFTFORCE) != 0);
    const bool            doPotential   = ((flags & GMX_NONBONDED_DO_POTENTIAL) != 0);

    // Extract data from interaction_const_t
    const real            elecEpsilonFactor        = interactionParameters.epsfac;
    const real            rCoulomb                 = interactionParameters.rcoulomb;
    const real            reactionFieldCoefficient = interactionParameters.reactionFieldCoefficient;
    const real gmx_unused reactionFieldShift       = interactionParameters.reactionFieldShift;
    const real gmx_unused shLjEwald                = interactionParameters.sh_lj_ewald;
    const real            rVdw                     = interactionParameters.rvdw;
    const real            dispersionShift          = interactionParameters.dispersion_shift.cpot;
    const real gmx_unused dispersionShift2         = interactionParameters.dispersion_shift.c2;
    const real gmx_unused dispersionShift3         = interactionParameters.dispersion_shift.c3;
    const real            repulsionShift           = interactionParameters.repulsion_shift.cpot;
    const real gmx_unused repulsionShift2          = interactionParameters.repulsion_shift.c2;
    const real gmx_unused repulsionShift3          = interactionParameters.repulsion_shift.c3;
    const real            ewaldBeta                = interactionParameters.ewaldcoeff_q;
    real gmx_unused       ewaldLJCoeffSq;
    real gmx_unused       ewaldLJCoeffSixDivSix;
    if constexpr (ljKernelType == LJKernelType::Ewald)
    {
        ewaldLJCoeffSq = interactionParameters.ewaldcoeff_lj * interactionParameters.ewaldcoeff_lj;
        ewaldLJCoeffSixDivSix = ewaldLJCoeffSq * ewaldLJCoeffSq * ewaldLJCoeffSq / six;
    }

    // Note that the nbnxm kernels do not support Coulomb potential switching at all
    GMX_ASSERT(interactionParameters.coulomb_modifier != InteractionModifiers::PotSwitch,
               "Potential switching is not supported for Coulomb with FEP");

    const real      rVdwSwitch = interactionParameters.rvdw_switch;
    real gmx_unused vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    if constexpr (ljKernelType == LJKernelType::PotentialSwitch)
    {
        const real d = rVdw - rVdwSwitch;
        vdw_swV3     = -10.0_real / (d * d * d);
        vdw_swV4     = 15.0_real / (d * d * d * d);
        vdw_swV5     = -6.0_real / (d * d * d * d * d);
        vdw_swF2     = -30.0_real / (d * d * d);
        vdw_swF3     = 60.0_real / (d * d * d * d);
        vdw_swF4     = -30.0_real / (d * d * d * d * d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        vdw_swV3 = vdw_swV4 = vdw_swV5 = vdw_swF2 = vdw_swF3 = vdw_swF4 = zero;
    }

    real gmx_unused minusDispersionShift2Div3, minusRepulsionShift2Div3;
    real gmx_unused minusDispersionShift3Div4, minusRepulsionShift3Div4;
    if constexpr (ljKernelType == LJKernelType::ForceSwitch)
    {
        constexpr real minusThird  = -1.0_real / 3.0_real;
        constexpr real minusFourth = -1.0_real / 4.0_real;
        minusDispersionShift2Div3  = dispersionShift2 * minusThird;
        minusDispersionShift3Div4  = dispersionShift3 * minusFourth;
        minusRepulsionShift2Div3   = repulsionShift2 * minusThird;
        minusRepulsionShift3Div4   = repulsionShift3 * minusFourth;
    }
    else
    {
        /* Avoid warnings from stupid compilers */
        minusDispersionShift2Div3 = minusRepulsionShift2Div3 = zero;
        minusDispersionShift3Div4 = minusRepulsionShift3Div4 = zero;
    }

    NbkernelElecType coulombInteractionType;
    if (interactionParameters.eeltype == CoulombInteractionType::Cut
        || usingRF(interactionParameters.eeltype))
    {
        coulombInteractionType = NbkernelElecType::ReactionField;
    }
    else
    {
        coulombInteractionType = NbkernelElecType::None;
    }

    real rCutoffMaxSq = std::max(interactionParameters.rcoulomb, interactionParameters.rvdw);
    rCutoffMaxSq      = rCutoffMaxSq * rCutoffMaxSq;
    const real gmx_unused rCutoffCoul = interactionParameters.rcoulomb;

    real gmx_unused sh_ewald = zero;
    if constexpr (elecInteractionTypeIsEwald || ljKernelType == LJKernelType::Ewald)
    {
        sh_ewald = interactionParameters.sh_ewald;
    }

    /* For Ewald/PME interactions we cannot easily apply the soft-core component to
     * reciprocal space. When we use non-switched Ewald interactions, we
     * assume the soft-coring does not significantly affect the grid contribution
     * and apply the soft-core only to the full 1/r (- shift) pair contribution.
     */

    const RealType            minDistanceSquared(c_minDistanceSquared);
    const RealType            maxRInvSix(c_maxRInvSix);
    const RealType gmx_unused floatMin(GMX_FLOAT_MIN);

    RealType dvdlCoul{ zero };
    RealType dvdlVdw{ zero };

    /* Lambda factor for state A, 1-lambda*/
    real lambdaFactorCoul[NSTATES], lambdaFactorVdw[NSTATES];
    lambdaFactorCoul[STATE_A] = one - lambdaCoul;
    lambdaFactorVdw[STATE_A]  = one - lambdaVdw;

    /* Lambda factor for state B, lambda*/
    lambdaFactorCoul[STATE_B] = lambdaCoul;
    lambdaFactorVdw[STATE_B]  = lambdaVdw;

    /*derivative of the lambda factor for state A and B */
    real dLambdaFactor[NSTATES];
    dLambdaFactor[STATE_A] = -one;
    dLambdaFactor[STATE_B] = one;

    real gmx_unused softcoreLambdaFactorCoul[NSTATES], softcoreDlFactorCoul[NSTATES],
            softcoreLambdaFactorVdw[NSTATES], softcoreDlFactorVdw[NSTATES];
    constexpr real softcoreRPower = six;
    for (int i = 0; i < NSTATES; i++)
    {
        softcoreLambdaFactorCoul[i] =
                (lambdaPower == 2 ? (1 - lambdaFactorCoul[i]) * (1 - lambdaFactorCoul[i])
                                  : (1 - lambdaFactorCoul[i]));
        softcoreDlFactorCoul[i] = dLambdaFactor[i] * lambdaPower / softcoreRPower
                                  * (lambdaPower == 2 ? (1 - lambdaFactorCoul[i]) : 1);
        softcoreLambdaFactorVdw[i] =
                (lambdaPower == 2 ? (1 - lambdaFactorVdw[i]) * (1 - lambdaFactorVdw[i])
                                  : (1 - lambdaFactorVdw[i]));
        softcoreDlFactorVdw[i] = dLambdaFactor[i] * lambdaPower / softcoreRPower
                                 * (lambdaPower == 2 ? (1 - lambdaFactorVdw[i]) : 1);
    }

    // We need pointers to real for SIMD access
    const real* gmx_restrict x = coords.paddedConstArrayRef().data()[0];
    real* gmx_restrict       forceRealPtr;
    if constexpr (computeForces)
    {
        GMX_ASSERT(nri == 0 || !threadForceBuffer.empty(), "need a valid threadForceBuffer");

        forceRealPtr = threadForceBuffer.paddedArrayRef().data()[0];

        if (doShiftForces)
        {
            GMX_ASSERT(threadForceShiftBuffer != nullptr, "need a valid threadForceShiftBuffer");
        }
    }

    // Used with reaction-field only
    BoolType haveExcludedPairsBeyondCutoff = false;

    for (int n = 0; n < nri; n++)
    {
        bool havePairsWithinCutoff = false;

        const int  is   = shift[n];
        const real shX  = shiftvec[is][XX];
        const real shY  = shiftvec[is][YY];
        const real shZ  = shiftvec[is][ZZ];
        const int  nj0  = jindex[n];
        const int  nj1  = jindex[n + 1];
        const int  ii   = iinr[n];
        const int  ii3  = 3 * ii;
        const real ix   = shX + x[ii3 + 0];
        const real iy   = shY + x[ii3 + 1];
        const real iz   = shZ + x[ii3 + 2];
        const real iqA  = elecEpsilonFactor * chargeA[ii];
        const real iqB  = elecEpsilonFactor * chargeB[ii];
        const int  ntiA = ntype * typeA[ii];
        const int  ntiB = ntype * typeB[ii];
        RealType   vCoulTot(zero);
        RealType   vVdwTot(zero);
        RealType   fIX(zero);
        RealType   fIY(zero);
        RealType   fIZ(zero);

#if GMX_SIMD_HAVE_REAL
        alignas(GMX_SIMD_ALIGNMENT) int            preloadIi[DataTypes::simdRealWidth];
        alignas(GMX_SIMD_ALIGNMENT) int gmx_unused preloadIs[DataTypes::simdRealWidth];
#else
        int            preloadIi[DataTypes::simdRealWidth];
        int gmx_unused preloadIs[DataTypes::simdRealWidth];
#endif
        for (int i = 0; i < DataTypes::simdRealWidth; i++)
        {
            preloadIi[i] = ii;
            preloadIs[i] = shift[n];
        }
        IntType ii_s = gmx::load<IntType>(preloadIi);

        for (int k = nj0; k < nj1; k += DataTypes::simdRealWidth)
        {
            RealType r, rInv;

#if GMX_SIMD_HAVE_REAL
            alignas(GMX_SIMD_ALIGNMENT) real    preloadPairIsValid[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real    preloadPairIncluded[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) int32_t preloadJnr[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) int32_t typeIndices[NSTATES][DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real    preloadQq[NSTATES][DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real gmx_unused preloadSigma6[NSTATES][DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real gmx_unused preloadAlphaVdwEff[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real gmx_unused preloadAlphaCoulEff[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT)
                    real gmx_unused preloadGapsysScaleLinpointVdW[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT)
                    real gmx_unused preloadGapsysScaleLinpointCoul[DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT)
                    real gmx_unused preloadGapsysSigma6VdW[NSTATES][DataTypes::simdRealWidth];
            alignas(GMX_SIMD_ALIGNMENT) real preloadLjPmeC6Grid[NSTATES][DataTypes::simdRealWidth];
#else
            real            preloadPairIsValid[DataTypes::simdRealWidth];
            real            preloadPairIncluded[DataTypes::simdRealWidth];
            int             preloadJnr[DataTypes::simdRealWidth];
            int             typeIndices[NSTATES][DataTypes::simdRealWidth];
            real            preloadQq[NSTATES][DataTypes::simdRealWidth];
            real gmx_unused preloadSigma6[NSTATES][DataTypes::simdRealWidth];
            real gmx_unused preloadAlphaVdwEff[DataTypes::simdRealWidth];
            real gmx_unused preloadAlphaCoulEff[DataTypes::simdRealWidth];
            real gmx_unused preloadGapsysScaleLinpointVdW[DataTypes::simdRealWidth];
            real gmx_unused preloadGapsysScaleLinpointCoul[DataTypes::simdRealWidth];
            real gmx_unused preloadGapsysSigma6VdW[NSTATES][DataTypes::simdRealWidth];
            real            preloadLjPmeC6Grid[NSTATES][DataTypes::simdRealWidth];
#endif
            for (int j = 0; j < DataTypes::simdRealWidth; j++)
            {
                if (k + j < nj1)
                {
                    preloadPairIsValid[j] = true;
                    /* Check if this pair on the exclusions list.*/
                    preloadPairIncluded[j]  = (nlist.excl_fep.empty() || nlist.excl_fep[k + j]);
                    const int jnr           = jjnr[k + j];
                    preloadJnr[j]           = jnr;
                    typeIndices[STATE_A][j] = ntiA + typeA[jnr];
                    typeIndices[STATE_B][j] = ntiB + typeB[jnr];
                    preloadQq[STATE_A][j]   = iqA * chargeA[jnr];
                    preloadQq[STATE_B][j]   = iqB * chargeB[jnr];

                    for (int i = 0; i < NSTATES; i++)
                    {
                        if constexpr (ljKernelType == LJKernelType::Ewald)
                        {
                            preloadLjPmeC6Grid[i][j] = nbfp_grid[2 * typeIndices[i][j]];
                        }
                        else
                        {
                            preloadLjPmeC6Grid[i][j] = 0;
                        }
                        if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                        {
                            const real c6  = nbfp[2 * typeIndices[i][j]];
                            const real c12 = nbfp[2 * typeIndices[i][j] + 1];
                            if (c6 > 0 && c12 > 0)
                            {
                                /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                                preloadSigma6[i][j] = 0.5_real * c12 / c6;
                                if (preloadSigma6[i][j]
                                    < sigma6Minimum) /* for disappearing coul and vdw with soft core at the same time */
                                {
                                    preloadSigma6[i][j] = sigma6Minimum;
                                }
                            }
                            else
                            {
                                preloadSigma6[i][j] = sigma6WithInvalidSigma;
                            }
                        }
                        if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                        {
                            const real c6  = nbfp[2 * typeIndices[i][j]];
                            const real c12 = nbfp[2 * typeIndices[i][j] + 1];
                            if (c6 > 0 && c12 > 0)
                            {
                                /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                                preloadGapsysSigma6VdW[i][j] = 0.5_real * c12 / c6;
                            }
                            else
                            {
                                preloadGapsysSigma6VdW[i][j] = gapsysSigma6VdW;
                            }
                        }
                    }
                    if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                    {
                        /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
                        const real c12A = nbfp[2 * typeIndices[STATE_A][j] + 1];
                        const real c12B = nbfp[2 * typeIndices[STATE_B][j] + 1];
                        if (c12A > 0 && c12B > 0)
                        {
                            preloadAlphaVdwEff[j]  = 0;
                            preloadAlphaCoulEff[j] = 0;
                        }
                        else
                        {
                            preloadAlphaVdwEff[j]  = alphaVdw;
                            preloadAlphaCoulEff[j] = alphaCoulomb;
                        }
                    }
                    if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                    {
                        /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
                        const real c12A = nbfp[2 * typeIndices[STATE_A][j] + 1];
                        const real c12B = nbfp[2 * typeIndices[STATE_B][j] + 1];
                        if (c12A > 0 && c12B > 0)
                        {
                            preloadGapsysScaleLinpointVdW[j]  = 0;
                            preloadGapsysScaleLinpointCoul[j] = 0;
                        }
                        else
                        {
                            preloadGapsysScaleLinpointVdW[j]  = gapsysScaleLinpointVdW;
                            preloadGapsysScaleLinpointCoul[j] = gapsysScaleLinpointCoul;
                        }
                    }
                }
                else
                {
                    preloadJnr[j]                     = jjnr[k];
                    preloadPairIsValid[j]             = false;
                    preloadPairIncluded[j]            = false;
                    preloadAlphaVdwEff[j]             = 0;
                    preloadAlphaCoulEff[j]            = 0;
                    preloadGapsysScaleLinpointVdW[j]  = 0;
                    preloadGapsysScaleLinpointCoul[j] = 0;

                    typeIndices[STATE_A][j] = ntiA + typeA[jjnr[k]];
                    typeIndices[STATE_B][j] = ntiB + typeB[jjnr[k]];
                    for (int i = 0; i < NSTATES; i++)
                    {
                        preloadLjPmeC6Grid[i][j]     = 0;
                        preloadQq[i][j]              = 0;
                        preloadSigma6[i][j]          = 0;
                        preloadGapsysSigma6VdW[i][j] = 0;
                    }
                }
            }

            RealType jx, jy, jz;
            gmx::gatherLoadUTranspose<3>(reinterpret_cast<const real*>(x), preloadJnr, &jx, &jy, &jz);

            const RealType pairIsValid   = gmx::load<RealType>(preloadPairIsValid);
            const RealType pairIncluded  = gmx::load<RealType>(preloadPairIncluded);
            const BoolType bPairIncluded = (pairIncluded != zero);
            const BoolType bPairExcluded = (pairIncluded == zero && pairIsValid != zero);

            const RealType dX  = ix - jx;
            const RealType dY  = iy - jy;
            const RealType dZ  = iz - jz;
            RealType       rSq = dX * dX + dY * dY + dZ * dZ;

            BoolType withinCutoffMask = (rSq < rCutoffMaxSq);

            if (!gmx::anyTrue(withinCutoffMask || bPairExcluded))
            {
                /* We save significant time by skipping all code below.
                 * Note that with soft-core interactions, the actual cut-off
                 * check might be different. But since the soft-core distance
                 * is always larger than r, checking on r here is safe.
                 * Exclusions outside the cutoff can not be skipped as
                 * when using Ewald: the reciprocal-space
                 * Ewald component still needs to be subtracted.
                 */
                continue;
            }
            else
            {
                havePairsWithinCutoff = true;
            }

            const IntType  jnr_s    = gmx::load<IntType>(preloadJnr);
            const BoolType bIiEqJnr = gmx::cvtIB2B(ii_s == jnr_s);

            RealType            c6[NSTATES];
            RealType            c12[NSTATES];
            RealType gmx_unused sigma6[NSTATES];
            RealType            qq[NSTATES];
            RealType gmx_unused ljPmeC6Grid[NSTATES];
            RealType gmx_unused alphaVdwEff;
            RealType gmx_unused alphaCoulEff;
            RealType gmx_unused gapsysScaleLinpointVdWEff;
            RealType gmx_unused gapsysScaleLinpointCoulEff;
            RealType gmx_unused gapsysSigma6VdWEff[NSTATES];
            for (int i = 0; i < NSTATES; i++)
            {
                gmx::gatherLoadTranspose<2>(nbfp.data(), typeIndices[i], &c6[i], &c12[i]);
                qq[i]          = gmx::load<RealType>(preloadQq[i]);
                ljPmeC6Grid[i] = gmx::load<RealType>(preloadLjPmeC6Grid[i]);
                if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                {
                    sigma6[i] = gmx::load<RealType>(preloadSigma6[i]);
                }
                if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                {
                    gapsysSigma6VdWEff[i] = gmx::load<RealType>(preloadGapsysSigma6VdW[i]);
                }
            }
            if constexpr (softcoreType == KernelSoftcoreType::Beutler)
            {
                alphaVdwEff  = gmx::load<RealType>(preloadAlphaVdwEff);
                alphaCoulEff = gmx::load<RealType>(preloadAlphaCoulEff);
            }
            if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
            {
                gapsysScaleLinpointVdWEff  = gmx::load<RealType>(preloadGapsysScaleLinpointVdW);
                gapsysScaleLinpointCoulEff = gmx::load<RealType>(preloadGapsysScaleLinpointCoul);
            }

            // Avoid overflow of r^-12 at distances near zero
            rSq  = gmx::max(rSq, minDistanceSquared);
            rInv = gmx::invsqrt(rSq);
            r    = rSq * rInv;

            RealType gmx_unused rp, rpm2;
            if constexpr (softcoreType == KernelSoftcoreType::Beutler)
            {
                rpm2 = rSq * rSq;  /* r4 */
                rp   = rpm2 * rSq; /* r6 */
            }
            else
            {
                /* The soft-core power p will not affect the results
                 * with not using soft-core, so we use power of 0 which gives
                 * the simplest math and cheapest code.
                 */
                rpm2 = rInv * rInv;
                rp   = one;
            }

            RealType scalarForcePerDistance(0);

            /* The following block is masked to only calculate values having bPairIncluded. If
             * bPairIncluded is true then withinCutoffMask must also be true. */
            if (gmx::anyTrue(withinCutoffMask && bPairIncluded))
            {
                RealType gmx_unused scalarForcePerDistanceCoul[NSTATES],
                        scalarForcePerDistanceVdw[NSTATES];
                RealType vCoul[NSTATES], vVdw[NSTATES];
                for (int i = 0; i < NSTATES; i++)
                {
                    scalarForcePerDistanceCoul[i] = zero;
                    scalarForcePerDistanceVdw[i]  = zero;
                    vCoul[i]                      = zero;
                    vVdw[i]                       = zero;

                    RealType gmx_unused rInvC, rInvV, rC, rV, rPInvC, rPInvV;

                    /* The following block is masked to require (qq[i] != 0 || c6[i] != 0 || c12[i]
                     * != 0) in addition to bPairIncluded, which in turn requires withinCutoffMask. */
                    BoolType nonZeroState = ((qq[i] != zero || c6[i] != zero || c12[i] != zero)
                                             && bPairIncluded && withinCutoffMask);
                    if (gmx::anyTrue(nonZeroState))
                    {
                        if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                        {
                            RealType divisor =
                                    (alphaCoulEff * softcoreLambdaFactorCoul[i] * sigma6[i] + rp);
                            rPInvC = gmx::inv(divisor);
                            sixthRoot(rPInvC, &rInvC, &rC);

                            if constexpr (scLambdasOrAlphasDiffer)
                            {
                                RealType divisor =
                                        (alphaVdwEff * softcoreLambdaFactorVdw[i] * sigma6[i] + rp);
                                rPInvV = gmx::inv(divisor);
                                sixthRoot(rPInvV, &rInvV, &rV);
                            }
                            else
                            {
                                /* We can avoid one expensive pow and one / operation */
                                rPInvV = rPInvC;
                                rInvV  = rInvC;
                                rV     = rC;
                            }
                        }
                        else
                        {
                            rPInvC = one;
                            rInvC  = rInv;
                            rC     = r;

                            rPInvV = one;
                            rInvV  = rInv;
                            rV     = r;
                        }

                        /* Only process the coulomb interactions if we either
                         * include all entries in the list (no cutoff
                         * used in the kernel), or if we are within the cutoff.
                         */
                        BoolType computeElecInteraction;
                        if constexpr (elecInteractionTypeIsEwald)
                        {
                            computeElecInteraction = (r < rCoulomb && qq[i] != zero && bPairIncluded);
                        }
                        else
                        {
                            computeElecInteraction = (rC < rCoulomb && qq[i] != zero && bPairIncluded);
                        }
                        if (gmx::anyTrue(computeElecInteraction))
                        {
                            if constexpr (elecInteractionTypeIsEwald)
                            {
                                vCoul[i] = ewaldPotential(qq[i], rInvC, sh_ewald);
                                if constexpr (computeScalarForce)
                                {
                                    scalarForcePerDistanceCoul[i] = ewaldScalarForce(qq[i], rInvC);
                                }

                                if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                                {
                                    ewaldQuadraticPotential<computeScalarForce>(
                                            qq[i],
                                            elecEpsilonFactor,
                                            rC,
                                            rCutoffCoul,
                                            lambdaFactorCoul[i],
                                            dLambdaFactor[i],
                                            gapsysScaleLinpointCoulEff,
                                            sh_ewald,
                                            &scalarForcePerDistanceCoul[i],
                                            &vCoul[i],
                                            &dvdlCoul,
                                            computeElecInteraction);
                                }
                            }
                            else
                            {
                                vCoul[i] = reactionFieldPotential(
                                        qq[i], rInvC, rC, reactionFieldCoefficient, reactionFieldShift);
                                if constexpr (computeScalarForce)
                                {
                                    scalarForcePerDistanceCoul[i] = reactionFieldScalarForce(
                                            qq[i], rInvC, rC, reactionFieldCoefficient, two);
                                }

                                if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                                {
                                    reactionFieldQuadraticPotential<computeScalarForce>(
                                            qq[i],
                                            elecEpsilonFactor,
                                            rC,
                                            rCutoffCoul,
                                            lambdaFactorCoul[i],
                                            dLambdaFactor[i],
                                            gapsysScaleLinpointCoulEff,
                                            reactionFieldCoefficient,
                                            reactionFieldShift,
                                            &scalarForcePerDistanceCoul[i],
                                            &vCoul[i],
                                            &dvdlCoul,
                                            computeElecInteraction);
                                }
                            }

                            vCoul[i] = gmx::selectByMask(vCoul[i], computeElecInteraction);
                            if constexpr (computeScalarForce)
                            {
                                scalarForcePerDistanceCoul[i] = gmx::selectByMask(
                                        scalarForcePerDistanceCoul[i], computeElecInteraction);
                            }
                        }

                        /* Only process the VDW interactions if we either
                         * include all entries in the list (no cutoff used
                         * in the kernel), or if we are within the cutoff.
                         */
                        BoolType computeVdwInteraction;
                        if constexpr (ljKernelType == LJKernelType::Ewald)
                        {
                            computeVdwInteraction =
                                    (r < rVdw && (c6[i] != zero || c12[i] != zero) && bPairIncluded);
                        }
                        else
                        {
                            computeVdwInteraction =
                                    (rV < rVdw && (c6[i] != zero || c12[i] != zero) && bPairIncluded);
                        }
                        if (gmx::anyTrue(computeVdwInteraction))
                        {
                            RealType rInv6;
                            if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                            {
                                rInv6 = rPInvV;
                            }
                            else
                            {
                                rInv6 = calculateRinv6(rInvV);
                            }
                            // Avoid overflow at short distance for masked exclusions and
                            // for foreign energy calculations at a hard core end state.
                            // Note that we should limit r^-6, and thus also r^-12, and
                            // not only r^-12, as that could lead to erroneously low instead
                            // of very high foreign energies.
                            rInv6 = gmx::min(rInv6, maxRInvSix);

                            // Notes on soft-core handling of Lennard-Jones interaction:
                            // Beutler soft-core "automatically" works, as we pass in
                            // soft-core modified rV and the force is used for computing
                            // dV/dlambda.
                            // For Gapsys soft-core, we compute the soft-core effect
                            // on plain, potential-shifted LJ. This also works for switched
                            // forces/potentials, as long at the Gapsys soft-core range
                            // is <= rvdw-switch.
                            if constexpr (ljKernelType != LJKernelType::ForceSwitch)
                            {
                                RealType vVdw6  = calculateVdw6(c6[i], rInv6);
                                RealType vVdw12 = calculateVdw12(c12[i], rInv6);

                                vVdw[i] = lennardJonesPotential(
                                        vVdw6, vVdw12, c6[i], c12[i], repulsionShift, dispersionShift, oneSixth, oneTwelfth);
                                if constexpr (computeScalarForce)
                                {
                                    scalarForcePerDistanceVdw[i] = lennardJonesScalarForce(vVdw6, vVdw12);
                                }
                            }
                            else
                            {
                                // LJ force switch
                                RealType rSwitched;
                                RealType rSwitchedSquared;
                                RealType rSwitchedSquaredTimesR;
                                computeForceSwitchVariables(
                                        rV, rVdwSwitch, &rSwitched, &rSwitchedSquared, &rSwitchedSquaredTimesR);

                                vVdw[i] = -c6[i]
                                          * (oneSixth * (rInv6 + dispersionShift)
                                             + forceSwitchPotentialMod(rSwitched,
                                                                       rSwitchedSquared,
                                                                       minusDispersionShift2Div3,
                                                                       minusDispersionShift3Div4,
                                                                       computeVdwInteraction));

                                vVdw[i] =
                                        vVdw[i]
                                        + c12[i]
                                                  * (oneTwelfth * (rInv6 * rInv6 + RealType(repulsionShift))
                                                     + forceSwitchPotentialMod(rSwitched,
                                                                               rSwitchedSquared,
                                                                               minusRepulsionShift2Div3,
                                                                               minusRepulsionShift3Div4,
                                                                               computeVdwInteraction));

                                if constexpr (computeScalarForce)
                                {
                                    scalarForcePerDistanceVdw[i] =
                                            -c6[i]
                                            * forceSwitchScalarForceMod(rInv6,
                                                                        rSwitched,
                                                                        rSwitchedSquaredTimesR,
                                                                        dispersionShift2,
                                                                        dispersionShift3,
                                                                        computeVdwInteraction);

                                    scalarForcePerDistanceVdw[i] =
                                            scalarForcePerDistanceVdw[i]
                                            + c12[i]
                                                      * forceSwitchScalarForceMod(rInv6 * rInv6,
                                                                                  rSwitched,
                                                                                  rSwitchedSquaredTimesR,
                                                                                  repulsionShift2,
                                                                                  repulsionShift3,
                                                                                  computeVdwInteraction);
                                }
                            }

                            if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
                            {
                                lennardJonesQuadraticPotential<computeForces>(
                                        c6[i],
                                        c12[i],
                                        r,
                                        rSq,
                                        lambdaFactorVdw[i],
                                        dLambdaFactor[i],
                                        gapsysSigma6VdWEff[i],
                                        gapsysScaleLinpointVdWEff,
                                        repulsionShift,
                                        dispersionShift,
                                        &scalarForcePerDistanceVdw[i],
                                        &vVdw[i],
                                        &dvdlVdw,
                                        computeVdwInteraction);
                            }

                            if constexpr (ljKernelType == LJKernelType::Ewald)
                            {
                                /* Subtract the grid potential at the cut-off */
                                vVdw[i] = vVdw[i]
                                          + gmx::selectByMask(ewaldLennardJonesGridSubtract(
                                                                      ljPmeC6Grid[i], shLjEwald, oneSixth),
                                                              computeVdwInteraction);
                            }

                            if constexpr (ljKernelType == LJKernelType::PotentialSwitch)
                            {
                                RealType d             = rV - rVdwSwitch;
                                BoolType zeroMask      = zero < d;
                                BoolType potSwitchMask = rV < rVdw;
                                d                      = gmx::selectByMask(d, zeroMask);
                                const RealType d2      = d * d;
                                const RealType sw =
                                        one + d2 * d * (vdw_swV3 + d * (vdw_swV4 + d * vdw_swV5));

                                if constexpr (computeScalarForce)
                                {
                                    const RealType dsw = d2 * (vdw_swF2 + d * (vdw_swF3 + d * vdw_swF4));
                                    scalarForcePerDistanceVdw[i] = potSwitchScalarForceMod(
                                            scalarForcePerDistanceVdw[i], vVdw[i], sw, rV, dsw, potSwitchMask);
                                }
                                vVdw[i] = potSwitchPotentialMod(vVdw[i], sw, potSwitchMask);
                            }

                            vVdw[i] = gmx::selectByMask(vVdw[i], computeVdwInteraction);
                            if constexpr (computeScalarForce)
                            {
                                scalarForcePerDistanceVdw[i] = gmx::selectByMask(
                                        scalarForcePerDistanceVdw[i], computeVdwInteraction);
                            }
                        }

                        if constexpr (computeScalarForce)
                        {
                            /* scalarForcePerDistanceCoul (and scalarForcePerDistanceVdw) now contain: dV/drC * rC
                             * Now we multiply by rC^-6, so it will be: dV/drC * rC^-5
                             * Further down we first multiply by r^4 and then by
                             * the vector r, which in total gives: dV/drC * (r/rC)^-5
                             */
                            scalarForcePerDistanceCoul[i] = scalarForcePerDistanceCoul[i] * rPInvC;
                            scalarForcePerDistanceVdw[i]  = scalarForcePerDistanceVdw[i] * rPInvV;
                        }
                    } // end of block requiring nonZeroState
                }     // end for (int i = 0; i < NSTATES; i++)

                /* Assemble A and B states. */
                BoolType assembleStates = (bPairIncluded && withinCutoffMask);
                if (gmx::anyTrue(assembleStates))
                {
                    for (int i = 0; i < NSTATES; i++)
                    {
                        vCoulTot = vCoulTot + lambdaFactorCoul[i] * vCoul[i];
                        vVdwTot  = vVdwTot + lambdaFactorVdw[i] * vVdw[i];

                        if constexpr (computeForces)
                        {
                            scalarForcePerDistance =
                                    scalarForcePerDistance
                                    + lambdaFactorCoul[i] * scalarForcePerDistanceCoul[i] * rpm2;
                            scalarForcePerDistance =
                                    scalarForcePerDistance
                                    + lambdaFactorVdw[i] * scalarForcePerDistanceVdw[i] * rpm2;
                        }

                        if constexpr (softcoreType == KernelSoftcoreType::Beutler)
                        {
                            dvdlCoul = dvdlCoul + vCoul[i] * dLambdaFactor[i]
                                       + lambdaFactorCoul[i] * alphaCoulEff * softcoreDlFactorCoul[i]
                                                 * scalarForcePerDistanceCoul[i] * sigma6[i];
                            dvdlVdw = dvdlVdw + vVdw[i] * dLambdaFactor[i]
                                      + lambdaFactorVdw[i] * alphaVdwEff * softcoreDlFactorVdw[i]
                                                * scalarForcePerDistanceVdw[i] * sigma6[i];
                        }
                        else
                        {
                            dvdlCoul = dvdlCoul + vCoul[i] * dLambdaFactor[i];
                            dvdlVdw  = dvdlVdw + vVdw[i] * dLambdaFactor[i];
                        }
                    }
                }
            } // end of block requiring bPairIncluded && withinCutoffMask
            /* In the following block bPairIncluded should be false in the masks. */
            if constexpr (!elecInteractionTypeIsEwald)
            {
                if (coulombInteractionType == NbkernelElecType::ReactionField)
                {
                    // With RF do not allow excluded pairs beyond the Coulomb cut-off, check this here.
                    // We'd like to use !withinCutoffMask, but there is no negation operator for SimdFBool.
                    // We need to use <= as this is the exact negation of the cutoff check.
                    const BoolType beyondCutoff = (rCutoffCoul * rCutoffCoul <= rSq);
                    haveExcludedPairsBeyondCutoff =
                            haveExcludedPairsBeyondCutoff || (bPairExcluded && beyondCutoff);

                    const BoolType computeReactionField = bPairExcluded;

                    if (gmx::anyTrue(computeReactionField))
                    {
                        /* For excluded pairs we don't use soft-core.
                         * As there is no singularity, there is no need for soft-core.
                         */
                        const RealType FF = -two * reactionFieldCoefficient;
                        RealType       VV = reactionFieldCoefficient * rSq - reactionFieldShift;

                        /* If ii == jnr the i particle (ii) has itself (jnr)
                         * in its neighborlist. This corresponds to a self-interaction
                         * that will occur twice. Scale it down by 50% to only include
                         * it once.
                         */
                        VV = VV * gmx::blend(one, half, bIiEqJnr);

                        for (int i = 0; i < NSTATES; i++)
                        {
                            vCoulTot = vCoulTot
                                       + gmx::selectByMask(lambdaFactorCoul[i] * qq[i] * VV,
                                                           computeReactionField);
                            scalarForcePerDistance = scalarForcePerDistance
                                                     + gmx::selectByMask(lambdaFactorCoul[i] * qq[i] * FF,
                                                                         computeReactionField);
                            dvdlCoul = dvdlCoul
                                       + gmx::selectByMask(dLambdaFactor[i] * qq[i] * VV,
                                                           computeReactionField);
                        }
                    }
                }
            }

            const BoolType computeElecEwaldInteraction = (bPairExcluded || r < rCoulomb);
            if (elecInteractionTypeIsEwald && gmx::anyTrue(computeElecEwaldInteraction))
            {
                /* See comment in the preamble. When using Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to 1/r (vanilla coulomb short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire electrostatic interaction,
                 * including the reciprocal-space component.
                 */
                RealType v_lr, f_lr;

                pmeCoulombCorrectionVF<computeForces>(rSq, ewaldBeta, &v_lr, &f_lr);
                if constexpr (computeForces)
                {
                    f_lr = f_lr * rInv * rInv;
                }

                /* Note that any possible Ewald shift has already been applied in
                 * the normal interaction part above.
                 */

                /* If ii == jnr the i particle (ii) has itself (jnr)
                 * in its neighborlist. This corresponds to a self-interaction
                 * that will occur twice. Scale it down by 50% to only include
                 * it once.
                 */
                v_lr = v_lr * gmx::blend(one, half, bIiEqJnr);

                for (int i = 0; i < NSTATES; i++)
                {
                    vCoulTot = vCoulTot
                               - gmx::selectByMask(lambdaFactorCoul[i] * qq[i] * v_lr,
                                                   computeElecEwaldInteraction);
                    if constexpr (computeForces)
                    {
                        scalarForcePerDistance = scalarForcePerDistance
                                                 - gmx::selectByMask(lambdaFactorCoul[i] * qq[i] * f_lr,
                                                                     computeElecEwaldInteraction);
                    }
                    dvdlCoul = dvdlCoul
                               - gmx::selectByMask(dLambdaFactor[i] * qq[i] * v_lr,
                                                   computeElecEwaldInteraction);
                }
            }

            const BoolType computeVdwEwaldInteraction = (bPairExcluded || r < rVdw);
            if (ljKernelType == LJKernelType::Ewald && gmx::anyTrue(computeVdwEwaldInteraction))
            {
                /* See comment in the preamble. When using LJ-Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to r^-6 (vanilla LJ6 short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire VdW interaction,
                 * including the reciprocal-space component.
                 */

                RealType v_lr, f_lr;
                pmeLJCorrectionVF<computeForces>(
                        rInv, rSq, ewaldLJCoeffSq, ewaldLJCoeffSixDivSix, &v_lr, &f_lr, computeVdwEwaldInteraction, bIiEqJnr);
                v_lr = v_lr * oneSixth;

                for (int i = 0; i < NSTATES; i++)
                {
                    vVdwTot = vVdwTot
                              + gmx::selectByMask(lambdaFactorVdw[i] * ljPmeC6Grid[i] * v_lr,
                                                  computeVdwEwaldInteraction);
                    if constexpr (computeForces)
                    {
                        scalarForcePerDistance =
                                scalarForcePerDistance
                                + gmx::selectByMask(lambdaFactorVdw[i] * ljPmeC6Grid[i] * f_lr,
                                                    computeVdwEwaldInteraction);
                    }
                    dvdlVdw = dvdlVdw
                              + gmx::selectByMask(dLambdaFactor[i] * ljPmeC6Grid[i] * v_lr,
                                                  computeVdwEwaldInteraction);
                }
            }

            if (computeForces && gmx::anyTrue(scalarForcePerDistance != zero))
            {
                const RealType tX = scalarForcePerDistance * dX;
                const RealType tY = scalarForcePerDistance * dY;
                const RealType tZ = scalarForcePerDistance * dZ;
                fIX               = fIX + tX;
                fIY               = fIY + tY;
                fIZ               = fIZ + tZ;

                gmx::transposeScatterDecrU<3>(forceRealPtr, preloadJnr, tX, tY, tZ);
            }
        } // end for (int k = nj0; k < nj1; k += DataTypes::simdRealWidth)

        if (havePairsWithinCutoff)
        {
            if constexpr (computeForces)
            {
                gmx::transposeScatterIncrU<3>(forceRealPtr, preloadIi, fIX, fIY, fIZ);

                if (doShiftForces)
                {
                    gmx::transposeScatterIncrU<3>(
                            reinterpret_cast<real*>(threadForceShiftBuffer), preloadIs, fIX, fIY, fIZ);
                }
            }
            if (doPotential)
            {
                int ggid = gid[n];
                threadVCoul[ggid] += gmx::reduce(vCoulTot);
                threadVVdw[ggid] += gmx::reduce(vVdwTot);
            }
        }
    } // end for (int n = 0; n < nri; n++)

    if (gmx::anyTrue(dvdlCoul != zero))
    {
        threadDvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)] += gmx::reduce(dvdlCoul);
    }
    if (gmx::anyTrue(dvdlVdw != zero))
    {
        threadDvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)] += gmx::reduce(dvdlVdw);
    }

    /* Estimate flops, average for free energy stuff:
     * 12  flops per outer iteration
     * 150 flops per inner iteration
     * TODO: Update the number of flops and/or use different counts for different code paths.
     */
    atomicNrnbIncrement(nrnb, eNR_NBKERNEL_FREE_ENERGY, nlist.nri * 12 + nlist.jindex[nri] * 150);

    if (coulombInteractionType == NbkernelElecType::ReactionField
        && gmx::anyTrue(haveExcludedPairsBeyondCutoff))
    {
        GMX_THROW(gmx::InvalidInputError(
                "One or more excluded and perturbed atom pairs are beyond the Coulomb cut-off, "
                "which is not allowed with reaction-field."));
    }
}

typedef void (*KernelFunction)(const t_nblist&                                  nlist,
                               const gmx::ArrayRefWithPadding<const gmx::RVec>& coords,
                               const int                                        ntype,
                               const interaction_const_t&          interactionParameters,
                               gmx::ArrayRef<const gmx::RVec>      shiftvec,
                               gmx::ArrayRef<const real>           nbfp,
                               gmx::ArrayRef<const real>           nbfp_grid,
                               gmx::ArrayRef<const real>           chargeA,
                               gmx::ArrayRef<const real>           chargeB,
                               gmx::ArrayRef<const int>            typeA,
                               gmx::ArrayRef<const int>            typeB,
                               int                                 flags,
                               gmx::ArrayRef<const real>           lambda,
                               t_nrnb* gmx_restrict                nrnb,
                               gmx::ArrayRefWithPadding<gmx::RVec> threadForceBuffer,
                               rvec*                               threadForceShiftBuffer,
                               gmx::ArrayRef<real>                 threadVCoul,
                               gmx::ArrayRef<real>                 threadVVdw,
                               gmx::ArrayRef<real>                 threadDvdl);

template<KernelSoftcoreType softcoreType, bool scLambdasOrAlphasDiffer, bool elecInteractionTypeIsEwald, LJKernelType ljKernelType, bool computeForces>
static KernelFunction dispatchKernelOnUseSimd(const bool useSimd)
{
    if (useSimd)
    {
#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS && GMX_USE_SIMD_KERNELS
        return (nb_free_energy_kernel<SimdDataTypes, softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces>);
#else
        return (nb_free_energy_kernel<ScalarDataTypes, softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces>);
#endif
    }
    else
    {
        return (nb_free_energy_kernel<ScalarDataTypes, softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces>);
    }
}

template<KernelSoftcoreType softcoreType, bool scLambdasOrAlphasDiffer, bool elecInteractionTypeIsEwald, LJKernelType ljKernelType>
static KernelFunction dispatchKernelOnComputeForces(const bool computeForces, const bool useSimd)
{
    if (computeForces)
    {
        return (dispatchKernelOnUseSimd<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, true>(
                useSimd));
    }
    else
    {
        return (dispatchKernelOnUseSimd<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, false>(
                useSimd));
    }
}

template<KernelSoftcoreType softcoreType, bool scLambdasOrAlphasDiffer, bool elecInteractionTypeIsEwald>
static KernelFunction dispatchKernelOnLJType(const LJKernelType ljKernelType,
                                             const bool         computeForces,
                                             const bool         useSimd)
{
    switch (ljKernelType)
    {
        case LJKernelType::Cutoff:
            return (dispatchKernelOnComputeForces<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, LJKernelType::Cutoff>(
                    computeForces, useSimd));
            break;
        case LJKernelType::ForceSwitch:
            return (dispatchKernelOnComputeForces<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, LJKernelType::ForceSwitch>(
                    computeForces, useSimd));
            break;
        case LJKernelType::PotentialSwitch:
            return (dispatchKernelOnComputeForces<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, LJKernelType::PotentialSwitch>(
                    computeForces, useSimd));
            break;
        case LJKernelType::Ewald:
            return (dispatchKernelOnComputeForces<softcoreType, scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, LJKernelType::Ewald>(
                    computeForces, useSimd));
            break;
        default: GMX_THROW(gmx::InternalError("Unimplemented LJ kernel type"));
    }
}

template<KernelSoftcoreType softcoreType, bool scLambdasOrAlphasDiffer>
static KernelFunction dispatchKernelOnElecInteractionType(const bool elecInteractionTypeIsEwald,
                                                          const LJKernelType ljKernelType,
                                                          const bool         computeForces,
                                                          const bool         useSimd)
{
    if (elecInteractionTypeIsEwald)
    {
        return (dispatchKernelOnLJType<softcoreType, scLambdasOrAlphasDiffer, true>(
                ljKernelType, computeForces, useSimd));
    }
    else
    {
        return (dispatchKernelOnLJType<softcoreType, scLambdasOrAlphasDiffer, false>(
                ljKernelType, computeForces, useSimd));
    }
}

template<KernelSoftcoreType softcoreType>
static KernelFunction dispatchKernelOnScLambdasOrAlphasDifference(const bool scLambdasOrAlphasDiffer,
                                                                  const bool elecInteractionTypeIsEwald,
                                                                  const LJKernelType ljKernelType,
                                                                  const bool         computeForces,
                                                                  const bool         useSimd)
{
    if (scLambdasOrAlphasDiffer)
    {
        return (dispatchKernelOnElecInteractionType<softcoreType, true>(
                elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
    }
    else
    {
        return (dispatchKernelOnElecInteractionType<softcoreType, false>(
                elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
    }
}

static KernelFunction dispatchKernel(const bool                 scLambdasOrAlphasDiffer,
                                     const bool                 elecInteractionTypeIsEwald,
                                     const LJKernelType         ljKernelType,
                                     const bool                 computeForces,
                                     const bool                 useSimd,
                                     const interaction_const_t& interactionParameters)
{
    const auto& scParams = *interactionParameters.softCoreParameters;
    if (scParams.softcoreType == SoftcoreType::Beutler)
    {
        if (scParams.alphaCoulomb == 0 && scParams.alphaVdw == 0)
        {
            return (dispatchKernelOnScLambdasOrAlphasDifference<KernelSoftcoreType::None>(
                    scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
        }
        return (dispatchKernelOnScLambdasOrAlphasDifference<KernelSoftcoreType::Beutler>(
                scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
    }
    else // Gapsys
    {
        if (scParams.gapsysScaleLinpointCoul == 0 && scParams.gapsysScaleLinpointVdW == 0)
        {
            return (dispatchKernelOnScLambdasOrAlphasDifference<KernelSoftcoreType::None>(
                    scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
        }
        return (dispatchKernelOnScLambdasOrAlphasDifference<KernelSoftcoreType::Gapsys>(
                scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd));
    }
}


void gmx_nb_free_energy_kernel(const t_nblist&                                  nlist,
                               const gmx::ArrayRefWithPadding<const gmx::RVec>& coords,
                               const bool                                       useSimd,
                               const int                                        ntype,
                               const interaction_const_t&          interactionParameters,
                               gmx::ArrayRef<const gmx::RVec>      shiftvec,
                               gmx::ArrayRef<const real>           nbfp,
                               gmx::ArrayRef<const real>           nbfp_grid,
                               gmx::ArrayRef<const real>           chargeA,
                               gmx::ArrayRef<const real>           chargeB,
                               gmx::ArrayRef<const int>            typeA,
                               gmx::ArrayRef<const int>            typeB,
                               int                                 flags,
                               gmx::ArrayRef<const real>           lambda,
                               t_nrnb*                             nrnb,
                               gmx::ArrayRefWithPadding<gmx::RVec> threadForceBuffer,
                               rvec*                               threadForceShiftBuffer,
                               gmx::ArrayRef<real>                 threadVCoul,
                               gmx::ArrayRef<real>                 threadVVdw,
                               gmx::ArrayRef<real>                 threadDvdl)
{
    GMX_ASSERT(usingPmeOrEwald(interactionParameters.eeltype)
                       || interactionParameters.eeltype == CoulombInteractionType::Cut
                       || usingRF(interactionParameters.eeltype),
               "Unsupported eeltype with free energy");
    GMX_ASSERT(interactionParameters.softCoreParameters, "We need soft-core parameters");

    // Not all SIMD implementations need padding, but we provide padding anyhow so we can assert
    GMX_ASSERT(!GMX_SIMD_HAVE_REAL || threadForceBuffer.empty()
                       || threadForceBuffer.size() > threadForceBuffer.unpaddedArrayRef().ssize(),
               "We need actual padding with at least one element for SIMD scatter operations");

    const auto&  scParams                   = *interactionParameters.softCoreParameters;
    const bool   elecInteractionTypeIsEwald = (usingPmeOrEwald(interactionParameters.eeltype));
    LJKernelType ljKernelType = LJKernelType::Cutoff; // Just to make sure it is initialized.
    if (usingLJPme(interactionParameters.vdwtype))
    {
        ljKernelType = LJKernelType::Ewald;

        GMX_ASSERT(interactionParameters.vdw_modifier == InteractionModifiers::PotShift
                           || interactionParameters.vdw_modifier == InteractionModifiers::None,
                   "No force or potential modifiers are supported with LJ-Ewald");
    }
    else if (interactionParameters.vdw_modifier == InteractionModifiers::PotShift
             || interactionParameters.vdw_modifier == InteractionModifiers::None)
    {
        ljKernelType = LJKernelType::Cutoff;
    }
    else if (interactionParameters.vdw_modifier == InteractionModifiers::ForceSwitch)
    {
        ljKernelType = LJKernelType::ForceSwitch;
    }
    else if (interactionParameters.vdw_modifier == InteractionModifiers::PotSwitch)
    {
        ljKernelType = LJKernelType::PotentialSwitch;
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported LJ interaction type");
    }
    const bool computeForces           = ((flags & GMX_NONBONDED_DO_FORCE) != 0);
    bool       scLambdasOrAlphasDiffer = true;

    if (scParams.alphaCoulomb == 0 && scParams.alphaVdw == 0)
    {
        scLambdasOrAlphasDiffer = false;
    }
    else
    {
        if (lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]
                    == lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)]
            && scParams.alphaCoulomb == scParams.alphaVdw)
        {
            scLambdasOrAlphasDiffer = false;
        }
    }

    KernelFunction kernelFunc;
    kernelFunc = dispatchKernel(
            scLambdasOrAlphasDiffer, elecInteractionTypeIsEwald, ljKernelType, computeForces, useSimd, interactionParameters);
    kernelFunc(nlist,
               coords,
               ntype,
               interactionParameters,
               shiftvec,
               nbfp,
               nbfp_grid,
               chargeA,
               chargeB,
               typeA,
               typeB,
               flags,
               lambda,
               nrnb,
               threadForceBuffer,
               threadForceShiftBuffer,
               threadVCoul,
               threadVVdw,
               threadDvdl);
}
