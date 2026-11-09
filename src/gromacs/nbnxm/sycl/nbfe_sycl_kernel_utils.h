/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Pair-interaction maths shared by the SYCL nonbonded free-energy kernels
 *
 * The routines here are the SYCL counterparts of the bodies of the CUDA kernels
 * in \c nbnxm/cuda/nbfe_cuda_kernel.cuh and \c nbnxm/cuda/nbfe_foreign_cuda_kernel.cuh.
 * They are factored out so that the current-lambda (force + energy) kernel and the
 * foreign-lambda (energy-only) kernel share exactly the same interaction maths.
 *
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_UTILS_H
#define GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_UTILS_H

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/nbnxm/nbnxm_kernel_utils.h"

namespace gmx
{

//! \brief The exponent of r in the soft-core interaction, matching the CPU kernels.
static constexpr float sc_softcoreRPower = 6.0F;

//! \brief Derivative of the lambda factor w.r.t. lambda for the two end states.
static constexpr float sc_dLambdaFactor[2] = { -1.0F, 1.0F };

//! \brief Reciprocal cube root, i.e. \f$ x^{-1/3} \f$.
static inline float inverseCbrt(const float x)
{
    return 1.0F / sycl::cbrt(x);
}

/*! \brief Convert LJ C6,C12 parameters to sigma^6, clamped from below by \p sigma6Min.
 *
 * \param c6         The C6 parameter.
 * \param c12        The C12 parameter.
 * \param sigma6Min  Lower bound on the returned sigma^6.
 * \param sigma6Def  Value returned when sigma cannot be computed from \p c6 and \p c12.
 */
static inline float convertC6C12ToSigma6(const float c6, const float c12, const float sigma6Min, const float sigma6Def)
{
    if ((c6 > 0.0F) && (c12 > 0.0F))
    {
        return sycl::max(0.5F * c12 / c6, sigma6Min);
    }
    return sigma6Def;
}

/*! \brief Lambda-dependent prefactors used by the free-energy kernels.
 *
 * Index 0 refers to end state A, index 1 to end state B.
 */
struct NbfeLambdaFactors
{
    //! Lambda factor for the electrostatic interactions
    float lambdaFactorCoul[2];
    //! Lambda factor for the Van der Waals interactions
    float lambdaFactorVdw[2];
    //! Soft-core lambda factor for the electrostatic interactions
    float softcoreLambdaFactorCoul[2];
    //! Soft-core lambda factor for the Van der Waals interactions
    float softcoreLambdaFactorVdw[2];
    //! Derivative of \c softcoreLambdaFactorCoul w.r.t. lambda, scaled by 1/rPower
    float softcoreDlFactorCoul[2];
    //! Derivative of \c softcoreLambdaFactorVdw w.r.t. lambda, scaled by 1/rPower
    float softcoreDlFactorVdw[2];
};

/*! \brief Compute the lambda-dependent prefactors for one pair of lambda values.
 *
 * \param lambdaCoul  Lambda for the electrostatic interactions.
 * \param lambdaVdw   Lambda for the Van der Waals interactions.
 * \param lambdaPower Exponent of the soft-core lambda dependence (1 or 2).
 */
static inline NbfeLambdaFactors makeNbfeLambdaFactors(const float lambdaCoul,
                                                      const float lambdaVdw,
                                                      const int   lambdaPower)
{
    NbfeLambdaFactors lf;
    lf.lambdaFactorCoul[0] = 1.0F - lambdaCoul;
    lf.lambdaFactorCoul[1] = lambdaCoul;
    lf.lambdaFactorVdw[0]  = 1.0F - lambdaVdw;
    lf.lambdaFactorVdw[1]  = lambdaVdw;

    for (int k = 0; k < 2; k++)
    {
        const float oneMinusCoul = 1.0F - lf.lambdaFactorCoul[k];
        const float oneMinusVdw  = 1.0F - lf.lambdaFactorVdw[k];

        lf.softcoreLambdaFactorCoul[k] = (lambdaPower == 2) ? oneMinusCoul * oneMinusCoul : oneMinusCoul;
        lf.softcoreLambdaFactorVdw[k] = (lambdaPower == 2) ? oneMinusVdw * oneMinusVdw : oneMinusVdw;
        lf.softcoreDlFactorCoul[k] = sc_dLambdaFactor[k] * lambdaPower / sc_softcoreRPower
                                     * ((lambdaPower == 2) ? oneMinusCoul : 1.0F);
        lf.softcoreDlFactorVdw[k] = sc_dLambdaFactor[k] * lambdaPower / sc_softcoreRPower
                                    * ((lambdaPower == 2) ? oneMinusVdw : 1.0F);
    }
    return lf;
}

/*! \brief Lambda-independent interaction parameters for the free-energy kernels. */
struct NbfeInteractionParams
{
    //! Cut-off squared for the electrostatic interactions
    float rCoulombSq;
    //! Cut-off squared for the Van der Waals interactions
    float rVdwSq;
    //! Twice the reaction-field coefficient
    float twoKRf;
    //! Reaction-field potential shift
    float cRF;
    //! Ewald splitting coefficient
    float ewaldBeta;
    //! Ewald potential shift
    float ewaldShift;
    //! Dispersion shift constants
    shift_consts_t dispersionShift;
    //! Repulsion shift constants
    shift_consts_t repulsionShift;
    //! Potential-switch constants
    switch_consts_t vdwSwitch;
    //! Start of the Van der Waals switching region
    float rVdwSwitch;
    //! Scale of the tabulated Ewald force
    float coulombTabScale;
    //! Soft-core alpha for the electrostatic interactions
    float alphaCoulomb;
    //! Soft-core alpha for the Van der Waals interactions
    float alphaVdw;
    //! Value used for sigma^6 when it cannot be derived from C6 and C12
    float sigma6WithInvalidSigma;
    //! Lower bound on sigma^6
    float sigma6Minimum;
};

//! \brief Collect the lambda-independent interaction parameters from \p nbp.
static inline NbfeInteractionParams makeNbfeInteractionParams(const NBParamGpu& nbp)
{
    NbfeInteractionParams params;
    params.rCoulombSq             = nbp.rcoulomb_sq;
    params.rVdwSq                 = nbp.rvdw_sq;
    params.twoKRf                 = nbp.two_k_rf;
    params.cRF                    = nbp.c_rf;
    params.ewaldBeta              = nbp.ewald_beta;
    params.ewaldShift             = nbp.sh_ewald;
    params.dispersionShift        = nbp.dispersion_shift;
    params.repulsionShift         = nbp.repulsion_shift;
    params.vdwSwitch              = nbp.vdw_switch;
    params.rVdwSwitch             = nbp.rvdw_switch;
    params.coulombTabScale        = nbp.coulomb_tab_scale;
    params.alphaCoulomb           = nbp.alphaCoul;
    params.alphaVdw               = nbp.alphaVdw;
    params.sigma6WithInvalidSigma = nbp.sigma6WithInvalidSigma;
    params.sigma6Minimum          = nbp.sigma6Minimum;
    return params;
}

//! \brief Lambda-independent per-pair data for the free-energy kernels.
struct NbfePairData
{
    //! Distance squared, clamped from below to avoid overflow of r^-12
    float r2;
    //! Inverse distance
    float rInv;
    //! Inverse distance squared
    float r2Inv;
    //! Product of the i- and j-charges (times epsfac) for both end states
    float qq[2];
    //! C6 for both end states
    float c6AB[2];
    //! C12 for both end states
    float c12AB[2];
    //! sigma^6 for both end states; only meaningful with the soft-core interaction
    float sigma6[2];
    //! Whether this pair is not a perturbed exclusion
    bool pairIncluded;
    //! Whether this pair is within the largest cut-off
    bool withinCutoff;
    //! Whether the i- and j-atom are the same atom (a self-interaction)
    bool selfPair;
    //! Effective soft-core alpha for the electrostatic interactions for this pair
    float alphaCoulombEff;
    //! Effective soft-core alpha for the Van der Waals interactions for this pair
    float alphaVdwEff;
};

//! Output accumulators of the free-energy pair interaction.
struct NbfeAccumulators
{
    //! Scalar force divided by the distance, i.e. F/r
    float scalarForcePerDistance = 0.0F;
    //! Electrostatic energy
    float energyElec = 0.0F;
    //! Van der Waals energy
    float energyVdw = 0.0F;
    //! Derivative of the electrostatic energy w.r.t. lambda
    float dvdlElec = 0.0F;
    //! Derivative of the Van der Waals energy w.r.t. lambda
    float dvdlVdw = 0.0F;
};

/*! \brief Fill in the C6, C12 and sigma^6 values of \p pair for both end states.
 *
 * \tparam useSoftCore Whether the soft-core interaction is used.
 * \tparam vdwType Van der Waals interaction type, selecting the parameter source.
 * \param  params  Lambda-independent interaction parameters.
 * \param  gm_nbfp Non-bonded parameter table; used iff \c vdwType has no combination rule.
 * \param  numTypes Number of atom types in \p gm_nbfp.
 * \param  typeiAB A- and B-state type of the i-atom; used iff \c vdwType has no combination rule.
 * \param  typejAB A- and B-state type of the j-atom; used iff \c vdwType has no combination rule.
 * \param  ljcpiAB A- and B-state LJ combination parameters of the i-atom.
 * \param  ljcpjAB A- and B-state LJ combination parameters of the j-atom.
 * \param  pair    Pair data to fill in.
 */
template<bool useSoftCore, enum VdwType vdwType>
static inline void nbfeComputeLjParams(const NbfeInteractionParams& params,
                                       const Float2* __restrict__ gm_nbfp,
                                       const int     numTypes,
                                       const int     typeiAB[2],
                                       const int     typejAB[2],
                                       const Float2  ljcpiAB[2],
                                       const Float2  ljcpjAB[2],
                                       NbfePairData* pair)
{
    static constexpr EnergyFunctionProperties<ElecType::Count, vdwType> props;

    for (int k = 0; k < 2; k++)
    {
        if constexpr (!props.vdwComb)
        {
            const Float2 c6c12 = gm_nbfp[numTypes * typeiAB[k] + typejAB[k]];
            pair->c6AB[k]      = c6c12[0];
            pair->c12AB[k]     = c6c12[1];
            if constexpr (useSoftCore)
            {
                pair->sigma6[k] = convertC6C12ToSigma6(
                        pair->c6AB[k], pair->c12AB[k], params.sigma6Minimum, params.sigma6WithInvalidSigma);
            }
        }
        else if constexpr (props.vdwCombGeom)
        {
            pair->c6AB[k]  = ljcpiAB[k][0] * ljcpjAB[k][0];
            pair->c12AB[k] = ljcpiAB[k][1] * ljcpjAB[k][1];
            if constexpr (useSoftCore)
            {
                pair->sigma6[k] = convertC6C12ToSigma6(
                        pair->c6AB[k], pair->c12AB[k], params.sigma6Minimum, params.sigma6WithInvalidSigma);
            }
        }
        else
        {
            float sigmaAB = ljcpiAB[k][0] + ljcpjAB[k][0];
            if (ljcpiAB[k][0] == 0.0F || ljcpjAB[k][0] == 0.0F)
            {
                sigmaAB = 0.0F;
            }
            const float  epsilonAB = ljcpiAB[k][1] * ljcpjAB[k][1];
            const Float2 c6c12     = convertSigmaEpsilonToC6C12(sigmaAB, epsilonAB);
            pair->c6AB[k]          = c6c12[0];
            pair->c12AB[k]         = c6c12[1];
            if constexpr (useSoftCore)
            {
                if ((pair->c6AB[k] > 0.0F) && (pair->c12AB[k] > 0.0F))
                {
                    const float sigma2 = sigmaAB * sigmaAB;
                    pair->sigma6[k] = sycl::max(sigma2 * sigma2 * sigma2 * 0.5F, params.sigma6Minimum);
                }
                else
                {
                    pair->sigma6[k] = params.sigma6WithInvalidSigma;
                }
            }
        }
    }

    if constexpr (useSoftCore)
    {
        // Only use softcore if one of the states has a zero endstate
        if ((pair->c12AB[0] > 0.0F) && (pair->c12AB[1] > 0.0F))
        {
            pair->alphaVdwEff     = 0.0F;
            pair->alphaCoulombEff = 0.0F;
        }
        else
        {
            pair->alphaVdwEff     = params.alphaVdw;
            pair->alphaCoulombEff = params.alphaCoulomb;
        }
    }
}

/*! \brief Accumulate the free-energy interaction of one atom pair at one lambda.
 *
 * \tparam doCalcEnergies Whether energies and dV/dlambda are accumulated.
 * \tparam useSoftCore    Whether the soft-core interaction is used.
 * \tparam elecType       Electrostatics interaction type.
 * \tparam vdwType        Van der Waals interaction type.
 * \param  params         Lambda-independent interaction parameters.
 * \param  lf             Lambda-dependent prefactors.
 * \param  pair           Lambda-independent pair data.
 * \param  gm_coulombTab  Tabulated Ewald force; used iff \c elecType is tabulated Ewald.
 * \param  acc            Accumulators to add the contributions of this pair to.
 */
template<bool doCalcEnergies, bool useSoftCore, enum ElecType elecType, enum VdwType vdwType>
static inline void nbfeAccumulatePair(const NbfeInteractionParams& params,
                                      const NbfeLambdaFactors&     lf,
                                      const NbfePairData&          pair,
                                      const float* __restrict__ gm_coulombTab,
                                      NbfeAccumulators* acc)
{
    static constexpr EnergyFunctionProperties<elecType, vdwType> props;

    /* Note that, unlike in the main NBNXM kernel, the LJ-Ewald grid correction is not
     * applied here: the perturbed long-range LJ part is handled outside the kernel.
     * Hence VdwType::EwaldGeom and VdwType::EwaldLB behave like VdwType::Cut. */
    static_assert(!props.vdwEwald || (!props.vdwComb && !props.vdwFSwitch && !props.vdwPSwitch));

    //! Whether the LJ potential is needed even when not calculating energies
    constexpr bool needVdwPotential = doCalcEnergies || props.vdwPSwitch;

    const float r2    = pair.r2;
    const float rInv  = pair.rInv;
    const float r2Inv = pair.r2Inv;

    if (pair.pairIncluded && pair.withinCutoff)
    {
        const float rpm2 = useSoftCore ? r2 * r2 : r2Inv;
        const float rp   = useSoftCore ? r2 * r2 * r2 : 1.0F;

        float scalarForcePerDistanceCoul[2] = { 0.0F, 0.0F };
        float scalarForcePerDistanceVdw[2]  = { 0.0F, 0.0F };
        float vCoul[2]                      = { 0.0F, 0.0F };
        float vVdw[2]                       = { 0.0F, 0.0F };

        for (int k = 0; k < 2; k++)
        {
            const bool nonZeroState =
                    ((pair.qq[k] != 0.0F) || (pair.c6AB[k] != 0.0F) || (pair.c12AB[k] != 0.0F));
            if (!nonZeroState)
            {
                continue;
            }

            float rPInvC;
            float rPInvV;
            float r2C;
            float rInvC;
            float r2V   = r2;
            float rInvV = rInv;
            if constexpr (useSoftCore)
            {
                rPInvC = 1.0F
                         / (pair.alphaCoulombEff * lf.softcoreLambdaFactorCoul[k] * pair.sigma6[k] + rp);
                r2C   = inverseCbrt(rPInvC);
                rInvC = sycl::rsqrt(r2C);

                // Equivalent to scLambdasOrAlphasDiffer
                if ((pair.alphaCoulombEff != pair.alphaVdwEff)
                    || (lf.softcoreLambdaFactorVdw[k] != lf.softcoreLambdaFactorCoul[k]))
                {
                    rPInvV = 1.0F / (pair.alphaVdwEff * lf.softcoreLambdaFactorVdw[k] * pair.sigma6[k] + rp);
                    if constexpr (props.vdwFSwitch || props.vdwPSwitch)
                    {
                        r2V   = inverseCbrt(rPInvV);
                        rInvV = sycl::rsqrt(r2V);
                    }
                }
                else
                {
                    rPInvV = rPInvC;
                    if constexpr (props.vdwFSwitch || props.vdwPSwitch)
                    {
                        r2V   = r2C;
                        rInvV = rInvC;
                    }
                }
            }
            else
            {
                rPInvC = 1.0F;
                r2C    = r2;
                rInvC  = rInv;
                rPInvV = 1.0F;
            }

            if (pair.c6AB[k] != 0.0F || pair.c12AB[k] != 0.0F)
            {
                const float rInv6  = useSoftCore ? rPInvV : r2Inv * r2Inv * r2Inv;
                const float vVdw6  = pair.c6AB[k] * rInv6;
                const float vVdw12 = pair.c12AB[k] * rInv6 * rInv6;

                scalarForcePerDistanceVdw[k] = vVdw12 - vVdw6;
                if constexpr (needVdwPotential)
                {
                    vVdw[k] = ((vVdw12 + pair.c12AB[k] * params.repulsionShift.cpot) * c_oneTwelfth
                               - (vVdw6 + pair.c6AB[k] * params.dispersionShift.cpot) * c_oneSixth);
                }

                if constexpr (props.vdwFSwitch)
                {
                    ljForceSwitch<doCalcEnergies, true>(params.dispersionShift,
                                                        params.repulsionShift,
                                                        params.rVdwSwitch,
                                                        pair.c6AB[k],
                                                        pair.c12AB[k],
                                                        rInvV,
                                                        r2V,
                                                        &scalarForcePerDistanceVdw[k],
                                                        &vVdw[k]);
                }
                if constexpr (props.vdwPSwitch)
                {
                    ljPotentialSwitch<doCalcEnergies, true>(params.vdwSwitch,
                                                            params.rVdwSwitch,
                                                            rInvV,
                                                            r2V,
                                                            &scalarForcePerDistanceVdw[k],
                                                            &vVdw[k]);
                }

                if constexpr (props.elecEwaldTwin)
                {
                    const float vdwInRange = (r2 < params.rVdwSq) ? 1.0F : 0.0F;
                    scalarForcePerDistanceVdw[k] *= vdwInRange;
                    if constexpr (needVdwPotential)
                    {
                        vVdw[k] *= vdwInRange;
                    }
                }
            }

            if (pair.qq[k] != 0.0F)
            {
                if constexpr (props.elecRF)
                {
                    scalarForcePerDistanceCoul[k] = pair.qq[k] * (rInvC - params.twoKRf * r2C);
                }
                else
                {
                    scalarForcePerDistanceCoul[k] = pair.qq[k] * rInvC;
                }
                if constexpr (doCalcEnergies)
                {
                    if constexpr (props.elecCutoff)
                    {
                        vCoul[k] = pair.qq[k] * (rInvC - params.cRF);
                    }
                    if constexpr (props.elecRF)
                    {
                        vCoul[k] = pair.qq[k] * (rInvC + 0.5F * params.twoKRf * r2C - params.cRF);
                    }
                    if constexpr (props.elecEwald)
                    {
                        vCoul[k] = pair.qq[k] * (rInvC - params.ewaldShift);
                    }
                }
            }
            scalarForcePerDistanceCoul[k] *= rPInvC;
            scalarForcePerDistanceVdw[k] *= rPInvV;
        } // for k

        for (int k = 0; k < 2; k++)
        {
            if constexpr (doCalcEnergies)
            {
                acc->energyElec += lf.lambdaFactorCoul[k] * vCoul[k];
                acc->energyVdw += lf.lambdaFactorVdw[k] * vVdw[k];

                if constexpr (useSoftCore)
                {
                    acc->dvdlElec += vCoul[k] * sc_dLambdaFactor[k]
                                     + lf.lambdaFactorCoul[k] * pair.alphaCoulombEff
                                               * lf.softcoreDlFactorCoul[k]
                                               * scalarForcePerDistanceCoul[k] * pair.sigma6[k];
                    acc->dvdlVdw += vVdw[k] * sc_dLambdaFactor[k]
                                    + lf.lambdaFactorVdw[k] * pair.alphaVdwEff
                                              * lf.softcoreDlFactorVdw[k]
                                              * scalarForcePerDistanceVdw[k] * pair.sigma6[k];
                }
                else
                {
                    acc->dvdlElec += vCoul[k] * sc_dLambdaFactor[k];
                    acc->dvdlVdw += vVdw[k] * sc_dLambdaFactor[k];
                }
            }
            acc->scalarForcePerDistance += lf.lambdaFactorCoul[k] * scalarForcePerDistanceCoul[k] * rpm2;
            acc->scalarForcePerDistance += lf.lambdaFactorVdw[k] * scalarForcePerDistanceVdw[k] * rpm2;
        }
    } // if (pairIncluded && withinCutoff)

    // Exclusion correction for plain cut-off and reaction-field electrostatics
    if constexpr (props.elecCutoff || props.elecRF)
    {
        if (!pair.pairIncluded)
        {
            float ffExcl;
            float vvExcl;
            if constexpr (props.elecCutoff)
            {
                ffExcl = 0.0F;
                vvExcl = -params.cRF;
            }
            else
            {
                ffExcl = -params.twoKRf;
                vvExcl = 0.5F * params.twoKRf * r2 - params.cRF;
            }
            if (pair.selfPair)
            {
                vvExcl *= 0.5F;
            }
            for (int k = 0; k < 2; k++)
            {
                if constexpr (doCalcEnergies)
                {
                    acc->energyElec += lf.lambdaFactorCoul[k] * pair.qq[k] * vvExcl;
                    acc->dvdlElec += sc_dLambdaFactor[k] * pair.qq[k] * vvExcl;
                }
                acc->scalarForcePerDistance += lf.lambdaFactorCoul[k] * pair.qq[k] * ffExcl;
            }
        }
    }

    // Subtract the long-range part of the Ewald interaction
    if constexpr (props.elecEwald)
    {
        if (!pair.pairIncluded || r2 < params.rCoulombSq)
        {
            const float beta2 = params.ewaldBeta * params.ewaldBeta;
            float       vLr   = rInv > 0.0F ? rInv * sycl::erf(r2 * rInv * params.ewaldBeta)
                                            : 2.0F * params.ewaldBeta * c_oneOverSqrtPi;
            if (pair.selfPair)
            {
                vLr *= 0.5F;
            }
            float fLr;
            if constexpr (props.elecEwaldAna)
            {
                fLr = -pmeCorrF(beta2 * r2) * beta2 * params.ewaldBeta;
            }
            else
            {
                fLr = interpolateCoulombForceR(gm_coulombTab, params.coulombTabScale, r2 * rInv) * rInv;
            }
            for (int k = 0; k < 2; k++)
            {
                if constexpr (doCalcEnergies)
                {
                    acc->energyElec -= lf.lambdaFactorCoul[k] * pair.qq[k] * vLr;
                    acc->dvdlElec -= sc_dLambdaFactor[k] * pair.qq[k] * vLr;
                }
                acc->scalarForcePerDistance -= lf.lambdaFactorCoul[k] * pair.qq[k] * fLr;
            }
        }
    }
}

} // namespace gmx

#endif // GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_UTILS_H
