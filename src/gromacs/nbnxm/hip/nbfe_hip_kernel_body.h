/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 *  \brief
 *  HIP FEP non-bonded kernels using C++ template-based code generation
 *
 *  Based on Original CUDA version by Yiqi Chen <yiqi.chen@metax-tech.com>
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_HIP_NBFE_HIP_KERNELS_H
#define GMX_NBNXM_HIP_NBFE_HIP_KERNELS_H

#include <string>
#include <string_view>
#include <type_traits>

#include <rocprim/rocprim.hpp>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_kernel_utils.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_hip_kernel_utils.h"

namespace gmx
{


/*! \brief Generate kernel name for FEP kernel
 *
 * \tparam elecType     Type of electrostatics
 * \tparam vdwType      Type of VDW interactions
 * \tparam doCalcEnergies Whether to calculate energies
 * \tparam doCalcVirial  Whether to calculate shift contributions
 * \tparam useSoftCore   Whether soft core potentials are used
 * \tparam pairlistType    Which pair list is used, for parallel execution width
 */
template<ElecType elecType, VdwType vdwType, bool doCalcEnergies, bool doCalcVirial, bool useSoftCore, PairlistType pairlistType>
static const std::string getFepKernelName()
{
    static constexpr std::array<const std::string_view, c_numElecTypes> elecNames = {
        "ElecCut", "ElecRF", "ElecEwQSTab", "ElecEwQSTabTwinCut", "ElecEw", "ElecEwTwinCut", "", ""
    };
    static constexpr std::array<const std::string_view, c_numVdwTypes> vdwNames = {
        "_VdwLJ",    "_VdwLJCombGeom",   "_VdwLJCombLB",  "_VdwLJFsw",
        "_VdwLJPsw", "_VdwLJEwCombGeom", "_VdwLJEwCombLB"
    };

    static constexpr std::string_view baseName         = "nbfe";
    static constexpr std::string_view elecName         = elecNames[static_cast<int>(elecType)];
    static constexpr std::string_view vdwName          = vdwNames[static_cast<int>(vdwType)];
    static constexpr std::string_view calcEnergyName   = "_VF";
    static constexpr std::string_view noCalcEnergyName = "_F";
    static constexpr std::string_view calcVirialName   = "_shift";
    static constexpr std::string_view noCalcVirialName = "";
    static constexpr std::string_view softCoreName     = "_softcore";
    static constexpr std::string_view noSoftCoreName   = "_nosoftcore";

    constexpr bool                    isWave64 = sc_gpuParallelExecutionWidth(pairlistType) == 64;
    static constexpr std::string_view executionWidth64Name = "_wave64";
    static constexpr std::string_view executionWidth32Name = "_wave32";


#if !defined(_MSC_VER)
    return std::string(CompileTimeStringJoin_v<baseName,
                                               elecName,
                                               vdwName,
                                               (doCalcEnergies ? calcEnergyName : noCalcEnergyName),
                                               (doCalcVirial ? calcVirialName : noCalcVirialName),
                                               (useSoftCore ? softCoreName : noSoftCoreName),
                                               (isWave64 ? executionWidth64Name : executionWidth32Name)>);
#else
    std::string returnValue;
    returnValue.reserve(1024);
    return returnValue.append(baseName)
            .append(elecName)
            .append(vdwName)
            .append(doCalcEnergies ? calcEnergyName : noCalcEnergyName)
            .append(doCalcVirial ? calcVirialName : noCalcVirialName)
            .append(useSoftCore ? softCoreName : noSoftCoreName)
            .append(isWave64 ? executionWidth64Name : executionWidth32Name);
#endif
}

/*! \brief HIP FEP nonbonded kernel
 *
 * \tparam doCalcEnergies  Whether to calculate energies
 * \tparam elecType        Type of electrostatics
 * \tparam vdwType         Type of VDW interactions
 * \tparam doCalcEnergies  Whether to calculate energies
 * \tparam doCalcVirial    Whether to calculate shift contributions
 * \tparam useSoftCore     Whether soft core potentials are used
 * \tparam pairlistType    Which pair list is used, for parallel execution width
 */
template<ElecType elecType, VdwType vdwType, bool doCalcEnergies, bool doCalcVirial, bool useSoftCore, PairlistType pairlistType>
__global__ void nbfeKernel(const NBAtomDataGpu atdat, const NBParamGpu nbparam, const GpuFeplist feplist)
{
    constexpr EnergyFunctionProperties<elecType, vdwType> props;
    constexpr int c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);

    // Convenience variables
    const float alphaCoulomb           = nbparam.alphaCoul;
    const float alphaVdw               = nbparam.alphaVdw;
    float       alphaCoulombEff        = alphaCoulomb;
    float       alphaVdwEff            = alphaVdw;
    const float sigma6WithInvalidSigma = nbparam.sigma6WithInvalidSigma;
    const float sigma6Minimum          = nbparam.sigma6Minimum;
    const float lambdaCoul             = nbparam.lambdaCoul;
    const float oneMinusLambdaCoul     = 1.0F - lambdaCoul;
    const float lambdaVdw              = nbparam.lambdaVdw;
    const float oneMinusLambdaVdw      = 1.0F - lambdaVdw;

    float lambdaFactorCoul[2]         = { oneMinusLambdaCoul, lambdaCoul };
    float lambdaFactorVdw[2]          = { oneMinusLambdaVdw, lambdaVdw };
    float softcoreLambdaFactorCoul[2] = { lambdaCoul, oneMinusLambdaCoul };
    float softcoreLambdaFactorVdw[2]  = { lambdaVdw, oneMinusLambdaVdw };

    // Load atom data pointers
    AmdFastBuffer<const float4> gm_xq{ atdat.xq };
    AmdFastBuffer<const float4> gm_q4{ atdat.q4 };
    float3*                     gm_f = asFloat3(atdat.f);
    AmdFastBuffer<const float3> gm_shiftVec{ asFloat3(atdat.shiftVec) };

    float                 rCutoffCoulSq   = nbparam.rcoulomb_sq;
    float                 rCutoffMaxSq    = rCutoffCoulSq;
    const shift_consts_t  dispersionShift = nbparam.dispersion_shift;
    const shift_consts_t  repulsionShift  = nbparam.repulsion_shift;
    const switch_consts_t vdwSwitch       = nbparam.vdw_switch;
    const float           rVdwSwitch      = nbparam.rvdw_switch;
    const float           ewaldShift      = nbparam.sh_ewald;


    float rCutoffVdwSq = 0.0F;
    float vdw_in_range = 1.0F;
    if constexpr (props.elecEwaldTwin)
    {
        rCutoffVdwSq = nbparam.rvdw_sq;
        rCutoffMaxSq = fmaxf(rCutoffCoulSq, rCutoffVdwSq);
    }

    float two_k_rf = 0.0F;
    if constexpr (props.elecRF)
    {
        two_k_rf = nbparam.two_k_rf;
    }

    float beta2 = 0.0F, beta3 = 0.0F;
    if constexpr (props.elecEwaldAna)
    {
        beta2 = nbparam.ewald_beta * nbparam.ewald_beta;
        beta3 = nbparam.ewald_beta * beta2;
    }

    float beta = 0.0F, v_lr = 0.0F, f_lr = 0.0F;
    if constexpr (props.elecEwald)
    {
        beta = nbparam.ewald_beta;
    }

    float ewald_shift = 0.0F, c_rf = 0.0F;
    if constexpr (doCalcEnergies)
    {
        if constexpr (props.elecEwald)
        {
            ewald_shift = nbparam.sh_ewald;
        }
        else
        {
            c_rf = nbparam.c_rf;
        }
    }

    // Thread/block/warp IDs
    const int tid        = threadIdx.y * blockDim.x + threadIdx.x;
    const int tid_global = blockIdx.x * blockDim.x * blockDim.y * blockDim.z
                           + threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x
                           + threadIdx.x;
    const int                wid_global  = tid_global / c_parallelExecutionWidth;
    const int                tid_in_warp = tid % c_parallelExecutionWidth;
    const int                numiAtoms   = feplist.numiAtoms;
    AmdFastBuffer<const int> gm_iinr{ feplist.iinr };
    AmdFastBuffer<const int> gm_jIndex{ feplist.jIndex };
    AmdFastBuffer<const int> gm_jjnr{ feplist.jjnr };
    AmdFastBuffer<const int> gm_shift{ feplist.shift };

    float* gm_e_lj     = nullptr;
    float* gm_e_el     = nullptr;
    float* gm_dvdl_lj  = nullptr;
    float* gm_dvdl_el  = nullptr;
    int    lambdaPower = 0;
    float  dLambdaFactor[2];
    float  softcoreDlFactorCoul[2];
    float  softcoreDlFactorVdw[2];

    if constexpr (doCalcEnergies)
    {
        gm_e_lj     = atdat.eLJ;
        gm_e_el     = atdat.eElec;
        gm_dvdl_lj  = atdat.dvdlLJ;
        gm_dvdl_el  = atdat.dvdlElec;
        lambdaPower = nbparam.lambdaPower;

        dLambdaFactor[0] = -1.0F;
        dLambdaFactor[1] = 1.0F;

        constexpr float softcoreRPower = 6.0F;

        for (int k = 0; k < 2; k++)
        {
            softcoreLambdaFactorCoul[k] =
                    (lambdaPower == 2 ? (1.0F - lambdaFactorCoul[k]) * (1.0F - lambdaFactorCoul[k])
                                      : (1.0F - lambdaFactorCoul[k]));
            softcoreDlFactorCoul[k] = dLambdaFactor[k] * lambdaPower / softcoreRPower
                                      * (lambdaPower == 2 ? (1.0F - lambdaFactorCoul[k]) : 1.0F);
            softcoreLambdaFactorVdw[k] =
                    (lambdaPower == 2 ? (1.0F - lambdaFactorVdw[k]) * (1.0F - lambdaFactorVdw[k])
                                      : (1.0F - lambdaFactorVdw[k]));
            softcoreDlFactorVdw[k] = dLambdaFactor[k] * lambdaPower / softcoreRPower
                                     * (lambdaPower == 2 ? (1.0F - lambdaFactorVdw[k]) : 1.0F);
        }
    }

    float E_lj = 0.0F, E_el = 0.0F;
    float DVDL_lj = 0.0F, DVDL_el = 0.0F;
    // Accumulate energies and derivatives
    float scalarForcePerDistance = 0.0F;

    if (wid_global < numiAtoms)
    {
        // Only thread 0 reads from memory, then broadcasts to avoid out-of-bounds accesses
        int nj0_tmp      = (tid_in_warp == 0) ? gm_jIndex[wid_global] : 0;
        int nj1_tmp      = (tid_in_warp == 0) ? gm_jIndex[wid_global + 1] : 0;
        int ai_tmp       = (tid_in_warp == 0) ? gm_iinr[wid_global] : 0;
        int shiftIdx_tmp = (tid_in_warp == 0) ? gm_shift[wid_global] : 0;

        const int nj0      = __shfl(nj0_tmp, 0, c_parallelExecutionWidth);
        const int nj1      = __shfl(nj1_tmp, 0, c_parallelExecutionWidth);
        const int ai       = __shfl(ai_tmp, 0, c_parallelExecutionWidth);
        const int shiftIdx = __shfl(shiftIdx_tmp, 0, c_parallelExecutionWidth);

        // Now all threads read from the same (valid) indices
        const float3 shiftI = gm_shiftVec[shiftIdx];

        float4 xqbuf = gm_xq[ai];
        xqbuf.x += shiftI.x;
        xqbuf.y += shiftI.y;
        xqbuf.z += shiftI.z;
        const AmdPackedFloat3 xi = { xqbuf.x, xqbuf.y, xqbuf.z };

        const float4 q4_buf = gm_q4[ai];
        const float  qAi    = q4_buf.x * nbparam.epsfac;
        const float  qBi    = q4_buf.y * nbparam.epsfac;

        int4 atomTypes4_buf;
        int  typeiAB[2];
        if constexpr (!props.vdwComb)
        {
            atomTypes4_buf = atdat.atomTypes4[ai];
            typeiAB[0]     = atomTypes4_buf.x;
            typeiAB[1]     = atomTypes4_buf.y;
        }

        float4 ljComb4_buf;
        float2 ljcp_iAB[2];
        if constexpr (props.vdwComb)
        {
            ljComb4_buf = atdat.ljComb4[ai];
            ljcp_iAB[0] = make_float2(ljComb4_buf.x, ljComb4_buf.y);
            ljcp_iAB[1] = make_float2(ljComb4_buf.z, ljComb4_buf.w);
        }

        AmdPackedFloat3 fci = { 0.0F, 0.0F, 0.0F };

        for (int i = nj0; i < nj1; i += c_parallelExecutionWidth)
        {
            const int j = i + tid_in_warp;
            if (j < nj1)
            {
                scalarForcePerDistance = 0.0F;
                const int aj           = gm_jjnr[j];
                bool      pairIncluded = (feplist.exclFep == nullptr || feplist.exclFep[j]);

                const float4          xqj   = gm_xq[aj];
                const AmdPackedFloat3 xj    = { xqj.x, xqj.y, xqj.z };
                const float4          q4j   = gm_q4[aj];
                const float           qAj_f = q4j.x;
                const float           qBj_f = q4j.y;
                float                 qq[2] = { qAi * qAj_f, qBi * qBj_f };

                int typejAB[2];
                if constexpr (!props.vdwComb)
                {
                    const int4 atomTypesj = atdat.atomTypes4[aj];
                    typejAB[0]            = atomTypesj.x;
                    typejAB[1]            = atomTypesj.y;
                }

                float2 ljcp_jAB[2];
                if constexpr (props.vdwComb)
                {
                    const float4 ljCombj = atdat.ljComb4[aj];
                    ljcp_jAB[0]          = make_float2(ljCombj.x, ljCombj.y);
                    ljcp_jAB[1]          = make_float2(ljCombj.z, ljCombj.w);
                }

                const AmdPackedFloat3 rv = xi - xj;
                float                 r2 = rv.norm2();

                const bool withinCutoffMask = (r2 < rCutoffMaxSq);
                if (!withinCutoffMask && pairIncluded)
                {
                    continue;
                }

                r2                 = fmaxf(r2, c_nbnxnMinDistanceSquared);
                const float inv_r  = __frsqrt_rn(r2);
                const float inv_r2 = inv_r * inv_r;

                float rpm2, rp;
                if (pairIncluded)
                {
                    if constexpr (useSoftCore)
                    {
                        rpm2 = r2 * r2;   // r^4
                        rp   = rpm2 * r2; // r^6
                    }
                    else
                    {
                        rpm2 = inv_r * inv_r;
                        rp   = 1.0F;
                    }

                    float  sigma6[2];
                    float2 c6c12AB[2];
                    float  scalarForcePerDistanceCoul[2] = { 0.0F, 0.0F };
                    float  scalarForcePerDistanceVdw[2]  = { 0.0F, 0.0F };
                    float  Vcoul[2]                      = { 0.0F, 0.0F };
                    float  Vvdw[2]                       = { 0.0F, 0.0F };

                    // Get LJ parameters for both states
                    for (int k = 0; k < 2; k++)
                    {
                        if constexpr (!props.vdwComb)
                        {
                            const int numTypes = atdat.numTypes;
                            c6c12AB[k] =
                                    fetchNbfpC6C12(nbparam.nbfp, numTypes * typeiAB[k] + typejAB[k]);
                            if constexpr (useSoftCore)
                            {
                                sigma6[k] = convertC6C12ToSigma6(
                                        c6c12AB[k], sigma6Minimum, sigma6WithInvalidSigma);
                            }
                        }
                        else if constexpr (props.vdwCombGeom)
                        {
                            c6c12AB[k].x = ljcp_iAB[k].x * ljcp_jAB[k].x;
                            c6c12AB[k].y = ljcp_iAB[k].y * ljcp_jAB[k].y;
                            if constexpr (useSoftCore)
                            {
                                sigma6[k] = convertC6C12ToSigma6(
                                        c6c12AB[k], sigma6Minimum, sigma6WithInvalidSigma);
                            }
                        }
                        else if constexpr (props.vdwCombLB)
                        {
                            float sigmaAB = ljcp_iAB[k].x + ljcp_jAB[k].x;
                            if (ljcp_iAB[k].x == 0.0F || ljcp_jAB[k].x == 0.0F)
                            {
                                sigmaAB = 0.0F;
                            }
                            const float epsilonAB = ljcp_iAB[k].y * ljcp_jAB[k].y;
                            c6c12AB[k]            = convertSigmaEpsilonToC6C12(sigmaAB, epsilonAB);
                            if constexpr (useSoftCore)
                            {
                                if ((c6c12AB[k].x > 0.0F) && (c6c12AB[k].y > 0.0F))
                                {
                                    const float sigma2 = sigmaAB * sigmaAB;
                                    sigma6[k] = fmaxf(sigma2 * sigma2 * sigma2 * 0.5F, sigma6Minimum);
                                }
                                else
                                {
                                    sigma6[k] = sigma6WithInvalidSigma;
                                }
                            }
                        }
                    }

                    // Determine effective softcore parameters
                    if constexpr (useSoftCore)
                    {
                        if ((c6c12AB[0].y > 0.0F) && (c6c12AB[1].y > 0.0F))
                        {
                            alphaVdwEff     = 0.0F;
                            alphaCoulombEff = 0.0F;
                        }
                        else
                        {
                            alphaVdwEff     = alphaVdw;
                            alphaCoulombEff = alphaCoulomb;
                        }
                    }

                    // Calculate forces and energies for both states
                    for (int k = 0; k < 2; k++)
                    {
                        const bool nonZeroState =
                                ((qq[k] != 0.0F) || (c6c12AB[k].x != 0.0F) || (c6c12AB[k].y != 0.0F));
                        if (nonZeroState)
                        {
                            float rPInvC, r2C, rInvC;
                            float rPInvV, r2V, rInvV;

                            if constexpr (useSoftCore)
                            {
                                rPInvC = 1.0F
                                         / (alphaCoulombEff * softcoreLambdaFactorCoul[k] * sigma6[k] + rp);
                                // rcbrt(x) = 1/cbrt(x) = x^(-1/3), so we use 1.0F/cbrtf to match CUDA
                                r2C   = 1.0F / cbrtf(rPInvC);
                                rInvC = __frsqrt_rn(r2C);

                                if ((alphaCoulombEff != alphaVdwEff)
                                    || (softcoreLambdaFactorVdw[k] != softcoreLambdaFactorCoul[k]))
                                {
                                    rPInvV = 1.0F
                                             / (alphaVdwEff * softcoreLambdaFactorVdw[k] * sigma6[k] + rp);
                                    r2V   = 1.0F / cbrtf(rPInvV);
                                    rInvV = __frsqrt_rn(r2V);
                                }
                                else
                                {
                                    rPInvV = rPInvC;
                                    r2V    = r2C;
                                    rInvV  = rInvC;
                                }
                            }
                            else
                            {
                                rPInvC = 1.0F;
                                r2C    = r2;
                                rInvC  = inv_r;
                                rPInvV = 1.0F;
                                r2V    = r2;
                                rInvV  = inv_r;
                            }

                            // VDW interactions
                            if (c6c12AB[k].x != 0.0F || c6c12AB[k].y != 0.0F)
                            {
                                float rInv6;
                                if constexpr (!useSoftCore)
                                {
                                    rInv6 = inv_r2 * inv_r2 * inv_r2;
                                }
                                else
                                {
                                    rInv6 = rPInvV;
                                }
                                const float Vvdw6            = c6c12AB[k].x * rInv6;
                                const float Vvdw12           = c6c12AB[k].y * rInv6 * rInv6;
                                scalarForcePerDistanceVdw[k] = Vvdw12 - Vvdw6;

                                if constexpr (doCalcEnergies || props.vdwPSwitch)
                                {
                                    Vvdw[k] = ((Vvdw12 + c6c12AB[k].y * repulsionShift.cpot) * c_oneTwelfth
                                               - (Vvdw6 + c6c12AB[k].x * dispersionShift.cpot) * c_oneSixth);
                                }

                                if constexpr (props.vdwFSwitch)
                                {
                                    ljForceSwitch<doCalcEnergies>(dispersionShift,
                                                                  repulsionShift,
                                                                  c6c12AB[k],
                                                                  rVdwSwitch,
                                                                  rInvV,
                                                                  r2V,
                                                                  &scalarForcePerDistanceVdw[k],
                                                                  &Vvdw[k]);
                                }

                                if constexpr (props.vdwPSwitch)
                                {
                                    ljPotentialSwitch<doCalcEnergies, true>(vdwSwitch,
                                                                            rVdwSwitch,
                                                                            rInvV,
                                                                            r2V,
                                                                            &scalarForcePerDistanceVdw[k],
                                                                            &Vvdw[k]);
                                }

                                if constexpr (props.elecEwaldTwin)
                                {
                                    vdw_in_range = (r2 < rCutoffVdwSq) ? 1.0F : 0.0F;
                                    scalarForcePerDistanceVdw[k] *= vdw_in_range;
                                    if constexpr (doCalcEnergies)
                                    {
                                        Vvdw[k] *= vdw_in_range;
                                    }
                                }
                            }

                            // Electrostatic interactions
                            if (qq[k] != 0.0F)
                            {
                                if constexpr (props.elecCutoff)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
                                    if constexpr (doCalcEnergies)
                                    {
                                        Vcoul[k] = qq[k] * (rInvC - c_rf);
                                    }
                                }
                                else if constexpr (props.elecRF)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * (rInvC - two_k_rf * r2C);
                                    if constexpr (doCalcEnergies)
                                    {
                                        Vcoul[k] = qq[k] * (rInvC + 0.5F * two_k_rf * r2C - c_rf);
                                    }
                                }
                                else if constexpr (props.elecEwald)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
                                    if constexpr (doCalcEnergies)
                                    {
                                        Vcoul[k] = qq[k] * (rInvC - ewald_shift);
                                    }
                                }
                            }

                            scalarForcePerDistanceCoul[k] *= rPInvC;
                            scalarForcePerDistanceVdw[k] *= rPInvV;
                        }
                    }

                    for (int k = 0; k < 2; k++)
                    {
                        if constexpr (doCalcEnergies)
                        {
                            E_el += lambdaFactorCoul[k] * Vcoul[k];
                            E_lj += lambdaFactorVdw[k] * Vvdw[k];

                            if constexpr (useSoftCore)
                            {
                                DVDL_el += Vcoul[k] * dLambdaFactor[k]
                                           + lambdaFactorCoul[k] * alphaCoulombEff
                                                     * softcoreDlFactorCoul[k]
                                                     * scalarForcePerDistanceCoul[k] * sigma6[k];
                                DVDL_lj += Vvdw[k] * dLambdaFactor[k]
                                           + lambdaFactorVdw[k] * alphaVdwEff * softcoreDlFactorVdw[k]
                                                     * scalarForcePerDistanceVdw[k] * sigma6[k];
                            }
                            else
                            {
                                DVDL_el += Vcoul[k] * dLambdaFactor[k];
                                DVDL_lj += Vvdw[k] * dLambdaFactor[k];
                            }
                        }
                        scalarForcePerDistance +=
                                lambdaFactorCoul[k] * scalarForcePerDistanceCoul[k] * rpm2;
                        scalarForcePerDistance += lambdaFactorVdw[k] * scalarForcePerDistanceVdw[k] * rpm2;
                    }
                }

                // Handle exclusions for RF/Cutoff
                if constexpr (props.elecCutoff || props.elecRF)
                {
                    if (!pairIncluded)
                    {
                        float FF = 0.0F, VV = 0.0F;
                        if constexpr (props.elecCutoff)
                        {
                            FF = 0.0F;
                            VV = -nbparam.c_rf;
                        }
                        else if constexpr (props.elecRF)
                        {
                            FF = -two_k_rf;
                            VV = 0.5F * two_k_rf * r2 - nbparam.c_rf;
                        }

                        if (ai == aj)
                        {
                            VV *= 0.5F;
                        }
                        for (int k = 0; k < 2; k++)
                        {
                            if constexpr (doCalcEnergies)
                            {
                                E_el += lambdaFactorCoul[k] * qq[k] * VV;
                                DVDL_el += dLambdaFactor[k] * qq[k] * VV;
                            }
                            scalarForcePerDistance += lambdaFactorCoul[k] * qq[k] * FF;
                        }
                    }
                }

                // Handle Ewald exclusions
                if constexpr (props.elecEwald)
                {
                    if (!pairIncluded || r2 < rCutoffCoulSq)
                    {
                        v_lr = inv_r > 0.0F ? inv_r * erff(r2 * inv_r * beta)
                                            : 2.0F * beta * M_FLOAT_1_SQRTPI;
                        if (ai == aj)
                        {
                            v_lr *= 0.5F;
                        }

                        if constexpr (props.elecEwaldAna)
                        {
                            f_lr = inv_r > 0.0F ? -pmeCorrF(beta2 * r2) * beta3 : 0.0F;
                        }
                        else if constexpr (props.elecEwaldTab)
                        {
                            f_lr = inv_r > 0.0F ? interpolateCoulombForceR(nbparam.coulomb_tab,
                                                                           nbparam.coulomb_tab_scale,
                                                                           r2 * inv_r)
                                                          * inv_r
                                                : 0.0F;
                        }

                        for (int k = 0; k < 2; k++)
                        {
                            if constexpr (doCalcEnergies)
                            {
                                E_el -= lambdaFactorCoul[k] * qq[k] * v_lr;
                                DVDL_el -= dLambdaFactor[k] * qq[k] * v_lr;
                            }
                            scalarForcePerDistance -= lambdaFactorCoul[k] * qq[k] * f_lr;
                        }
                    }
                }

                // Accumulate forces
                if (scalarForcePerDistance != 0.0F)
                {
                    const AmdPackedFloat3 f_ij = rv * scalarForcePerDistance;

                    float3 fcj = make_float3(-f_ij[0], -f_ij[1], -f_ij[2]);
                    fci += f_ij;
                    staggeredAtomicAddForce(gm_f, fcj, aj, tid);
                }
            }
        }

        // Reduce i forces
        reduceFepForceWarpShuffle<pairlistType>(fci, gm_f, tid, ai);
        if (doCalcVirial && (gm_shift[wid_global] != gmx::c_centralShiftIndex))
        {
            float3* fShift = asFloat3(atdat.fShift);
            staggeredAtomicAddForce(fShift, make_float3(fci[0], fci[1], fci[2]), gm_shift[wid_global], tid);
        }

        // Reduce energies
        if constexpr (doCalcEnergies)
        {
            reduceEnergyWarpShuffle<pairlistType>(E_lj, E_el, gm_e_lj, gm_e_el, tid);
            reduceEnergyWarpShuffle<pairlistType>(DVDL_lj, DVDL_el, gm_dvdl_lj, gm_dvdl_el, tid);
        }
    }
}

/*! \brief Implementation of kernel launcher helper
 */
template<ElecType elecType, VdwType vdwType, bool doCalcEnergies, bool doCalcVirial, bool useSoftCore, PairlistType pairlistType, class... Args>
void launchNbfeKernel(const GpuFeplist*    feplist,
                      const DeviceContext& deviceContext,
                      const DeviceStream&  deviceStream,
                      Args*... args)
{
    constexpr int c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);

    KernelLaunchConfig config;
    config.blockSize[0] = 64;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;

    const int blockSize     = config.blockSize[0] * config.blockSize[1] * config.blockSize[2];
    const int nriPerBlock   = blockSize / c_parallelExecutionWidth;
    const int nblock        = gmx::divideRoundUp(feplist->numiAtoms, nriPerBlock);
    config.gridSize[0]      = numberOfKernelBlocksSanityCheck(nblock, deviceContext.deviceInfo());
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0; // FEP kernels don't use shared memory for i-atom data

    auto kernelPtr =
            nbfeKernel<elecType, vdwType, doCalcEnergies, doCalcVirial, useSoftCore, pairlistType>;
    auto kernelName =
            getFepKernelName<elecType, vdwType, doCalcEnergies, doCalcVirial, useSoftCore, pairlistType>();
    // Use variadic args pattern to match non-FEP kernel launch
    auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, args...);

    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, kernelName.c_str(), kernelArgs);
}

/*! \brief Kernel launcher helper for FEP foreign lambda energy kernels
 */
template<bool doCalcEnergies, bool doCalcVirial>
void launchNbfeKernelHelper(NbnxmGpu* nb, const InteractionLocality iloc)
{
    const NBParamGpu*    nbparam       = nb->nbparam;
    const NBAtomDataGpu* atdat         = nb->atdat;
    const ElecType       elecType      = nbparam->elecType;
    const VdwType        vdwType       = nbparam->vdwType;
    const bool           useSoftCore   = (nbparam->alphaVdw != 0.0F);
    auto*                feplist       = nb->feplist[iloc].get();
    const DeviceStream&  deviceStream  = *nb->deviceStreams[iloc];
    const DeviceContext& deviceContext = *nb->deviceContext_;
    std::visit(
            [&](auto&& pairlists)
            {
                auto* plist                 = pairlists[iloc].get();
                using T                     = std::decay_t<decltype(*plist)>;
                constexpr auto pairlistType = getPairlistTypeFromPairlist<T>();

                gmx::dispatchTemplatedFunction(
                        [&](auto elecType_, auto vdwType_, auto useSoftCore_)
                        {
                            launchNbfeKernel<elecType_, vdwType_, doCalcEnergies, doCalcVirial, useSoftCore_, pairlistType>(
                                    feplist, deviceContext, deviceStream, atdat, nbparam, feplist);
                        },
                        elecType,
                        vdwType,
                        useSoftCore);
            },
            nb->plist);
}

template<bool doCalcEnergies, bool doCalcVirial>
void launchNbnxmFepKernelHelper(NbnxmGpu* nb, const InteractionLocality iloc)
{
    launchNbfeKernelHelper<doCalcEnergies, doCalcVirial>(nb, iloc);
}

} // namespace gmx

#endif // GMX_NBNXM_HIP_NBFE_HIP_KERNELS_H
