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
 *  HIP FEP foreign energy non-bonded kernels using C++ template-based code generation
 *
 *  Based on Original CUDA version by Yiqi Chen <yiqi.chen@metax-tech.com>
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_HIP_NBFE_FOREIGN_HIP_KERNEL_BODY_H
#define GMX_NBNXM_HIP_NBFE_FOREIGN_HIP_KERNEL_BODY_H

#include <string>
#include <string_view>
#include <type_traits>

#include <rocprim/rocprim.hpp>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gpu_kernel_utils.h"
#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/hip/nbnxm_hip_kernel_utils.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/template_mp.h"

namespace gmx
{

/*! \brief Generate kernel name for FEP foreign energy kernel
 *
 * \tparam elecType     Type of electrostatics
 * \tparam vdwType      Type of VDW interactions
 * \tparam useSoftCore  Whether softcore potentials are in use
 * \tparam pairlistType What pairlist is used, important for execution width
 */
template<ElecType elecType, VdwType vdwType, bool useSoftCore, PairlistType pairlistType>
static const std::string getFepForeignKernelName()
{
    static constexpr std::array<const std::string_view, c_numElecTypes> elecNames = {
        "ElecCut", "ElecRF", "ElecEwQSTab", "ElecEwQSTabTwinCut", "ElecEw", "ElecEwTwinCut", "", ""
    };
    static constexpr std::array<const std::string_view, c_numVdwTypes> vdwNames = {
        "_VdwLJ",    "_VdwLJCombGeom",   "_VdwLJCombLB",  "_VdwLJFsw",
        "_VdwLJPsw", "_VdwLJEwCombGeom", "_VdwLJEwCombLB"
    };

    static constexpr std::string_view softCoreName   = "_softcore";
    static constexpr std::string_view noSoftCoreName = "_nosoftcore";

    constexpr bool                    isWave64 = sc_gpuParallelExecutionWidth(pairlistType) == 64;
    static constexpr std::string_view executionWidth64Name = "_wave64";
    static constexpr std::string_view executionWidth32Name = "_wave32";

    static constexpr std::string_view baseName = "nbfe_foreign_kernel_";
    static constexpr std::string_view elecName = elecNames[static_cast<int>(elecType)];
    static constexpr std::string_view vdwName  = vdwNames[static_cast<int>(vdwType)];
    static constexpr std::string_view type     = "_V";
#if !defined(_MSC_VER)
    return std::string(CompileTimeStringJoin_v<baseName,
                                               elecName,
                                               vdwName,
                                               (useSoftCore ? softCoreName : noSoftCoreName),
                                               (isWave64 ? executionWidth64Name : executionWidth32Name),
                                               type>);
#else
    std::string returnValue;
    returnValue.reserve(1024);
    return returnValue.append(baseName)
            .append(elecName)
            .append(vdwName)
            .append(useSoftCore ? softCoreName : noSoftCoreName)
            .append(isWave64 ? executionWidth64Name : executionWidth32Name)
            .append(type);
#endif
}

/*! \brief FEP foreign energy non-bonded kernel
 *
 * This kernel calculates energies for multiple lambda values to enable
 * efficient calculation of free energy derivatives at foreign lambda points.
 *
 * \tparam elecType        Type of electrostatics
 * \tparam vdwType         Type of VDW interactions
 * \tparam useSoftCore     Whether soft core potentials are in use.
 * \tparam pairlistType    Which specific pairlist is in use, used for reduction
 */
template<ElecType elecType, VdwType vdwType, bool useSoftCore, PairlistType pairlistType>
__global__ void nbfeForeignKernel(const NBAtomDataGpu atdat,
                                  const NBParamGpu    nbparam,
                                  const GpuFeplist    feplist,
                                  const int           nLambda)
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

    // Load atom data pointers
    AmdFastBuffer<const float4> gm_xq{ atdat.xq };
    AmdFastBuffer<const float4> gm_q4{ atdat.q4 };
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

    float beta = 0.0F, v_lr = 0.0F;
    if constexpr (props.elecEwald)
    {
        beta = nbparam.ewald_beta;
    }

    float ewald_shift = 0.0F, c_rf = 0.0F;
    if constexpr (props.elecEwald)
    {
        ewald_shift = nbparam.sh_ewald;
    }
    else
    {
        c_rf = nbparam.c_rf;
    }

    // Thread/block/warp IDs
    const int tid        = threadIdx.y * blockDim.x + threadIdx.x;
    const int tid_global = blockIdx.x * blockDim.x * blockDim.y * blockDim.z
                           + threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x
                           + threadIdx.x;
    const int wid_global  = tid_global / c_parallelExecutionWidth;
    const int tid_in_warp = tid % c_parallelExecutionWidth;
    const int block_size  = blockDim.x * blockDim.y * blockDim.z;
    const int numiAtoms   = feplist.numiAtoms;

    AmdFastBuffer<const int> gm_iinr{ feplist.iinr };
    AmdFastBuffer<const int> gm_jIndex{ feplist.jIndex };
    AmdFastBuffer<const int> gm_jjnr{ feplist.jjnr };
    AmdFastBuffer<const int> gm_shift{ feplist.shift };

    const int       lambdaPower    = nbparam.lambdaPower;
    constexpr float softcoreRPower = 6.0F;

    float dLambdaFactor[2];
    float softcoreDlFactorCoul[2];
    float softcoreDlFactorVdw[2];

    dLambdaFactor[0] = -1.0F;
    dLambdaFactor[1] = 1.0F;

    AmdFastBuffer<const float> gm_allLambdaCoul{ nbparam.allLambdaCoul };
    AmdFastBuffer<const float> gm_allLambdaVdw{ nbparam.allLambdaVdw };

    // Shared memory for lambda values
    extern __shared__ char sm_dynamicSharedMemory[];
    char*                  sm_nextSlotPtr = sm_dynamicSharedMemory;
    float*                 sm_lambdaCoul  = reinterpret_cast<float*>(sm_nextSlotPtr);
    sm_nextSlotPtr += (nLambda + 1) * sizeof(float);
    float* sm_lambdaVdw = reinterpret_cast<float*>(sm_nextSlotPtr);


    for (int idx = tid; idx <= nLambda; idx += block_size)
    {
        if (idx == 0)
        {
            sm_lambdaCoul[0] = nbparam.lambdaCoul;
            sm_lambdaVdw[0]  = nbparam.lambdaVdw;
        }
        else
        {
            sm_lambdaCoul[idx] = gm_allLambdaCoul[idx - 1];
            sm_lambdaVdw[idx]  = gm_allLambdaVdw[idx - 1];
        }
    }

    __syncthreads();

    float* gm_e_lj    = atdat.eLJForeign;
    float* gm_e_el    = atdat.eElecForeign;
    float* gm_dvdl_lj = atdat.dvdlLJForeign;
    float* gm_dvdl_el = atdat.dvdlElecForeign;

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
        const float3 xi = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

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

        for (int i = nj0; i < nj1; i += c_parallelExecutionWidth)
        {
            const int j = i + tid_in_warp;
            if (j >= nj1)
            {
                continue;
            }
            const int  aj           = gm_jjnr[j];
            const bool pairIncluded = (feplist.exclFep == nullptr || (feplist.exclFep[j]));

            const float4 xqj   = gm_xq[aj];
            const float3 xj    = make_float3(xqj.x, xqj.y, xqj.z);
            const float4 q4j   = gm_q4[aj];
            const float  qAj_f = q4j.x;
            const float  qBj_f = q4j.y;
            float        qq[2] = { qAi * qAj_f, qBi * qBj_f };

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

            const float3 rv               = make_float3(xi.x - xj.x, xi.y - xj.y, xi.z - xj.z);
            float        r2               = rv.x * rv.x + rv.y * rv.y + rv.z * rv.z;
            const bool   withinCutoffMask = (r2 < rCutoffMaxSq);
            r2                            = fmaxf(r2, c_nbnxnMinDistanceSquared);
            const float inv_r             = __frsqrt_rn(r2);
            const float inv_r2            = inv_r * inv_r;

            // Loop over all lambda values
            for (int lambdaIdx = 0; lambdaIdx <= nLambda; lambdaIdx++)
            {
                float E_lj    = 0.0F;
                float E_el    = 0.0F;
                float DVDL_lj = 0.0F;
                float DVDL_el = 0.0F;

                const float lambdaCoul         = sm_lambdaCoul[lambdaIdx];
                const float lambdaVdw          = sm_lambdaVdw[lambdaIdx];
                const float oneMinusLambdaCoul = 1.0F - lambdaCoul;
                const float oneMinusLambdaVdw  = 1.0F - lambdaVdw;

                float lambdaFactorCoul[2]         = { oneMinusLambdaCoul, lambdaCoul };
                float lambdaFactorVdw[2]          = { oneMinusLambdaVdw, lambdaVdw };
                float softcoreLambdaFactorCoul[2] = { lambdaCoul, oneMinusLambdaCoul };
                float softcoreLambdaFactorVdw[2]  = { lambdaVdw, oneMinusLambdaVdw };

                for (int k = 0; k < 2; k++)
                {
                    softcoreLambdaFactorCoul[k] =
                            (lambdaPower == 2
                                     ? (1.0F - lambdaFactorCoul[k]) * (1.0F - lambdaFactorCoul[k])
                                     : (1.0F - lambdaFactorCoul[k]));
                    softcoreDlFactorCoul[k] = dLambdaFactor[k] * lambdaPower / softcoreRPower
                                              * (lambdaPower == 2 ? (1.0F - lambdaFactorCoul[k]) : 1.0F);
                    softcoreLambdaFactorVdw[k] =
                            (lambdaPower == 2 ? (1.0F - lambdaFactorVdw[k]) * (1.0F - lambdaFactorVdw[k])
                                              : (1.0F - lambdaFactorVdw[k]));
                    softcoreDlFactorVdw[k] = dLambdaFactor[k] * lambdaPower / softcoreRPower
                                             * (lambdaPower == 2 ? (1.0F - lambdaFactorVdw[k]) : 1.0F);
                }

                float scalarForcePerDistanceVdw[2]  = { 0.0F, 0.0F };
                float scalarForcePerDistanceCoul[2] = { 0.0F, 0.0F };

                if (pairIncluded && withinCutoffMask)
                {
                    float rpm2, rp;
                    if constexpr (useSoftCore)
                    {
                        rpm2 = r2 * r2;
                        rp   = rpm2 * r2;
                    }
                    else
                    {
                        rpm2 = inv_r * inv_r;
                        rp   = 1.0F;
                    }

                    float  sigma6[2];
                    float2 c6c12AB[2];

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

                    float Vcoul[2] = { 0.0F, 0.0F };
                    float Vvdw[2]  = { 0.0F, 0.0F };

                    // Calculate energies for both states
                    for (int k = 0; k < 2; k++)
                    {
                        scalarForcePerDistanceVdw[k]  = 0.0F;
                        scalarForcePerDistanceCoul[k] = 0.0F;

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

                                Vvdw[k] = ((Vvdw12 + c6c12AB[k].y * repulsionShift.cpot) * c_oneTwelfth
                                           - (Vvdw6 + c6c12AB[k].x * dispersionShift.cpot) * c_oneSixth);

                                if constexpr (props.vdwFSwitch)
                                {
                                    ljForceSwitch<true>(dispersionShift,
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
                                    ljPotentialSwitch<true, true>(vdwSwitch,
                                                                  rVdwSwitch,
                                                                  rInvV,
                                                                  r2V,
                                                                  &scalarForcePerDistanceVdw[k],
                                                                  &Vvdw[k]);
                                }

                                if constexpr (props.elecEwaldTwin)
                                {
                                    vdw_in_range = (r2 < rCutoffVdwSq) ? 1.0F : 0.0F;
                                    Vvdw[k] *= vdw_in_range;
                                }
                            }

                            // Electrostatic interactions
                            if (qq[k] != 0.0F)
                            {
                                if constexpr (props.elecCutoff)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
                                    Vcoul[k]                      = qq[k] * (rInvC - c_rf);
                                }
                                else if constexpr (props.elecRF)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * (rInvC - two_k_rf * r2C);
                                    Vcoul[k] = qq[k] * (rInvC + 0.5F * two_k_rf * r2C - c_rf);
                                }
                                else if constexpr (props.elecEwald)
                                {
                                    scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
                                    Vcoul[k]                      = qq[k] * (rInvC - ewald_shift);
                                }
                            }

                            scalarForcePerDistanceCoul[k] *= rPInvC;
                            scalarForcePerDistanceVdw[k] *= rPInvV;
                        }
                    }

                    // Accumulate energies and derivatives
                    for (int k = 0; k < 2; k++)
                    {
                        E_el += lambdaFactorCoul[k] * Vcoul[k];
                        E_lj += lambdaFactorVdw[k] * Vvdw[k];

                        if constexpr (useSoftCore)
                        {
                            DVDL_el += Vcoul[k] * dLambdaFactor[k]
                                       + lambdaFactorCoul[k] * alphaCoulombEff * softcoreDlFactorCoul[k]
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
                }

                // Handle exclusions for RF/Cutoff
                if constexpr (props.elecCutoff || props.elecRF)
                {
                    if (!pairIncluded)
                    {
                        float VV = 0.0F;
                        if constexpr (props.elecCutoff)
                        {
                            VV = -nbparam.c_rf;
                        }
                        else if constexpr (props.elecRF)
                        {
                            VV = 0.5F * two_k_rf * r2 - nbparam.c_rf;
                        }

                        if (ai == aj)
                        {
                            VV *= 0.5F;
                        }
                        for (int k = 0; k < 2; k++)
                        {
                            E_el += lambdaFactorCoul[k] * qq[k] * VV;
                            DVDL_el += (dLambdaFactor[k] * qq[k]) * VV;
                        }
                    }
                }

                // Handle Ewald exclusions
                if constexpr (props.elecEwald)
                {
                    if ((!pairIncluded || r2 < rCutoffCoulSq))
                    {
                        v_lr = inv_r > 0.0F ? inv_r * erff(r2 * inv_r * beta)
                                            : 2.0F * beta * M_FLOAT_1_SQRTPI;
                        if (ai == aj)
                        {
                            v_lr *= 0.5F;
                        }
                        for (int k = 0; k < 2; k++)
                        {
                            E_el -= lambdaFactorCoul[k] * qq[k] * v_lr;
                            DVDL_el -= (dLambdaFactor[k] * qq[k]) * v_lr;
                        }
                    }
                }

                // Reduce energies for this lambda value
                reduceEnergyWarpShuffle<pairlistType>(
                        E_lj, E_el, gm_e_lj + lambdaIdx, gm_e_el + lambdaIdx, tid);
                reduceEnergyWarpShuffle<pairlistType>(
                        DVDL_lj, DVDL_el, gm_dvdl_lj + lambdaIdx, gm_dvdl_el + lambdaIdx, tid);
            }
        }
    }
}

/*! \brief Implementation of foreign kernel launcher helper
 *
 * This function is explicitly instantiated in the compilation units
 * to avoid compiling all kernel variants in a single translation unit.
 */
template<ElecType elecType, VdwType vdwType, bool useSoftCore, PairlistType pairlistType>
void launchNbfeForeignKernel(const NBAtomDataGpu* atdat,
                             const NBParamGpu*    nbparam,
                             const GpuFeplist*    feplist,
                             int                  nLambda,
                             const DeviceContext& deviceContext,
                             const DeviceStream&  deviceStream)
{
    constexpr int      c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);
    KernelLaunchConfig config;
    config.blockSize[0] = 64;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    const int blockSize = config.blockSize[0] * config.blockSize[1] * config.blockSize[2];
    // one wave per nri.
    const int nriPerBlock   = blockSize / c_parallelExecutionWidth;
    const int nblock        = gmx::divideRoundUp(feplist->numiAtoms, nriPerBlock);
    config.gridSize[0]      = numberOfKernelBlocksSanityCheck(nblock, deviceContext.deviceInfo());
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = (nLambda + 1) * 2 * sizeof(float);

    auto kernelPtr  = nbfeForeignKernel<elecType, vdwType, useSoftCore, pairlistType>;
    auto kernelName = getFepForeignKernelName<elecType, vdwType, useSoftCore, pairlistType>();
    auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, atdat, nbparam, feplist, &nLambda);

    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, kernelName.c_str(), kernelArgs);
}

/*! \brief Kernel launcher helper for FEP foreign lambda energy kernels
 */
void launchNbfeForeignKernelHelper(NbnxmGpu* nb, const InteractionLocality iloc)
{
    const NBParamGpu*    nbparam       = nb->nbparam;
    const NBAtomDataGpu* atdat         = nb->atdat;
    const ElecType       elecType      = nbparam->elecType;
    const VdwType        vdwType       = nbparam->vdwType;
    const bool           useSoftCore   = (nbparam->alphaVdw != 0.0F);
    auto*                feplist       = nb->feplist[iloc].get();
    const int            nLambda       = nb->fephostdata->allLambdaCoul.size();
    const DeviceStream&  deviceStream  = *nb->deviceStreams[iloc];
    const DeviceContext& deviceContext = *nb->deviceContext_;

    std::visit(
            [&](auto&& pairlists)
            {
                // only needed to get type
                auto* plist                 = pairlists[iloc].get();
                using T                     = std::decay_t<decltype(*plist)>;
                constexpr auto pairlistType = getPairlistTypeFromPairlist<T>();

                gmx::dispatchTemplatedFunction(
                        [&](auto elecType_, auto vdwType_, auto useSoftCore_)
                        {
                            launchNbfeForeignKernel<elecType_, vdwType_, useSoftCore_, pairlistType>(
                                    atdat, nbparam, feplist, nLambda, deviceContext, deviceStream);
                        },
                        elecType,
                        vdwType,
                        useSoftCore);
            },
            nb->plist);
}

} // namespace gmx

#endif // GMX_NBNXM_HIP_NBFE_FOREIGN_HIP_KERNEL_BODY_H
