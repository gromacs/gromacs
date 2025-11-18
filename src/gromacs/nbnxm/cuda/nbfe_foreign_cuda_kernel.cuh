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
 *  CUDA FEP foreign energy non-bonded kernel used through preprocessor-based code generation
 *  of multiple kernel flavors, see nbfe_foreign_cuda_kernels.cuh.
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 *
 *  \author Yiqi Chen <yiqi.chen@metax-tech.com>
 *
 */

#include <cub/cub.cuh>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/pbcutil/ishift.h"

namespace gmx
{
/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#if defined EL_EWALD_ANA || defined EL_EWALD_TAB
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#    define EL_EWALD_ANY
#endif

#if defined EL_EWALD_ANY || defined EL_RF || defined LJ_EWALD \
        || (defined EL_CUTOFF && defined CALC_ENERGIES)
/* Macro to control the calculation of exclusion forces in the kernel
 * We do that with Ewald (elec/vdw) and RF. Cut-off only has exclusion
 * energy terms.
 *
 * Note: convenience macro, needs to be undef-ed at the end of the file.
 */
#    define EXCLUSION_FORCES
#endif

#if defined LJ_EWALD_COMB_GEOM || defined LJ_EWALD_COMB_LB
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#    define LJ_EWALD
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
#    define LJ_COMB
#endif

/*
    Each warp is in charge of one entry in numiAtoms.

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */


__global__ void NB_FOREIGN_FEP_KERNEL_FUNC_NAME(nbfe_foreign_kernel, _V_cuda)(const NBAtomDataGpu atdat,
                                                                              const NBParamGpu nbparam,
                                                                              const GpuFeplist feplist,
                                                                              const int nLambda)
#ifdef FUNCTION_DECLARATION_ONLY
        ; /* Only do function declaration, omit the function body. */
#else
{
    /* convenience variables */
    const float alphaCoulomb           = nbparam.alphaCoul;
    const float alphaVdw               = nbparam.alphaVdw;
    float       alphaCoulombEff        = alphaCoulomb;
    float       alphaVdwEff            = alphaVdw;
    const bool  useSoftCore            = (alphaVdw != 0.0F);
    const float sigma6WithInvalidSigma = nbparam.sigma6WithInvalidSigma;
    const float sigma6Minimum          = nbparam.sigma6Minimum;

#    ifndef LJ_COMB
    const int4* gm_atomTypes4 = atdat.atomTypes4;
    int         numTypes      = atdat.numTypes;
    int         typeiAB[2], typejAB[2];
#    else
    const float4* gm_ljComb4 = atdat.ljComb4;
#    endif

    float rInvC, r2C, rPInvC, rPInvV, rInv6;
#    if defined LJ_POT_SWITCH || defined LJ_FORCE_SWITCH
    float rInvV, r2V;
#    endif
    float sigma6[2], c6AB[2], c12AB[2];
    float qq[2];
    float scalarForcePerDistanceCoul[2], scalarForcePerDistanceVdw[2];

    float Vcoul[2];
    float Vvdw[2];

#    ifdef LJ_COMB_LB
    float sigmaAB[2];
    float epsilonAB[2];
#    endif

    const float4* gm_xq       = atdat.xq;
    const float4* gm_q4       = atdat.q4;
    const float3* gm_shiftVec = asFloat3(atdat.shiftVec);

    float rCutoffCoulSq = nbparam.rcoulomb_sq;
    float rCutoffMaxSq  = rCutoffCoulSq;
#    ifdef VDW_CUTOFF_CHECK
    float rCutoffVdwSq = nbparam.rvdw_sq;
    float vdw_in_range;
    rCutoffMaxSq = max(rCutoffCoulSq, rCutoffVdwSq);
#    endif
#    ifdef EL_RF
    float two_k_rf = nbparam.two_k_rf;
#    endif

#    ifdef EL_EWALD_ANY
    float beta = nbparam.ewald_beta;
    float v_lr;
#    endif

// #    ifdef CALC_ENERGIES
#    ifdef EL_EWALD_ANY
    float ewald_shift = nbparam.sh_ewald;
#    else
    float c_rf = nbparam.c_rf;
#    endif /* EL_EWALD_ANY */

    /* thread/block/warp id-s */
    int tid        = threadIdx.y * blockDim.x + threadIdx.x;
    int tid_global = blockIdx.x * blockDim.x * blockDim.y * blockDim.z
                     + threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    int wid_global = tid_global / warp_size;
    // thread Id within a warp
    int tid_in_warp = tid % warp_size;

#    ifndef LJ_COMB
#    else
    float2 ljcp_iAB[2], ljcp_jAB[2];
#    endif
    float qAi, qAj_f, r2, rpm2, rp, inv_r, inv_r2;
    float qBi, qBj_f;

    // float  int_bit;
    float E_lj, E_el;
    float DVDL_lj, DVDL_el;

    float4 xqbuf, q4_buf;
#    ifndef LJ_COMB
    int4 atomTypes4_buf;
#    else
    float4 ljComb4_buf;
#    endif

    float3 xi, xj, rv;

    // Extract pair list data
    const int  numiAtoms = feplist.numiAtoms;
    const int* gm_iinr   = feplist.iinr;
    const int* gm_jIndex = feplist.jIndex;
    const int* gm_jjnr   = feplist.jjnr;
    const int* gm_shift  = feplist.shift;

    const int       lambdaPower    = nbparam.lambdaPower;
    constexpr float softcoreRPower = 6.0F;

    float dLambdaFactor[2];
    float softcoreDlFactorCoul[2];
    float softcoreDlFactorVdw[2];

    /*derivative of the lambda factor for state A and B */
    dLambdaFactor[0] = -1.0F;
    dLambdaFactor[1] = 1.0F;

    const float* gm_allLambdaCoul = nbparam.allLambdaCoul;
    const float* gm_allLambdaVdw  = nbparam.allLambdaVdw;
    float        lambdaCoul, oneMinusLambdaCoul, lambdaVdw, oneMinusLambdaVdw;
    float        lambdaFactorCoul[2], lambdaFactorVdw[2], softcoreLambdaFactorCoul[2],
            softcoreLambdaFactorVdw[2];

    extern __shared__ float sm_dynamicShmem[];
    float*                  sm_lambdaCoul = sm_dynamicShmem;
    float*                  sm_lambdaVdw  = sm_lambdaCoul + nLambda + 1;

    if (tid == 0)
    {
        sm_lambdaCoul[0] = nbparam.lambdaCoul;
        sm_lambdaVdw[0]  = nbparam.lambdaVdw;
    }
    else if (tid <= nLambda)
    {
        sm_lambdaCoul[tid] = gm_allLambdaCoul[tid - 1];
        sm_lambdaVdw[tid]  = gm_allLambdaVdw[tid - 1];
    }

    __syncthreads();

    float* gm_e_lj    = atdat.eLJForeign;
    float* gm_e_el    = atdat.eElecForeign;
    float* gm_dvdl_lj = atdat.dvdlLJForeign;
    float* gm_dvdl_el = atdat.dvdlElecForeign;

    // Each warp calculates one ri
    if (wid_global < numiAtoms)
    {
        const int nj0 = __shfl_sync(c_fullWarpMask, gm_jIndex[wid_global], 0, warp_size);
        const int nj1 = __shfl_sync(c_fullWarpMask, gm_jIndex[wid_global + 1], 0, warp_size);
        const int ai  = __shfl_sync(c_fullWarpMask, gm_iinr[wid_global], 0, warp_size);

        const float3 shiftI =
                cub::ShuffleIndex<warp_size>(gm_shiftVec[gm_shift[wid_global]], 0, c_fullWarpMask);

        xqbuf = cub::ShuffleIndex<warp_size>(gm_xq[ai], 0, c_fullWarpMask);
        xqbuf = xqbuf + make_float4(shiftI.x, shiftI.y, shiftI.z, 0.0F);
        xi    = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

        q4_buf = gm_q4[ai];
        qAi    = q4_buf.x * nbparam.epsfac;
        qBi    = q4_buf.y * nbparam.epsfac;

#    ifndef LJ_COMB
        atomTypes4_buf = gm_atomTypes4[ai];
        typeiAB[0]     = atomTypes4_buf.x;
        typeiAB[1]     = atomTypes4_buf.y;
#    else
        ljComb4_buf = gm_ljComb4[ai];
        ljcp_iAB[0] = make_float2(ljComb4_buf.x, ljComb4_buf.y);
        ljcp_iAB[1] = make_float2(ljComb4_buf.z, ljComb4_buf.w);
#    endif

        int  aj;
        bool pairIncluded;

        for (int i = nj0; i < nj1; i += warp_size)
        {
            int j        = i + tid_in_warp;
            aj           = gm_jjnr[j];
            pairIncluded = (feplist.exclFep == nullptr || (feplist.exclFep[j] && j < nj1));
            /* load j atom data */
            xqbuf = gm_xq[aj];
            xj    = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

            q4_buf = gm_q4[aj];
            qAj_f  = q4_buf.x;
            qBj_f  = q4_buf.y;

            qq[0] = qAi * qAj_f;
            qq[1] = qBi * qBj_f;
#    ifndef LJ_COMB
            atomTypes4_buf = gm_atomTypes4[aj];
            typejAB[0]     = atomTypes4_buf.x;
            typejAB[1]     = atomTypes4_buf.y;
#    else
            ljComb4_buf = gm_ljComb4[aj];
            ljcp_jAB[0] = make_float2(ljComb4_buf.x, ljComb4_buf.y);
            ljcp_jAB[1] = make_float2(ljComb4_buf.z, ljComb4_buf.w);
#    endif

            /* distance between i and j atoms */
            rv = xi - xj;
            r2 = gmxDeviceNorm2(rv);

            bool withinCutoffMask = (r2 < rCutoffMaxSq);
            // Ensure distance do not become so small that r^-12 overflows
            r2     = max(r2, c_nbnxnMinDistanceSquared);
            inv_r  = rsqrt(r2);
            inv_r2 = inv_r * inv_r;

            for (int lambdaIdx = 0; lambdaIdx <= nLambda; lambdaIdx++)
            {
                E_lj    = 0.0F;
                E_el    = 0.0F;
                DVDL_lj = 0.0F;
                DVDL_el = 0.0F;

                lambdaCoul         = sm_lambdaCoul[lambdaIdx];
                lambdaVdw          = sm_lambdaVdw[lambdaIdx];
                oneMinusLambdaCoul = 1.0F - lambdaCoul;
                oneMinusLambdaVdw  = 1.0F - lambdaVdw;

                lambdaFactorCoul[0]         = oneMinusLambdaCoul;
                lambdaFactorCoul[1]         = lambdaCoul;
                lambdaFactorVdw[0]          = oneMinusLambdaVdw;
                lambdaFactorVdw[1]          = lambdaVdw;
                softcoreLambdaFactorCoul[0] = lambdaCoul;
                softcoreLambdaFactorCoul[1] = oneMinusLambdaCoul;
                softcoreLambdaFactorVdw[0]  = lambdaVdw;
                softcoreLambdaFactorVdw[1]  = oneMinusLambdaVdw;

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

                scalarForcePerDistanceVdw[0] = scalarForcePerDistanceVdw[1] = 0.0F;
                scalarForcePerDistanceCoul[0] = scalarForcePerDistanceCoul[1] = 0.0F;

                if (pairIncluded && withinCutoffMask)
                {
                    if (useSoftCore)
                    {
                        rpm2 = r2 * r2;
                        rp   = rpm2 * r2;
                    }
                    else
                    {
                        rpm2 = inv_r * inv_r;
                        rp   = 1.0F;
                    }

                    for (int k = 0; k < 2; k++)
                    {
#    ifndef LJ_COMB
                        /* LJ 6*C6 and 12*C12 */
                        fetch_nbfp_c6_c12(c6AB[k], c12AB[k], nbparam, numTypes * typeiAB[k] + typejAB[k]);
                        if (useSoftCore)
                        {
                            convert_c6_c12_to_sigma6(
                                    c6AB[k], c12AB[k], &(sigma6[k]), sigma6Minimum, sigma6WithInvalidSigma);
                        }
#    else
#        ifdef LJ_COMB_GEOM
                        c6AB[k]  = ljcp_iAB[k].x * ljcp_jAB[k].x;
                        c12AB[k] = ljcp_iAB[k].y * ljcp_jAB[k].y;
                        if (useSoftCore)
                        {
                            convert_c6_c12_to_sigma6(
                                    c6AB[k], c12AB[k], &(sigma6[k]), sigma6Minimum, sigma6WithInvalidSigma);
                        }
#        else
                        /* LJ 2^(1/6)*sigma and 12*epsilon */
                        sigmaAB[k] = ljcp_iAB[k].x + ljcp_jAB[k].x;
                        if (ljcp_iAB[k].x == 0.0F || ljcp_jAB[k].x == 0.0F)
                        {
                            sigmaAB[k] = 0.0F;
                        }
                        epsilonAB[k] = ljcp_iAB[k].y * ljcp_jAB[k].y;
                        convert_sigma_epsilon_to_c6_c12(sigmaAB[k], epsilonAB[k], &(c6AB[k]), &(c12AB[k]));
                        if (useSoftCore)
                        {
                            if ((c6AB[k] > 0.0F) && (c12AB[k] > 0.0F))
                            {
                                float sigma2 = sigmaAB[k] * sigmaAB[k];
                                sigma6[k]    = min(sigma2 * sigma2 * sigma2 * 0.5F, sigma6Minimum);
                            }
                            else
                            {
                                sigma6[k] = sigma6WithInvalidSigma;
                            }
                        }

#        endif /* LJ_COMB_GEOM */
#    endif     /* LJ_COMB */
                    }
                    if (useSoftCore)
                    {
                        if ((c12AB[0] > 0.0F) && (c12AB[1] > 0.0F))
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

                    for (int k = 0; k < 2; k++)
                    {
                        scalarForcePerDistanceVdw[k]  = 0.0F;
                        scalarForcePerDistanceCoul[k] = 0.0F;
                        Vcoul[k]                      = 0.0F;
                        Vvdw[k]                       = 0.0F;

                        bool nonZeroState = ((qq[k] != 0.0F) || (c6AB[k] != 0.0F) || (c12AB[k] != 0.0F));
                        if (nonZeroState)
                        {
                            if (useSoftCore)
                            {
                                rPInvC = 1.0F
                                         / (alphaCoulombEff * softcoreLambdaFactorCoul[k] * sigma6[k] + rp);
                                r2C   = rcbrt(rPInvC);
                                rInvC = rsqrt(r2C);

                                // equivalent to scLambdasOrAlphasDiffer
                                if ((alphaCoulombEff != alphaVdwEff)
                                    || (softcoreLambdaFactorVdw[k] != softcoreLambdaFactorCoul[k]))
                                {
                                    rPInvV = 1.0F
                                             / (alphaVdwEff * softcoreLambdaFactorVdw[k] * sigma6[k] + rp);
#    if defined LJ_POT_SWITCH || defined LJ_FORCE_SWITCH
                                    r2V   = rcbrt(rPInvV);
                                    rInvV = rsqrt(r2V);
#    endif
                                }
                                else
                                {
                                    /* We can avoid one expensive pow and one / operation */
                                    rPInvV = rPInvC;
#    if defined LJ_POT_SWITCH || defined LJ_FORCE_SWITCH
                                    r2V   = r2C;
                                    rInvV = rInvC;
#    endif
                                }
                            }
                            else
                            {
                                rPInvC = 1.0F;
                                r2C    = r2;
                                rInvC  = inv_r;
                                rPInvV = 1.0F;
#    if defined LJ_POT_SWITCH || defined LJ_FORCE_SWITCH
                                r2V   = r2;
                                rInvV = inv_r;
#    endif
                            }

                            if (c6AB[k] != 0.0F || c12AB[k] != 0.0F)
                            {
                                if (!useSoftCore)
                                {
                                    rInv6 = inv_r2 * inv_r2 * inv_r2;
                                }
                                else
                                {
                                    rInv6 = rPInvV;
                                }
                                float Vvdw6                  = c6AB[k] * rInv6;
                                float Vvdw12                 = c12AB[k] * rInv6 * rInv6;
                                scalarForcePerDistanceVdw[k] = Vvdw12 - Vvdw6;

                                Vvdw[k] = ((Vvdw12 + c12AB[k] * nbparam.repulsion_shift.cpot) * c_oneTwelfth
                                           - (Vvdw6 + c6AB[k] * nbparam.dispersion_shift.cpot) * c_oneSixth);

#    ifdef LJ_FORCE_SWITCH
                                calculate_force_switch_Fr_E(nbparam,
                                                            c6AB[k],
                                                            c12AB[k],
                                                            rInvV,
                                                            r2V,
                                                            &(scalarForcePerDistanceVdw[k]),
                                                            &(Vvdw[k]));
#    endif /* LJ_POT_SWITCH */

#    ifdef LJ_POT_SWITCH
                                calculate_potential_switch_Fr_E(
                                        nbparam, rInvV, r2V, &(scalarForcePerDistanceVdw[k]), &(Vvdw[k]));
#    endif /* LJ_POT_SWITCH */

#    ifdef VDW_CUTOFF_CHECK
                                /* Separate VDW cut-off check to enable twin-range cut-offs
                                 * (rvdw < rcoulomb <= rlist)
                                 */
                                vdw_in_range = (r2 < rCutoffVdwSq) ? 1.0F : 0.0F;
                                Vvdw[k] *= vdw_in_range;
#    endif /* VDW_CUTOFF_CHECK */
                            }

                            if (qq[k] != 0.0F)
                            {
#    ifdef EL_CUTOFF
#        ifdef EXCLUSION_FORCES
                                scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
#        else
                                scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
#        endif
#    endif
#    ifdef EL_RF
                                scalarForcePerDistanceCoul[k] = qq[k] * (rInvC - two_k_rf * r2C);
#    endif
#    if defined EL_EWALD_ANY
                                scalarForcePerDistanceCoul[k] = qq[k] * rInvC;
#    endif /* EL_EWALD_ANA/TAB */
#    ifdef EL_CUTOFF
                                Vcoul[k] = qq[k] * (rInvC - c_rf);
#    endif
#    ifdef EL_RF
                                Vcoul[k] = qq[k] * (rInvC + 0.5F * two_k_rf * r2C - c_rf);
#    endif
#    ifdef EL_EWALD_ANY
                                /* 1.0f - erff is faster than erfcf */
                                Vcoul[k] = qq[k] * (rInvC - ewald_shift);
#    endif /* EL_EWALD_ANY */
                            }
                            scalarForcePerDistanceCoul[k] *= rPInvC;
                            scalarForcePerDistanceVdw[k] *= rPInvV;
                        }
                    } // end for (int k = 0; k < 2; k++)

                    for (int k = 0; k < 2; k++)
                    {
                        E_el += lambdaFactorCoul[k] * Vcoul[k];
                        E_lj += lambdaFactorVdw[k] * Vvdw[k];

                        if (useSoftCore)
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
                } // if (pairIncluded && withinCutoffMask)

// ELEC REACTIONFIELD part
#    if defined EL_CUTOFF || defined EL_RF
                if (!pairIncluded && j < nj1)
                {
#        if defined EL_CUTOFF
                    float VV = -nbparam.c_rf;
#        else
                    float VV = 0.5F * two_k_rf * r2 - nbparam.c_rf;
#        endif

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
#    endif

#    ifdef EL_EWALD_ANY
                if ((!pairIncluded || r2 < rCutoffCoulSq) && j < nj1)
                {
                    v_lr = inv_r > 0.0F ? inv_r * erff(r2 * inv_r * beta) : 2.0F * beta * M_FLOAT_1_SQRTPI;
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
#    endif

                reduce_energy_warp_shfl(
                        E_lj, E_el, gm_e_lj + lambdaIdx, gm_e_el + lambdaIdx, tid, c_fullWarpMask);
                reduce_energy_warp_shfl(
                        DVDL_lj, DVDL_el, gm_dvdl_lj + lambdaIdx, gm_dvdl_el + lambdaIdx, tid, c_fullWarpMask);
            } // end for lambdaIdx in (0, nLambda)
            //} // end if (j < nj1)
        } // end for (int i = nj0; i < nj1; i += warp_size)
    } // end if (wid_global < numiAtoms)
}
#endif /* FUNCTION_DECLARATION_ONLY */

#undef EL_EWALD_ANY
#undef EXCLUSION_FORCES
#undef LJ_EWALD

#undef LJ_COMB
} // namespace gmx
