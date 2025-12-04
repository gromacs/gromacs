/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  NBNXM HIP kernels
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_HIP_NBNXM_HIP_KERNEL_BODY_H
#define GMX_NBNXM_HIP_NBNXM_HIP_KERNEL_BODY_H

#include <memory>
#include <type_traits>
#include <utility>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_hip_kernel.h"
#include "nbnxm_hip_kernel_utils.h"
#include "nbnxm_hip_types.h"

namespace gmx
{

//! Lookup kernel name based on launch configuration
template<ElecType elecType, VdwType vdwType, bool doPrune, bool doCalcEnergies, PairlistType pairlistType>
static const std::string getKernelName()
{

    // using std::array here instead of enumeration array because enumeration array [] is not constexpr
    static constexpr std::array<const std::string_view, c_numElecTypes> elecNames = {
        "_el-cut",           "_el-rf", "_el-ewald-tab", "_el-ewald-tab-twin", "_el_ewald-ana",
        "_el-ewald-ana-twin"
    };
    // using std::array here instead of enumeration array because enumeration array [] is not constexpr
    static constexpr std::array<const std::string_view, c_numVdwTypes> vdwNames = {
        "_lj-cut",        "_lj-comb-geom",       "_lj-comb-lb",      "_lj-force-switch",
        "_lj-pot-switch", "_lj-ewald-comb-geom", "_lj-ewald-comb-lb"
    };
    static constexpr std::string_view baseName = "nbnxm";
    static constexpr std::string_view elecName = elecNames[static_cast<int>(elecType)];
    static constexpr std::string_view vdwName  = vdwNames[static_cast<int>(vdwType)];

    static constexpr std::string_view pruneName            = "_prune";
    static constexpr std::string_view noPruneName          = "_noprune";
    static constexpr std::string_view calcEnergyName       = "_VF";
    static constexpr std::string_view noCalcEnergyName     = "_F";
    static constexpr std::string_view executionWidth64Name = "_wave64";
    static constexpr std::string_view executionWidth32Name = "_wave32";

    constexpr bool isWave64 = sc_gpuParallelExecutionWidth(pairlistType) == 64;

#if !defined(_MSC_VER)
    return std::string(CompileTimeStringJoin_v<baseName,
                                               elecName,
                                               vdwName,
                                               (doPrune ? pruneName : noPruneName),
                                               (doCalcEnergies ? calcEnergyName : noCalcEnergyName),
                                               (isWave64 ? executionWidth64Name : executionWidth32Name)>);
#else
    std::string returnValue;
    returnValue.reserve(1024);
    return returnValue.append(baseName)
            .append(elecName)
            .append(vdwName)
            .append(doPrune ? pruneName : noPruneName)
            .append(doCalcEnergies ? calcEnergyName : noCalcEnergyName)
            .append(isWave64 ? executionWidth64Name : executionWidth32Name);
#endif
}


//! Helper method to calculate launch bounds
static constexpr int minBlocksPerMp(const bool hasLargeRegisterPool, const bool doCalcEnergies)
{
    return hasLargeRegisterPool ? 1 : doCalcEnergies ? 6 : 8;
}


/*! \brief Main kernel for NBNXM.
 *
 */
template<bool hasLargeRegisterPool, bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType, int nthreadZ, PairlistType pairlistType>
__launch_bounds__(c_clSizeSq<pairlistType>* nthreadZ,
                  minBlocksPerMp(hasLargeRegisterPool, doCalcEnergies)) __global__
        static void nbnxmKernel(NBAtomDataGpu atdat, NBParamGpu nbparam, GpuPairlist<pairlistType> plist, bool doCalcShift)
{
    static constexpr EnergyFunctionProperties<elecType, vdwType> props;

    constexpr int c_clusterPerSuperCluster = sc_gpuClusterPerSuperCluster(pairlistType);
    constexpr int c_gpuJGroupSize          = sc_gpuJgroupSize(pairlistType);
    constexpr int c_clSize                 = sc_gpuClusterSize(pairlistType);
    constexpr int c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);

    /* thread/block/warp id-s */
    const int tidxi = threadIdx.x;
    const int tidxj = threadIdx.y;
    const int tidx  = tidxj * c_clSize + tidxi;
    const int tidxz = nthreadZ == 1 ? 0 : threadIdx.z;
    const int widx =
            (c_clSizeSq<pairlistType> == c_parallelExecutionWidth) ? 0 : tidx / c_subWarp<pairlistType>;
    const int bidx       = blockIdx.x;
    const int tidxInWarp = tidx & (c_parallelExecutionWidth - 1);

    const unsigned int energyIndexBase =
            1 + (bidx & (sc_energyVirialNumElementsSeparateDeviceReduction - 1));

    using NbnxmExcl     = nbnxn_excl_t<pairlistType>;
    using NbnxmCjPacked = nbnxn_cj_packed_t<pairlistType>;

    AmdFastBuffer<const float4> gm_xq{ atdat.xq };
    float3*                     gm_f             = asFloat3(atdat.f);
    float3*                     gm_shiftVec      = asFloat3(atdat.shiftVec);
    float3*                     gm_fShift        = asFloat3(atdat.fShift);
    float*                      gm_energyElec    = atdat.eElec + energyIndexBase;
    float*                      gm_energyVdw     = atdat.eLJ + energyIndexBase;
    NbnxmCjPacked*              gm_plistCJPacked = plist.cjPacked;
    AmdFastBuffer<const nbnxn_sci_t> gm_plistSci{ doPruneNBL ? plist.sci : plist.sorting.sciSorted };
    int* gm_plistSciHistogram = plist.sorting.sciHistogram;
    int* gm_sciCount          = plist.sorting.sciCount;


    AmdFastBuffer<const NbnxmExcl> gm_plistExcl{ plist.excl };
    AmdFastBuffer<const Float2>    gm_ljComb{ atdat.ljComb };       /* used iff ljComb<vdwType> */
    AmdFastBuffer<const int>       gm_atomTypes{ atdat.atomTypes }; /* used iff !ljComb<vdwType> */
    const Float2*                  gm_nbfp     = nbparam.nbfp;      /* used iff !ljComb<vdwType> */
    const Float2*                  gm_nbfpComb = nbparam.nbfp_comb; /* used iff ljEwald<vdwType> */
    const float* gm_coulombTab = nbparam.coulomb_tab; /* used iff elecEwaldTab<elecType> */

    const int             numTypes        = atdat.numTypes;
    const float           rCoulombSq      = nbparam.rcoulomb_sq;
    const float           rVdwSq          = nbparam.rvdw_sq;
    const float           twoKRf          = nbparam.two_k_rf;
    const float           ewaldBeta       = nbparam.ewald_beta;
    const float           rlistOuterSq    = nbparam.rlistOuter_sq;
    const float           ewaldShift      = nbparam.sh_ewald;
    const float           epsFac          = nbparam.epsfac;
    const float           ewaldCoeffLJ_2  = nbparam.ewaldcoeff_lj * nbparam.ewaldcoeff_lj;
    const float           cRF             = nbparam.c_rf;
    const shift_consts_t  dispersionShift = nbparam.dispersion_shift;
    const shift_consts_t  repulsionShift  = nbparam.repulsion_shift;
    const switch_consts_t vdwSwitch       = nbparam.vdw_switch;
    const float           rVdwSwitch      = nbparam.rvdw_switch;
    const float           ljEwaldShift    = nbparam.sh_lj_ewald;
    const float           coulombTabScale = nbparam.coulomb_tab_scale;

    // We use the regular extern __shared__ declaration for the shared memory buffer instead of
    // the legacy HIP_DYNAMIC_SHARED macro, as our minimum required version already properly
    // supports the regular syntax in the amd llvm compiler
    extern __shared__ char sm_reductionBuffer[];
    char*                  sm_nextSlotPtr = sm_reductionBuffer;

    // We rely on the flat shared memory space here, and abuse undefined behavior of changing
    // the data type based on what we expect on the specific offset
    float4* sm_xqBufferPtr = reinterpret_cast<float4*>(sm_nextSlotPtr);
    sm_nextSlotPtr += incrementSharedMemorySlotPtr<pairlistType, float4>();
    auto sm_atomTypeI = [&]()
    {
        int* temp = nullptr;
        if constexpr (!props.vdwComb)
        {
            temp = reinterpret_cast<int*>(sm_nextSlotPtr);
            sm_nextSlotPtr += incrementSharedMemorySlotPtr<pairlistType, int>();
        }
        return temp;
    }();

    auto sm_ljCombI = [&]()
    {
        float2* temp = nullptr;
        if constexpr (props.vdwComb)
        {
            temp = reinterpret_cast<float2*>(sm_nextSlotPtr);
            sm_nextSlotPtr += incrementSharedMemorySlotPtr<pairlistType, float2>();
        }
        return temp;
    }();

    auto sm_prunedPairCount = reinterpret_cast<int*>(sm_nextSlotPtr);
    sm_nextSlotPtr += incrementSharedMemorySlotPtr<pairlistType, int>();


    AmdPackedFloat3 fCiBuffer[c_clusterPerSuperCluster]; // i force buffer
    for (int i = 0; i < c_clusterPerSuperCluster; i++)
    {
        fCiBuffer[i] = { 0.0F, 0.0F, 0.0F };
    }

    const nbnxn_sci_t nbSci          = gm_plistSci[bidx];
    const int         sci            = nbSci.sci;
    const int         cijPackedBegin = nbSci.cjPackedBegin;
    const int         cijPackedEnd   = nbSci.cjPackedEnd;

    /*! i-cluster interaction mask for a super-cluster with all c_clusterPerSuperCluster bits set */
    constexpr int perSuperClInteractionMask = ((1U << c_clusterPerSuperCluster) - 1U);

    // Only needed if props.elecEwaldAna
    const float beta2 = ewaldBeta * ewaldBeta;
    const float beta3 = ewaldBeta * ewaldBeta * ewaldBeta;

    // We may need only a subset of threads active for preloading i-atoms
    // depending on the super-cluster and cluster / thread-block size.
    constexpr bool c_loadUsingAllXYThreads = (c_clSize == c_clusterPerSuperCluster);
    if (tidxz == 0 && (c_loadUsingAllXYThreads || tidxj < c_clusterPerSuperCluster))
    {
        /* Pre-load i-atom x and q into shared memory */
        const int    ci       = sci * c_clusterPerSuperCluster + tidxj;
        const int    ai       = ci * c_clSize + tidxi;
        const float3 shift    = gm_shiftVec[nbSci.shift];
        const int    cacheIdx = tidxj * c_clSize + tidxi;
        float4       xqbuf    = gm_xq[ai];
        // TODO: Remove `-` and reverse operators in `xi + xj` and `+- f_ij` when it's fixed.
        // For some reason the compiler does not generate v_pk_add_f32 and v_sub_f32 for `xi -
        // xj` but generates 3 v_sub_f32. Hence all this mess with signs.
        xqbuf.x = -(xqbuf.x + shift.x);
        xqbuf.y = -(xqbuf.y + shift.y);
        xqbuf.z = -(xqbuf.z + shift.z);
        xqbuf.w *= epsFac;

        sm_xqBufferPtr[cacheIdx] = xqbuf;

        if constexpr (!props.vdwComb)
        {
            // Pre-load the i-atom types into shared memory
            sm_atomTypeI[cacheIdx] = gm_atomTypes[ai];
        }
        else
        {
            // Pre-load the LJ combination parameters into shared memory
            sm_ljCombI[cacheIdx] = gm_ljComb[ai];
        }
    }

    if constexpr (doPruneNBL)
    {
        if (tidx == 0)
        {
            sm_prunedPairCount[0] = 0;
        }
    }
    int prunedPairCount = 0;

    __syncthreads();

    float ewaldCoeffLJ_6_6; // Only needed if (props.vdwEwald)
    if constexpr (props.vdwEwald)
    {
        ewaldCoeffLJ_6_6 = ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * c_oneSixth;
    }

    float energyVdw, energyElec; // Only needed if (doCalcEnergies)
    if constexpr (doCalcEnergies)
    {
        energyVdw = energyElec = 0.0F;
    }
    if constexpr (doCalcEnergies && props.doExclusionForces)
    {
        if (nbSci.shift == c_centralShiftIndex
            && gm_plistCJPacked[cijPackedBegin].cj[0] == sci * c_clusterPerSuperCluster)
        {
            // we have the diagonal: add the charge and LJ self interaction energy term
            for (int i = 0; i < c_clusterPerSuperCluster; i++)
            {
                // TODO: Are there other options?
                if constexpr (props.elecEwald || props.elecRF || props.elecCutoff)
                {
                    const float qi = sm_xqBufferPtr[i * c_clSize + tidxi].w;
                    energyElec += qi * qi;
                }
                if constexpr (props.vdwEwald)
                {
                    energyVdw += LDG(reinterpret_cast<const float*>(
                            &gm_nbfp[gm_atomTypes[(sci * c_clusterPerSuperCluster + i) * c_clSize + tidxi]
                                     * (numTypes + 1)]));
                }
            }
            /* divide the self term(s) equally over the j-threads, then multiply with the coefficients. */
            if constexpr (props.vdwEwald)
            {
                energyVdw /= c_clSize * nthreadZ;
                energyVdw *= 0.5F * c_oneSixth * ewaldCoeffLJ_6_6; // c_OneTwelfth?
            }
            if constexpr (props.elecRF || props.elecCutoff)
            {
                // Correct for epsfac^2 due to adding qi^2 */
                energyElec /= epsFac * c_clSize * nthreadZ;
                energyElec *= -0.5F * cRF;
            }
            if constexpr (props.elecEwald)
            {
                // Correct for epsfac^2 due to adding qi^2 */
                energyElec /= epsFac * c_clSize * nthreadZ;
                energyElec *= -ewaldBeta * c_oneOverSqrtPi; /* last factor 1/sqrt(pi) */
            }
        } // (nbSci.shift == c_centralShiftIndex && a_plistCJPacked[cijPackedBegin].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
    } // (doCalcEnergies && doExclusionForces)

    // Only needed if (doExclusionForces)
    const bool nonSelfInteraction = !(nbSci.shift == c_centralShiftIndex & tidxj <= tidxi);

    // loop over the j clusters = seen by any of the atoms in the current super-cluster
    for (int jPacked = cijPackedBegin; jPacked < cijPackedEnd; ++jPacked)
    {
        int imask = gm_plistCJPacked[jPacked].imei[widx].imask;
        imask     = (c_clSizeSq<pairlistType> == c_parallelExecutionWidth)
                            ? __builtin_amdgcn_readfirstlane(imask)
                            : imask;
        if (!doPruneNBL && !imask)
        {
            continue;
        }
        const int wexclIdx = gm_plistCJPacked[jPacked].imei[widx].excl_ind;

        const int wexcl = gm_plistExcl[wexclIdx].pair[tidx & (c_subWarp<pairlistType> - 1)];
#pragma unroll c_gpuJGroupSize
        for (int jm = 0; jm < c_gpuJGroupSize; jm++)
        {
            const bool maskSet = imask & (perSuperClInteractionMask << (jm * c_clusterPerSuperCluster));
            if (!maskSet)
            {
                continue;
            }
            int       maskJI = (1U << (jm * c_clusterPerSuperCluster));
            const int cj     = gm_plistCJPacked[jPacked].cj[jm];
            const int aj     = cj * c_clSize + tidxj;

            // load j atom data
            const Float4 xqjbuf = gm_xq[aj];

            const AmdPackedFloat3 xj(xqjbuf.x, xqjbuf.y, xqjbuf.z);
            const float           qj = xqjbuf.w;
            int                   atomTypeJ; // Only needed if (!props.vdwComb)
            float2                ljCombJ;   // Only needed if (props.vdwComb)
            if constexpr (props.vdwComb)
            {
                ljCombJ = gm_ljComb[aj];
            }
            else
            {
                atomTypeJ = gm_atomTypes[aj];
            }

            AmdPackedFloat3 fCjBuf(0.0F, 0.0F, 0.0F);

#pragma unroll c_clusterPerSuperCluster
            for (int i = 0; i < c_clusterPerSuperCluster; i++)
            {
                if (imask & maskJI)
                {
                    // i cluster index
                    const int ci = sci * c_clusterPerSuperCluster + i;
                    // all threads load an atom from i cluster ci into shmem!
                    const Float4          xqibuf = sm_xqBufferPtr[i * c_clSize + tidxi];
                    const AmdPackedFloat3 xi(xqibuf.x, xqibuf.y, xqibuf.z);

                    // distance between i and j atoms
                    const AmdPackedFloat3 rv = xi + xj;
                    float                 r2 = rv.norm2();

                    if constexpr (doPruneNBL)
                    {
                        /* If _none_ of the atoms pairs are in cutoff range,
                         * the bit corresponding to the current
                         * cluster-pair in imask gets set to 0. */
                        if (!nb_any_internal<pairlistType>(r2 < rlistOuterSq, widx))
                        {
                            imask &= ~maskJI;
                        }
                    }
                    const float pairExclMask = (wexcl >> (jm * c_clusterPerSuperCluster + i)) & 1;

                    // cutoff & exclusion check

                    const bool notExcluded = props.doExclusionForces
                                                     ? (ci != (nonSelfInteraction ? -1 : cj))
                                                     : (1 * pairExclMask);

                    if ((r2 < rCoulombSq) && notExcluded)
                    {
                        const float qi = xqibuf.w;
                        int         atomTypeI; // Only needed if (!props.vdwComb)
                        float       sigma, epsilon;
                        Float2      c6c12;

                        if constexpr (!props.vdwComb)
                        {
                            /* LJ 6*C6 and 12*C12 */
                            atomTypeI = sm_atomTypeI[i * c_clSize + tidxi];
                            if constexpr (sc_gpuParallelExecutionWidth(pairlistType) == 64)
                            {
                                c6c12 = fetchNbfpC6C12(gm_nbfp, __mul24(numTypes, atomTypeI) + atomTypeJ);
                            }
                            else
                            {
                                c6c12 = fetchNbfpC6C12(gm_nbfp, numTypes * atomTypeI + atomTypeJ);
                            }
                        }
                        else
                        {
                            const float2 ljCombI = sm_ljCombI[i * c_clSize + tidxi];
                            if constexpr (props.vdwCombGeom)
                            {
                                c6c12 = ljCombI * ljCombJ;
                            }
                            else
                            {
                                static_assert(props.vdwCombLB);
                                // LJ 2^(1/6)*sigma and 12*epsilon
                                sigma   = ljCombI.x + ljCombJ.x;
                                epsilon = ljCombI.y * ljCombJ.y;
                                if constexpr (doCalcEnergies)
                                {
                                    c6c12 = convertSigmaEpsilonToC6C12(sigma, epsilon);
                                }
                            } // props.vdwCombGeom
                        } // !props.vdwComb

                        // c6 and c12 are unused and garbage iff props.vdwCombLB && !doCalcEnergies
                        const float c6  = c6c12.x;
                        const float c12 = c6c12.y;

                        // Ensure distance do not become so small that r^-12 overflows
                        r2                = fmax(r2, c_nbnxnMinDistanceSquared);
                        const float rInv  = __frsqrt_rn(r2);
                        const float r2Inv = rInv * rInv;
                        float       r6Inv, fInvR, energyLJPair;
                        if constexpr (!props.vdwCombLB || doCalcEnergies)
                        {
                            r6Inv = r2Inv * r2Inv * r2Inv;
                            if constexpr (props.doExclusionForces)
                            {
                                /* We could mask r2Inv, but with Ewald masking both
                                 * r6Inv and fInvR is faster */
                                r6Inv *= pairExclMask;
                            }
                            fInvR = r6Inv * (c12 * r6Inv - c6) * r2Inv;
                        }
                        else
                        {
                            float sig_r  = sigma * rInv;
                            float sig_r2 = sig_r * sig_r;
                            float sig_r6 = sig_r2 * sig_r2 * sig_r2;
                            if constexpr (props.doExclusionForces)
                            {
                                sig_r6 *= pairExclMask;
                            }
                            fInvR = epsilon * sig_r6 * (sig_r6 - 1.0F) * r2Inv;
                        } // (!props.vdwCombLB || doCalcEnergies)
                        if constexpr (doCalcEnergies || props.vdwPSwitch)
                        {
                            energyLJPair = pairExclMask
                                           * (c12 * (r6Inv * r6Inv + repulsionShift.cpot) * c_oneTwelfth
                                              - c6 * (r6Inv + dispersionShift.cpot) * c_oneSixth);
                        }
                        if constexpr (props.vdwFSwitch)
                        {
                            ljForceSwitch<doCalcEnergies>(
                                    dispersionShift, repulsionShift, c6c12, rVdwSwitch, rInv, r2, &fInvR, &energyLJPair);
                        }
                        if constexpr (props.vdwEwald)
                        {
                            ljEwaldComb<doCalcEnergies, vdwType>(gm_nbfpComb,
                                                                 ljEwaldShift,
                                                                 atomTypeI,
                                                                 atomTypeJ,
                                                                 r2,
                                                                 r2Inv,
                                                                 ewaldCoeffLJ_2,
                                                                 ewaldCoeffLJ_6_6,
                                                                 pairExclMask,
                                                                 &fInvR,
                                                                 &energyLJPair);
                        } // (props.vdwEwald)
                        if constexpr (props.vdwPSwitch)
                        {
                            ljPotentialSwitch<doCalcEnergies>(
                                    vdwSwitch, rVdwSwitch, rInv, r2, &fInvR, &energyLJPair);
                        }
                        if constexpr (props.elecEwaldTwin)
                        {
                            // Separate VDW cut-off check to enable twin-range cut-offs
                            // (rVdw < rCoulomb <= rList)
                            const float vdwInRange = (r2 < rVdwSq) ? 1.0F : 0.0F;
                            fInvR *= vdwInRange;
                            if constexpr (doCalcEnergies)
                            {
                                energyLJPair *= vdwInRange;
                            }
                        }
                        if constexpr (doCalcEnergies)
                        {
                            energyVdw += energyLJPair;
                        }

                        if constexpr (props.elecCutoff)
                        {
                            if constexpr (props.doExclusionForces)
                            {
                                fInvR += qi * qj * pairExclMask * r2Inv * rInv;
                            }
                            else
                            {
                                fInvR += qi * qj * r2Inv * rInv;
                            }
                        }
                        if constexpr (props.elecRF)
                        {
                            fInvR += qi * qj * (pairExclMask * r2Inv * rInv - twoKRf);
                        }
                        if constexpr (props.elecEwaldAna)
                        {
                            fInvR += qi * qj * (pairExclMask * r2Inv * rInv + pmeCorrF(beta2 * r2) * beta3);
                        }
                        if constexpr (props.elecEwaldTab)
                        {
                            fInvR += qi * qj
                                     * (pairExclMask * r2Inv
                                        - interpolateCoulombForceR(gm_coulombTab, coulombTabScale, r2 * rInv))
                                     * rInv;
                        }

                        if constexpr (doCalcEnergies)
                        {
                            if constexpr (props.elecCutoff)
                            {
                                energyElec += qi * qj * (pairExclMask * rInv - cRF);
                            }
                            if constexpr (props.elecRF)
                            {
                                energyElec += qi * qj * (pairExclMask * rInv + 0.5F * twoKRf * r2 - cRF);
                            }
                            if constexpr (props.elecEwald)
                            {
                                energyElec += qi * qj
                                              * (rInv * (pairExclMask - erff(r2 * rInv * ewaldBeta))
                                                 - pairExclMask * ewaldShift);
                            }
                        }

                        const AmdPackedFloat3 forceIJ = rv * fInvR;

                        /* accumulate j forces in registers */
                        fCjBuf += forceIJ;
                        /* accumulate i forces in registers */
                        fCiBuffer[i] -= forceIJ;
                    } // (r2 < rCoulombSq) && notExcluded
                } // (imask & maskJI)
                /* shift the mask bit by 1 */
                maskJI += maskJI;
            } // for (int i = 0; i < c_clusterPerSuperCluster; i++)
            /* reduce j forces */
            const float reducedForceJ = reduceForceJWarpShuffle<pairlistType>(fCjBuf, tidxi);
            if (tidxi < 3)
            {
                amdFastAtomicAddForce(gm_f, aj, tidxi, reducedForceJ);
            }
        } // for (int jm = 0; jm < c_gpuJGroupSize; jm++)
        if constexpr (doPruneNBL)
        {
            /* Update the imask with the new one which does not contain the
             * out of range clusters anymore. */
            gm_plistCJPacked[jPacked].imei[widx].imask = imask;
            prunedPairCount += __popc(imask);
        }
    } // for (int jPacked = cijPackedBegin; jPacked < cijPackedEnd; jPacked += 1)

    if (doCalcShift && nbSci.shift != c_centralShiftIndex)
    {
        reduceForceI<true, pairlistType>(fCiBuffer, gm_f, tidxi, tidxj, tidx, sci, gm_fShift, nbSci.shift);
    }
    else
    {
        reduceForceI<false, pairlistType>(fCiBuffer, gm_f, tidxi, tidxj, tidx, sci, gm_fShift, nbSci.shift);
    }

    if constexpr (doCalcEnergies)
    {
        reduceEnergyWarpShuffle<pairlistType>(energyVdw, energyElec, gm_energyVdw, gm_energyElec, tidx);
    }

    if constexpr (doPruneNBL)
    {
        if (tidxInWarp == 0)
        {
            atomicAdd(sm_prunedPairCount, prunedPairCount);
            __syncthreads();
            prunedPairCount = *sm_prunedPairCount;

            if (tidxi == 0 && tidxj == 0 && tidxz == 0)
            {
                int index = max(c_sciHistogramSize - prunedPairCount - 1, 0);
                atomicAdd(gm_plistSciHistogram + index, 1);
                gm_sciCount[bidx] = index;
            }
        }
    }
}

//! \brief NBNXM kernel launch code.
template<PairlistType pairlistType, bool hasLargeRegisterPool, bool doPrune, bool doCalcEnergies, ElecType elecType, VdwType vdwType, class... Args>
static void launchNbnxmKernel(const DeviceStream&      deviceStream,
                              const int                numSci,
                              const DeviceInformation& deviceInfo,
                              Args*... args)
{
    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */

    constexpr size_t   numThreadZ    = 1;
    constexpr int      c_clSize      = sc_gpuClusterSize(pairlistType);
    constexpr bool     isPruneKernel = false;
    KernelLaunchConfig config;
    config.gridSize[0]  = numberOfKernelBlocksSanityCheck(numSci, deviceInfo);
    config.blockSize[0] = c_clSize;
    config.blockSize[1] = c_clSize;
    config.blockSize[2] = numThreadZ;
    config.sharedMemorySize =
            requiredSharedMemorySize<isPruneKernel, numThreadZ, vdwType, pairlistType>();

    auto kernel =
            nbnxmKernel<hasLargeRegisterPool, doPrune, doCalcEnergies, elecType, vdwType, numThreadZ, pairlistType>;
    auto kernelName = getKernelName<elecType, vdwType, doPrune, doCalcEnergies, pairlistType>();

    const auto kernelArgs = prepareGpuKernelArguments(kernel, config, args...);


    launchGpuKernel(kernel, config, deviceStream, nullptr, kernelName.c_str(), kernelArgs);
}

//! \brief Select templated kernel and launch it.
template<PairlistType pairlistType, bool hasLargeRegisterPool, bool doPruneNBL, bool doCalcEnergies, class... Args>
void chooseAndLaunchNbnxmKernel(ElecType                 elecType,
                                VdwType                  vdwType,
                                const DeviceStream&      deviceStream,
                                const int                numSci,
                                const DeviceInformation& deviceInfo,
                                Args*... args)
{
    dispatchTemplatedFunction(
            [&](auto elecType_, auto vdwType_)
            {
                return launchNbnxmKernel<pairlistType, hasLargeRegisterPool, doPruneNBL, doCalcEnergies, elecType_, vdwType_>(
                        deviceStream, numSci, deviceInfo, args...);
            },
            elecType,
            vdwType);
}

template<bool hasLargeRegisterPool, bool doPruneNBL, bool doCalcEnergies>
void launchNbnxmKernelHelper(NbnxmGpu* nb, const StepWorkload& stepWork, const InteractionLocality iloc)
{
    NBAtomDataGpu*           adat         = nb->atdat;
    NBParamGpu*              nbp          = nb->nbparam;
    const DeviceStream&      deviceStream = *nb->deviceStreams[iloc];
    const DeviceInformation& deviceInfo   = nb->deviceContext_->deviceInfo();

    std::visit(
            [&](auto&& pairlists)
            {
                auto* plist                 = pairlists[iloc].get();
                using T                     = std::decay_t<decltype(*plist)>;
                constexpr auto pairlistType = getPairlistTypeFromPairlist<T>();

                GMX_ASSERT(doPruneNBL == (plist->haveFreshList && !nb->didPrune[iloc]),
                           "Wrong template called");
                GMX_ASSERT(doCalcEnergies == stepWork.computeEnergy, "Wrong template called");

                chooseAndLaunchNbnxmKernel<pairlistType, hasLargeRegisterPool, doPruneNBL, doCalcEnergies>(
                        nbp->elecType,
                        nbp->vdwType,
                        deviceStream,
                        plist->numSci,
                        deviceInfo,
                        adat,
                        nbp,
                        plist,
                        &stepWork.computeVirial);
            },
            nb->plist);
}

} // namespace gmx


#endif
