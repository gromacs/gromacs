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

/*! \internal \file
 *  \brief
 *  CUDA non-bonded kernel used through preprocessor-based code generation
 *  of multiple kernel flavors, see nbnxn_cuda_kernels.cuh.
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Berk Hess <hess@kth.se>
 *  \ingroup module_nbnxm
 */

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

#if defined LJ_EWALD_COMB_GEOM || defined LJ_EWALD_COMB_LB
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#    define LJ_EWALD
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

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
#    define LJ_COMB
#endif

/*
   Kernel launch parameters:
    - #blocks   = #pair lists, blockId = pair list Id
    - #threads  = NTHREAD_Z * c_clSize^2
    - shmem     = see nbnxn_cuda.cu:calc_shmem_required_nonbonded()

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */

/**@{*/
/*! \brief Compute capability dependent definition of kernel launch configuration parameters.
 *
 * NTHREAD_Z controls the number of j-clusters processed concurrently on NTHREAD_Z
 * warp-pairs per block.
 *
 * - On CC 3.0-3.5, and >=5.0 NTHREAD_Z == 1, translating to 64 th/block with 16
 * blocks/multiproc, is the fastest even though this setup gives low occupancy
 * (except on 6.0).
 * NTHREAD_Z > 1 results in excessive register spilling unless the minimum blocks
 * per multiprocessor is reduced proportionally to get the original number of max
 * threads in flight (and slightly lower performance).
 * - On CC 3.7 there are enough registers to double the number of threads; using
 * NTHREADS_Z == 2 is fastest with 16 blocks (TODO: test with RF and other kernels
 * with low-register use).
 *
 * Note that the current kernel implementation only supports NTHREAD_Z > 1 with
 * shuffle-based reduction, hence CC >= 3.0.
 *
 *
 * NOTEs on Volta / CUDA 9 extensions:
 *
 * - While active thread masks are required for the warp collectives
 *   (we use any and shfl), the kernel is designed such that all conditions
 *   (other than the inner-most distance check) including loop trip counts
 *   are warp-synchronous. Therefore, we don't need ballot to compute the
 *   active masks as these are all full-warp masks.
 *
 */

/* Kernel launch bounds for different compute capabilities. The value of NTHREAD_Z
 * determines the number of threads per block and it is chosen such that
 * 16 blocks/multiprocessor can be kept in flight.
 * - CC 3.0,3.5, and >=5.0: NTHREAD_Z=1, (64, 16) bounds
 * - CC 3.7:                NTHREAD_Z=2, (128, 16) bounds
 *
 * Note: convenience macros, need to be undef-ed at the end of the file.
 */
#if GMX_PTX_ARCH == 370
#    define NTHREAD_Z (2)
#    define MIN_BLOCKS_PER_MP (16)
#else
#    define NTHREAD_Z (1)
#    define MIN_BLOCKS_PER_MP (16)
#endif /* GMX_PTX_ARCH == 370 */
#define THREADS_PER_BLOCK (c_clSize * c_clSize * NTHREAD_Z)

#if GMX_PTX_ARCH >= 350
/**@}*/
__launch_bounds__(THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)
#else
__launch_bounds__(THREADS_PER_BLOCK)
#endif /* GMX_PTX_ARCH >= 350 */
#ifdef PRUNE_NBL
#    ifdef CALC_ENERGIES
        __global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_prune_cuda)
#    else
        __global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_prune_cuda)
#    endif /* CALC_ENERGIES */
#else
#    ifdef CALC_ENERGIES
        __global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_cuda)
#    else
        __global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_cuda)
#    endif /* CALC_ENERGIES */
#endif     /* PRUNE_NBL */
                (NBAtomDataGpu atdat, NBParamGpu nbparam, GpuPairlist plist, bool bCalcFshift)
#ifdef FUNCTION_DECLARATION_ONLY
                        ; /* Only do function declaration, omit the function body. */
#else
{
    /* convenience variables */
#    ifdef PRUNE_NBL
    /* we can't use the sorted plist in this call as we need to use this kernel to perform counts
     * which will be used in the sorting */
    const nbnxn_sci_t*  pl_sci      = plist.sci;
#    else
    /* the sorted list has been generated using data from a previous call to this kernel */
    const nbnxn_sci_t* pl_sci = plist.sorting.sciSorted;
    const
#    endif
    nbnxn_cj_packed_t*  pl_cjPacked = plist.cjPacked;
    const nbnxn_excl_t* excl        = plist.excl;
#    ifndef LJ_COMB
    const int*          atom_types  = atdat.atomTypes;
    int                 ntypes      = atdat.numTypes;
#    else
    const float2* lj_comb = atdat.ljComb;
    float2        ljcp_i, ljcp_j;
#    endif
    const float4*       xq          = atdat.xq;
    float3*             f           = asFloat3(atdat.f);
    const float3*       shift_vec   = asFloat3(atdat.shiftVec);
    float               rcoulomb_sq = nbparam.rcoulomb_sq;
#    ifdef VDW_CUTOFF_CHECK
    float               rvdw_sq     = nbparam.rvdw_sq;
    float               vdw_in_range;
#    endif
#    ifdef LJ_EWALD
    float               lje_coeff2, lje_coeff6_6;
#    endif
#    ifdef EL_RF
    float               two_k_rf    = nbparam.two_k_rf;
#    endif
#    ifdef EL_EWALD_ANA
    float               beta2       = nbparam.ewald_beta * nbparam.ewald_beta;
    float               beta3       = nbparam.ewald_beta * nbparam.ewald_beta * nbparam.ewald_beta;
#    endif
#    ifdef PRUNE_NBL
    float               rlist_sq    = nbparam.rlistOuter_sq;
#    endif

#    ifdef CALC_ENERGIES
#        ifdef EL_EWALD_ANY
    float               beta        = nbparam.ewald_beta;
    float               ewald_shift = nbparam.sh_ewald;
#        else
    float                reactionFieldShift = nbparam.c_rf;
#        endif /* EL_EWALD_ANY */
    float*              e_lj        = atdat.eLJ;
    float*              e_el        = atdat.eElec;
#    endif     /* CALC_ENERGIES */

    /* thread/block/warp id-s */
    unsigned int tidxi      = threadIdx.x;
    unsigned int tidxj      = threadIdx.y;
    unsigned int tidx       = threadIdx.y * blockDim.x + threadIdx.x;
#    if NTHREAD_Z == 1
    unsigned int tidxz      = 0;
#    else
    unsigned int  tidxz = threadIdx.z;
#    endif
    unsigned int bidx       = blockIdx.x;
    unsigned int widx       = tidx / warp_size; /* warp index */
#    ifdef PRUNE_NBL
    unsigned int tidxInWarp = tidx & (warp_size - 1);
#    endif
    int          sci, ci, cj, ai, aj, cijPackedBegin, cijPackedEnd;
#    ifndef LJ_COMB
    int          typei, typej;
#    endif
    int          i, jm, jPacked, wexcl_idx;
    float        qi, qj_f, r2, inv_r, inv_r2;
#    if !defined LJ_COMB_LB || defined CALC_ENERGIES
    float        inv_r6, c6, c12;
#    endif
#    ifdef LJ_COMB_LB
    float        sigma, epsilon;
#    endif
    float        int_bit, F_invr;
#    ifdef CALC_ENERGIES
    float        E_lj, E_el;
#    endif
#    if defined CALC_ENERGIES || defined LJ_POT_SWITCH
    float        E_lj_p;
#    endif
    unsigned int wexcl, imask, mask_ji;
    float4       xqbuf;
    float3       xi, xj, rv, f_ij, fcj_buf;
    float3       fci_buf[c_nbnxnGpuNumClusterPerSupercluster]; /* i force buffer */
    nbnxn_sci_t  nb_sci;

    /*! i-cluster interaction mask for a super-cluster with all c_nbnxnGpuNumClusterPerSupercluster=8 bits set */
    const unsigned superClInteractionMask = ((1U << c_nbnxnGpuNumClusterPerSupercluster) - 1U);

    // cj preload is off in the following cases:
    // - sm_70 (V100), sm_80 (A100), sm_86 (GA02)
    // - for future arch (> 8.6 at the time of writing) we assume it is better to keep it off
    // cj preload is left on for:
    // - sm_75: improvements +/- very small
    // - sm_61: tested and slower without preload
    // - sm_6x and earlier not tested to
    constexpr bool c_preloadCj = (GMX_PTX_ARCH < 700 || GMX_PTX_ARCH == 750);


    // Full or partial unroll on Ampere (and later) GPUs is beneficial given the increased L1
    // instruction cache. Tested with CUDA 11-12.
#    if GMX_PTX_ARCH >= 800
#        define DO_JM_UNROLL 1
#        if !defined CALC_ENERGIES && !defined PRUNE_NBL
#            if (defined EL_CUTOFF || defined EL_RF                                            \
                 || defined EL_EWALD_ANY && !defined LJ_FORCE_SWITCH && !defined LJ_POT_SWITCH \
                            && (defined LJ_COMB_GEOM || GMX_PTX_ARCH == 800))
    static constexpr int jmLoopUnrollFactor = 4;
#            else
    static constexpr int jmLoopUnrollFactor = 2;
#            endif
#        else // CALC_ENERGIES
#            if (defined EL_CUTOFF || defined EL_RF && !defined LJ_FORCE_SWITCH && !defined LJ_POT_SWITCH)
    static constexpr int jmLoopUnrollFactor = 2;
#            else
    static constexpr int jmLoopUnrollFactor = 1;
#            endif
#        endif
#    else
#        define DO_JM_UNROLL 0
#    endif

    /*********************************************************************
     * Set up shared memory pointers.
     * sm_nextSlotPtr should always be updated to point to the "next slot",
     * that is past the last point where data has been stored.
     */
    // NOLINTNEXTLINE(readability-redundant-declaration)
    extern __shared__ char sm_dynamicShmem[];
    char*                  sm_nextSlotPtr = sm_dynamicShmem;
    static_assert(sizeof(char) == 1,
                  "The shared memory offset calculation assumes that char is 1 byte");

    /* shmem buffer for i x+q pre-loading */
    float4* xqib = reinterpret_cast<float4*>(sm_nextSlotPtr);
    sm_nextSlotPtr += (c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(*xqib));

    /* shmem buffer for cj, for each warp separately */
    int* cjs = reinterpret_cast<int*>(sm_nextSlotPtr);
    if (c_preloadCj)
    {
        /* the cjs buffer's use expects a base pointer offset for pairs of warps in the j-concurrent execution */
        cjs += tidxz * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize;
        sm_nextSlotPtr += (NTHREAD_Z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(*cjs));
    }

#    ifndef LJ_COMB
    /* shmem buffer for i atom-type pre-loading */
    int* atib = reinterpret_cast<int*>(sm_nextSlotPtr);
    sm_nextSlotPtr += (c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(*atib));
#    else
    /* shmem buffer for i-atom LJ combination rule parameters */
    float2* ljcpib = reinterpret_cast<float2*>(sm_nextSlotPtr);
    sm_nextSlotPtr += (c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(*ljcpib));
#    endif
    /*********************************************************************/

    nb_sci         = pl_sci[bidx];         /* my i super-cluster's index = current bidx */
    sci            = nb_sci.sci;           /* super-cluster */
    cijPackedBegin = nb_sci.cjPackedBegin; /* first ...*/
    cijPackedEnd   = nb_sci.cjPackedEnd;   /* and last index of j clusters */

    // We may need only a subset of threads active for preloading i-atoms
    // depending on the super-cluster and cluster / thread-block size.
    constexpr bool c_loadUsingAllXYThreads = (c_clSize == c_nbnxnGpuNumClusterPerSupercluster);
    if (tidxz == 0 && (c_loadUsingAllXYThreads || tidxj < c_nbnxnGpuNumClusterPerSupercluster))
    {
        /* Pre-load i-atom x and q into shared memory */
        ci = sci * c_nbnxnGpuNumClusterPerSupercluster + tidxj;
        ai = ci * c_clSize + tidxi;

        const float* shiftptr = reinterpret_cast<const float*>(&shift_vec[nb_sci.shift]);
        xqbuf = xq[ai] + make_float4(LDG(shiftptr), LDG(shiftptr + 1), LDG(shiftptr + 2), 0.0F);
        xqbuf.w *= nbparam.epsfac;
        xqib[tidxj * c_clSize + tidxi] = xqbuf;

#    ifndef LJ_COMB
        /* Pre-load the i-atom types into shared memory */
        atib[tidxj * c_clSize + tidxi] = atom_types[ai];
#    else
        /* Pre-load the LJ combination parameters into shared memory */
        ljcpib[tidxj * c_clSize + tidxi] = lj_comb[ai];
#    endif
    }

#    ifdef PRUNE_NBL
    /* Initialise one int for reducing prunedPairCount over warps */
    int* sm_prunedPairCount = reinterpret_cast<int*>(sm_nextSlotPtr);
    sm_nextSlotPtr += sizeof(*sm_prunedPairCount);
    if (tidx == 0 && tidxz == 0)
    {
        *sm_prunedPairCount = 0;
    }
    int prunedPairCount = 0;
#    endif

    __syncthreads();

    for (i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
    {
        fci_buf[i] = make_float3(0.0F);
    }

#    ifdef LJ_EWALD
    /* TODO: we are trading registers with flops by keeping lje_coeff-s, try re-calculating it later */
    lje_coeff2   = nbparam.ewaldcoeff_lj * nbparam.ewaldcoeff_lj;
    lje_coeff6_6 = lje_coeff2 * lje_coeff2 * lje_coeff2 * c_oneSixth;
#    endif


#    ifdef CALC_ENERGIES
    E_lj         = 0.0F;
    E_el         = 0.0F;

#        ifdef EXCLUSION_FORCES /* Ewald or RF */
    if (nb_sci.shift == gmx::c_centralShiftIndex
        && pl_cjPacked[cijPackedBegin].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
    {
        /* we have the diagonal: add the charge and LJ self interaction energy term */
        for (i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
        {
#            if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
            qi = xqib[i * c_clSize + tidxi].w;
            E_el += qi * qi;
#            endif

#            ifdef LJ_EWALD
            // load only the first 4 bytes of the parameter pair (equivalent with nbfp[idx].x)
            E_lj += LDG(reinterpret_cast<float*>(
                    &nbparam.nbfp[atom_types[(sci * c_nbnxnGpuNumClusterPerSupercluster + i) * c_clSize + tidxi]
                                  * (ntypes + 1)]));
#            endif
        }

        /* divide the self term(s) equally over the j-threads, then multiply with the coefficients. */
#            ifdef LJ_EWALD
        E_lj /= c_clSize * NTHREAD_Z;
        E_lj *= 0.5F * c_oneSixth * lje_coeff6_6;
#            endif

#            if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
        /* Correct for epsfac^2 due to adding qi^2 */
        E_el /= nbparam.epsfac * c_clSize * NTHREAD_Z;
#                if defined EL_RF || defined EL_CUTOFF
        E_el *= -0.5F * reactionFieldShift;
#                else
        E_el *= -beta * M_FLOAT_1_SQRTPI; /* last factor 1/sqrt(pi) */
#                endif
#            endif /* EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF */
    }
#        endif     /* EXCLUSION_FORCES */

#    endif /* CALC_ENERGIES */

#    ifdef EXCLUSION_FORCES
    // Note that we use & instead of && for performance (benchmarked in 2017)
    const int nonSelfInteraction = !(nb_sci.shift == gmx::c_centralShiftIndex & tidxj <= tidxi);
#    endif

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster;
     * The loop stride NTHREAD_Z ensures that consecutive warps-pairs are assigned
     * consecutive jPacked's entries.
     */
    for (jPacked = cijPackedBegin + tidxz; jPacked < cijPackedEnd; jPacked += NTHREAD_Z)
    {
        wexcl_idx = pl_cjPacked[jPacked].imei[widx].excl_ind;
        imask     = pl_cjPacked[jPacked].imei[widx].imask;
        wexcl     = excl[wexcl_idx].pair[(tidx) & (warp_size - 1)];

#    ifndef PRUNE_NBL
        if (imask)
#    endif
        {
            if (c_preloadCj)
            {
                /* Pre-load cj into shared memory on both warps separately */
                if ((tidxj == 0 | tidxj == 4) & (tidxi < c_nbnxnGpuJgroupSize))
                {
                    cjs[tidxi + tidxj * c_nbnxnGpuJgroupSize / c_splitClSize] =
                            pl_cjPacked[jPacked].cj[tidxi];
                }
                __syncwarp(c_fullWarpMask);
            }

#    if DO_JM_UNROLL
#        pragma unroll jmLoopUnrollFactor
#    endif
            for (jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                if (imask & (superClInteractionMask << (jm * c_nbnxnGpuNumClusterPerSupercluster)))
                {
                    mask_ji = (1U << (jm * c_nbnxnGpuNumClusterPerSupercluster));

                    cj = c_preloadCj ? cjs[jm + (tidxj & 4) * c_nbnxnGpuJgroupSize / c_splitClSize]
                                     : cj = pl_cjPacked[jPacked].cj[jm];

                    aj = cj * c_clSize + tidxj;

                    /* load j atom data */
                    xqbuf = xq[aj];
                    xj    = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f  = xqbuf.w;
#    ifndef LJ_COMB
                    typej = atom_types[aj];
#    else
                    ljcp_j = lj_comb[aj];
#    endif

                    fcj_buf = make_float3(0.0F);

#    if !defined PRUNE_NBL
#        pragma unroll c_nbnxnGpuNumClusterPerSupercluster
#    endif
                    for (i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                    {
                        if (imask & mask_ji)
                        {
                            ci = sci * c_nbnxnGpuNumClusterPerSupercluster + i; /* i cluster index */

                            /* all threads load an atom from i cluster ci into shmem! */
                            xqbuf = xqib[i * c_clSize + tidxi];
                            xi    = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

                            /* distance between i and j atoms */
                            rv = xi - xj;
                            r2 = norm2(rv);

#    ifdef PRUNE_NBL
                            /* If _none_ of the atoms pairs are in cutoff range,
                               the bit corresponding to the current
                               cluster-pair in imask gets set to 0. */
                            if (!__any_sync(c_fullWarpMask, r2 < rlist_sq))
                            {
                                imask &= ~mask_ji;
                            }
#    endif

                            int_bit = (wexcl & mask_ji) ? 1.0F : 0.0F;

                            /* cutoff & exclusion check */
#    ifdef EXCLUSION_FORCES
                            if ((r2 < rcoulomb_sq) * (nonSelfInteraction | (ci != cj)))
#    else
                            if ((r2 < rcoulomb_sq) * int_bit)
#    endif
                            {
                                /* load the rest of the i-atom parameters */
                                qi = xqbuf.w;

#    ifndef LJ_COMB
                                /* LJ 6*C6 and 12*C12 */
                                typei = atib[i * c_clSize + tidxi];
                                fetch_nbfp_c6_c12(c6, c12, nbparam, ntypes * typei + typej);
#    else
                                ljcp_i       = ljcpib[i * c_clSize + tidxi];
#        ifdef LJ_COMB_GEOM
                                c6           = ljcp_i.x * ljcp_j.x;
                                c12          = ljcp_i.y * ljcp_j.y;
#        else
                                /* LJ 2^(1/6)*sigma and 12*epsilon */
                                sigma   = ljcp_i.x + ljcp_j.x;
                                epsilon = ljcp_i.y * ljcp_j.y;
#            if defined CALC_ENERGIES || defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                                convert_sigma_epsilon_to_c6_c12(sigma, epsilon, &c6, &c12);
#            endif
#        endif /* LJ_COMB_GEOM */
#    endif     /* LJ_COMB */

                                // Ensure distance do not become so small that r^-12 overflows
                                r2 = max(r2, c_nbnxnMinDistanceSquared);

                                inv_r  = rsqrt(r2);
                                inv_r2 = inv_r * inv_r;
#    if !defined LJ_COMB_LB || defined CALC_ENERGIES
                                inv_r6 = inv_r2 * inv_r2 * inv_r2;
#        ifdef EXCLUSION_FORCES
                                /* We could mask inv_r2, but with Ewald
                                 * masking both inv_r6 and F_invr is faster */
                                inv_r6 *= int_bit;
#        endif /* EXCLUSION_FORCES */

                                F_invr = inv_r6 * (c12 * inv_r6 - c6) * inv_r2;
#        if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                                E_lj_p = int_bit
                                         * (c12 * (inv_r6 * inv_r6 + nbparam.repulsion_shift.cpot) * c_oneTwelfth
                                            - c6 * (inv_r6 + nbparam.dispersion_shift.cpot) * c_oneSixth);
#        endif
#    else /* !LJ_COMB_LB || CALC_ENERGIES */
                                float sig_r  = sigma * inv_r;
                                float sig_r2 = sig_r * sig_r;
                                float sig_r6 = sig_r2 * sig_r2 * sig_r2;
#        ifdef EXCLUSION_FORCES
                                sig_r6 *= int_bit;
#        endif /* EXCLUSION_FORCES */

                                F_invr = epsilon * sig_r6 * (sig_r6 - 1.0F) * inv_r2;
#    endif     /* !LJ_COMB_LB || CALC_ENERGIES */

#    ifdef LJ_FORCE_SWITCH
#        ifdef CALC_ENERGIES
                                calculate_force_switch_F_E(nbparam, c6, c12, inv_r, r2, &F_invr, &E_lj_p);
#        else
                                calculate_force_switch_F(nbparam, c6, c12, inv_r, r2, &F_invr);
#        endif /* CALC_ENERGIES */
#    endif     /* LJ_FORCE_SWITCH */


#    ifdef LJ_EWALD
#        ifdef LJ_EWALD_COMB_GEOM
#            ifdef CALC_ENERGIES
                                calculate_lj_ewald_comb_geom_F_E(
                                        nbparam, typei, typej, r2, inv_r2, lje_coeff2, lje_coeff6_6, int_bit, &F_invr, &E_lj_p);
#            else
                                calculate_lj_ewald_comb_geom_F(
                                        nbparam, typei, typej, r2, inv_r2, lje_coeff2, lje_coeff6_6, &F_invr);
#            endif /* CALC_ENERGIES */
#        elif defined LJ_EWALD_COMB_LB
                                calculate_lj_ewald_comb_LB_F_E(nbparam,
                                                               typei,
                                                               typej,
                                                               r2,
                                                               inv_r2,
                                                               lje_coeff2,
                                                               lje_coeff6_6,
#            ifdef CALC_ENERGIES
                                                               int_bit,
                                                               &F_invr,
                                                               &E_lj_p
#            else
                                                               0,
                                                               &F_invr,
                                                               nullptr
#            endif /* CALC_ENERGIES */
                                );
#        endif     /* LJ_EWALD_COMB_GEOM */
#    endif         /* LJ_EWALD */

#    ifdef LJ_POT_SWITCH
#        ifdef CALC_ENERGIES
                                calculate_potential_switch_F_E(nbparam, inv_r, r2, &F_invr, &E_lj_p);
#        else
                                calculate_potential_switch_F(nbparam, inv_r, r2, &F_invr, &E_lj_p);
#        endif /* CALC_ENERGIES */
#    endif     /* LJ_POT_SWITCH */

#    ifdef VDW_CUTOFF_CHECK
                                /* Separate VDW cut-off check to enable twin-range cut-offs
                                 * (rvdw < rcoulomb <= rlist)
                                 */
                                vdw_in_range = (r2 < rvdw_sq) ? 1.0F : 0.0F;
                                F_invr *= vdw_in_range;
#        ifdef CALC_ENERGIES
                                E_lj_p *= vdw_in_range;
#        endif
#    endif /* VDW_CUTOFF_CHECK */

#    ifdef CALC_ENERGIES
                                E_lj += E_lj_p;
#    endif


#    ifdef EL_CUTOFF
#        ifdef EXCLUSION_FORCES
                                F_invr += qi * qj_f * int_bit * inv_r2 * inv_r;
#        else
                                F_invr += qi * qj_f * inv_r2 * inv_r;
#        endif
#    endif
#    ifdef EL_RF
                                F_invr += qi * qj_f * (int_bit * inv_r2 * inv_r - two_k_rf);
#    endif
#    if defined   EL_EWALD_ANA
                                F_invr += qi * qj_f
                                          * (int_bit * inv_r2 * inv_r + pmecorrF(beta2 * r2) * beta3);
#    elif defined EL_EWALD_TAB
                                F_invr += qi * qj_f
                                          * (int_bit * inv_r2
                                             - interpolate_coulomb_force_r(nbparam, r2 * inv_r))
                                          * inv_r;
#    endif /* EL_EWALD_ANA/TAB */

#    ifdef CALC_ENERGIES
#        ifdef EL_CUTOFF
                                E_el += qi * qj_f * (int_bit * inv_r - reactionFieldShift);
#        endif
#        ifdef EL_RF
                                E_el += qi * qj_f
                                        * (int_bit * inv_r + 0.5F * two_k_rf * r2 - reactionFieldShift);
#        endif
#        ifdef EL_EWALD_ANY
                                /* 1.0F - erff is faster than erfcf */
                                E_el += qi * qj_f
                                        * (inv_r * (int_bit - erff(r2 * inv_r * beta)) - int_bit * ewald_shift);
#        endif /* EL_EWALD_ANY */
#    endif
                                f_ij = rv * F_invr;

                                /* accumulate j forces in registers */
                                fcj_buf -= f_ij;

                                /* accumulate i forces in registers */
                                fci_buf[i] += f_ij;
                            }
                        }

                        /* shift the mask bit by 1 */
                        mask_ji += mask_ji;
                    }

                    /* reduce j forces */
                    reduce_force_j_warp_shfl(fcj_buf, f, tidxi, aj, c_fullWarpMask);
                }
            }
#    ifdef PRUNE_NBL
            /* Update the imask with the new one which does not contain the
               out of range clusters anymore. */
            pl_cjPacked[jPacked].imei[widx].imask = imask;
            prunedPairCount += __popc(imask);
#    endif
        }
        if (c_preloadCj)
        {
            // avoid shared memory WAR hazards on sm_cjs between loop iterations
            __syncwarp(c_fullWarpMask);
        }
    }

    /* skip central shifts when summing shift forces */
    if (nb_sci.shift == gmx::c_centralShiftIndex)
    {
        bCalcFshift = false;
    }

    float fshift_buf = 0.0F;

    /* reduce i forces */
    for (i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
    {
        ai = (sci * c_nbnxnGpuNumClusterPerSupercluster + i) * c_clSize + tidxi;
        reduce_force_i_warp_shfl(fci_buf[i], f, &fshift_buf, bCalcFshift, tidxj, ai, c_fullWarpMask);
    }

    /* add up local shift forces into global mem, tidxj indexes x,y,z */
    if (bCalcFshift && (tidxj & 3) < 3)
    {
        float3* fShift = asFloat3(atdat.fShift);
        atomicAdd(&(fShift[nb_sci.shift].x) + (tidxj & 3), fshift_buf);
    }

#    ifdef CALC_ENERGIES
    /* reduce the energies over warps and store into global memory */
    reduce_energy_warp_shfl(E_lj, E_el, e_lj, e_el, tidx, c_fullWarpMask);
#    endif


#    ifdef PRUNE_NBL
    /* aggregate neighbour counts, to be used in bucket sci sort */
    /* One thread in each warp contributes the count for that warp as soon as it reaches here.
     * Masks are calculated per warp in a warp synchronising operation, so no syncthreads
     * required here. */
    if (tidxInWarp == 0)
    {
        atomicAdd(sm_prunedPairCount, prunedPairCount);
    }
    __syncthreads();
    prunedPairCount = *sm_prunedPairCount;
    if (tidxi == 0 && tidxj == 0 && tidxz == 0)
    {
        /* one thread in the block writes the final count for this sci */
        int  index            = max(c_sciHistogramSize - prunedPairCount - 1, 0);
        int* pl_sci_histogram = plist.sorting.sciHistogram;
        atomicAdd(pl_sci_histogram + index, 1);
        int* pl_sci_count  = plist.sorting.sciCount;
        pl_sci_count[bidx] = index;
    }
#    endif
}
#endif /* FUNCTION_DECLARATION_ONLY */

#undef NTHREAD_Z
#undef MIN_BLOCKS_PER_MP
#undef THREADS_PER_BLOCK

#undef EL_EWALD_ANY
#undef EXCLUSION_FORCES
#undef LJ_EWALD

#undef LJ_COMB

} // namespace gmx
