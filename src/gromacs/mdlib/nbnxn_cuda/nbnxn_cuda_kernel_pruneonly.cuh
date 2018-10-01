/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

/*! \internal \file
 *  \brief
 *  CUDA non-bonded prune-only kernel.
 *
 *  Unlike the non-bonded interaction kernels, this is not preprocessor-generated,
 *  the two flavors achieved by templating.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Berk Hess <hess@kth.se>
 *  \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/math/utilities.h"
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_utils.cuh"
#include "gromacs/pbcutil/ishift.h"

#include "nbnxn_cuda_types.h"

/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

/**@{*/
/*! \brief Compute capability dependent definition of kernel launch configuration parameters.
 *
 * Kernel launch bounds for different compute capabilities. The value of NTHREAD_Z
 * represents the j-concurrency, hence it determines the number of threads per block.
 * It is chosen such that 100% occupancy is maintained (on Maxwell and later for any NTHREAD_Z,
 * requires >=4 warp/block, NTHREAD_Z>=2 on Kepler).
 *
 * Hence, values NTHREAD_Z >= 2 trade inter- for intra-block parallelism
 * which has the advantage of lowering the overhead of starting up a block, filling shmem
 * and registers, etc. Ideally we'd want to expose as much intra-block work as possible
 * As we also split lists to cater for the block-parallelization needed by the register-
 * limited non-bonded kernels, for very short j-loops large NTHREAD_Z will cause slowdown
 * as it leads to intra-block warp imbalance. Ideally, we'd want to auto-tune the choice
 * of NTHREAD_Z, but for now we instead pick a reasonable tradeoff-value.
 *
 * Note that given the above input size tradeoffs and that performance depends on
 * additional factors including GPU arch, #SM's, we'll accept performance tradeoffs
 * of using a fixed NTHREAD_Z=4. The following outliers have been observed:
 *   - up to 25% faster (rolling) prune kernels with NTHREAD_Z=8 in the regime where lists
 *     are not split (much), but the rolling chunks are small;
 *   - with large inputs NTHREAD_Z=1 is 2-3% faster (on CC>=5.0)
 */
#define NTHREAD_Z           (GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY)
#define THREADS_PER_BLOCK   (c_clSize*c_clSize*NTHREAD_Z)
// we want 100% occupancy, so max threads/block
#define MIN_BLOCKS_PER_MP   (GMX_CUDA_MAX_THREADS_PER_MP/THREADS_PER_BLOCK)
/**@}*/

/*! \brief Nonbonded list pruning kernel.
 *
 *  The \p haveFreshList template parameter defines the two flavors of the kernel; when
 *  true a new list from immediately after pair-list generation is pruned using rlistOuter,
 *  the pruned masks are stored in a separate buffer and the outer-list is pruned
 *  using the rlistInner distance; when false only the pruning with rlistInner is performed.
 *
 *  Kernel launch parameters:
 *   - #blocks   = #pair lists, blockId = pair list Id
 *   - #threads  = NTHREAD_Z * c_clSize^2
 *   - shmem     = see nbnxn_cuda.cu:calc_shmem_required_prune()
 *
 *   Each thread calculates an i-j atom distance..
 */
template <bool haveFreshList>
__launch_bounds__(THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)
__global__ void nbnxn_kernel_prune_cuda(const cu_atomdata_t atdat,
                                        const cu_nbparam_t  nbparam,
                                        const cu_plist_t    plist,
                                        int                 numParts,
                                        int                 part)
#ifdef FUNCTION_DECLARATION_ONLY
;     /* Only do function declaration, omit the function body. */
#else
{

    /* convenience variables */
    const nbnxn_sci_t  *pl_sci      = plist.sci;
    nbnxn_cj4_t        *pl_cj4      = plist.cj4;
    const float4       *xq          = atdat.xq;
    const float3       *shift_vec   = atdat.shift_vec;

    float               rlistOuter_sq = nbparam.rlistOuter_sq;
    float               rlistInner_sq = nbparam.rlistInner_sq;

    /* thread/block/warp id-s */
    unsigned int tidxi  = threadIdx.x;
    unsigned int tidxj  = threadIdx.y;
#if NTHREAD_Z == 1
    unsigned int tidxz  = 0;
#else
    unsigned int tidxz  = threadIdx.z;
#endif
    unsigned int bidx   = blockIdx.x;
    unsigned int widx   = (threadIdx.y * c_clSize) / warp_size; /* warp index */

    /*********************************************************************
     * Set up shared memory pointers.
     * sm_nextSlotPtr should always be updated to point to the "next slot",
     * that is past the last point where data has been stored.
     */
    extern __shared__  char sm_dynamicShmem[];
    char                   *sm_nextSlotPtr = sm_dynamicShmem;
    static_assert(sizeof(char) == 1, "The shared memory offset calculation assumes that char is 1 byte");

    /* shmem buffer for i x+q pre-loading */
    float4 *xib     = (float4 *)sm_nextSlotPtr;
    sm_nextSlotPtr += (c_numClPerSupercl * c_clSize * sizeof(*xib));

    /* shmem buffer for cj, for each warp separately */
    int *cjs        = (int *)(sm_nextSlotPtr);
    /* the cjs buffer's use expects a base pointer offset for pairs of warps in the j-concurrent execution */
    cjs            += tidxz * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize;
    sm_nextSlotPtr += (NTHREAD_Z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(*cjs));
    /*********************************************************************/


    nbnxn_sci_t nb_sci      = pl_sci[bidx*numParts + part]; /* my i super-cluster's index = sciOffset + current bidx * numParts + part */
    int         sci         = nb_sci.sci;                   /* super-cluster */
    int         cij4_start  = nb_sci.cj4_ind_start;         /* first ...*/
    int         cij4_end    = nb_sci.cj4_ind_end;           /* and last index of j clusters */

    if (tidxz == 0)
    {
        /* Pre-load i-atom x and q into shared memory */
        int ci = sci * c_numClPerSupercl + tidxj;
        int ai = ci * c_clSize + tidxi;

        /* We don't need q, but using float4 in shmem avoids bank conflicts.
           (but it also wastes L2 bandwidth). */
        float4 tmp = xq[ai];
        float4 xi  = tmp + shift_vec[nb_sci.shift];
        xib[tidxj * c_clSize + tidxi] = xi;
    }
    __syncthreads();

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster;
     * The loop stride NTHREAD_Z ensures that consecutive warps-pairs are assigned
     * consecutive j4's entries.
     */
    for (int j4 = cij4_start + tidxz; j4 < cij4_end; j4 += NTHREAD_Z)
    {
        unsigned int imaskFull, imaskCheck, imaskNew;

        if (haveFreshList)
        {
            /* Read the mask from the list transferred from the CPU */
            imaskFull = pl_cj4[j4].imei[widx].imask;
            /* We attempt to prune all pairs present in the original list */
            imaskCheck = imaskFull;
            imaskNew   = 0;
        }
        else
        {
            /* Read the mask from the "warp-pruned" by rlistOuter mask array */
            imaskFull = plist.imask[j4*c_nbnxnGpuClusterpairSplit + widx];
            /* Read the old rolling pruned mask, use as a base for new */
            imaskNew = pl_cj4[j4].imei[widx].imask;
            /* We only need to check pairs with different mask */
            imaskCheck = (imaskNew ^ imaskFull);
        }

        if (imaskCheck)
        {
            /* Pre-load cj into shared memory on both warps separately */
            if ((tidxj == 0 || tidxj == 4) && tidxi < c_nbnxnGpuJgroupSize)
            {
                cjs[tidxi + tidxj * c_nbnxnGpuJgroupSize/c_splitClSize] = pl_cj4[j4].cj[tidxi];
            }
            gmx_syncwarp(c_fullWarpMask);

#pragma unroll 4
            for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                if (imaskCheck & (superClInteractionMask << (jm * c_numClPerSupercl)))
                {
                    unsigned int mask_ji = (1U << (jm * c_numClPerSupercl));

                    int          cj      = cjs[jm + (tidxj & 4) * c_nbnxnGpuJgroupSize/c_splitClSize];
                    int          aj      = cj * c_clSize + tidxj;

                    /* load j atom data */
                    float4 tmp  = xq[aj];
                    float3 xj   = make_float3(tmp.x, tmp.y, tmp.z);

#pragma unroll 8
                    for (int i = 0; i < c_numClPerSupercl; i++)
                    {
                        if (imaskCheck & mask_ji)
                        {
                            /* load i-cluster coordinates from shmem */
                            float4 xi = xib[i * c_clSize + tidxi];


                            /* distance between i and j atoms */
                            float3 rv = make_float3(xi.x, xi.y, xi.z) - xj;
                            float  r2 = norm2(rv);

                            /* If _none_ of the atoms pairs are in rlistOuter
                               range, the bit corresponding to the current
                               cluster-pair in imask gets set to 0. */
                            if (haveFreshList && !gmx_any_sync(c_fullWarpMask, r2 < rlistOuter_sq))
                            {
                                imaskFull &= ~mask_ji;
                            }
                            /* If any atom pair is within range, set the bit
                               corresponding to the current cluster-pair. */
                            if (gmx_any_sync(c_fullWarpMask, r2 < rlistInner_sq))
                            {
                                imaskNew |= mask_ji;
                            }
                        }

                        /* shift the mask bit by 1 */
                        mask_ji += mask_ji;
                    }
                }
            }

            if (haveFreshList)
            {
                /* copy the list pruned to rlistOuter to a separate buffer */
                plist.imask[j4*c_nbnxnGpuClusterpairSplit + widx] = imaskFull;
            }
            /* update the imask with only the pairs up to rlistInner */
            plist.cj4[j4].imei[widx].imask = imaskNew;
        }
        // avoid shared memory WAR hazards between loop iterations
        gmx_syncwarp(c_fullWarpMask);
    }
}
#endif /* FUNCTION_DECLARATION_ONLY */

#undef NTHREAD_Z
#undef MIN_BLOCKS_PER_MP
#undef THREADS_PER_BLOCK
