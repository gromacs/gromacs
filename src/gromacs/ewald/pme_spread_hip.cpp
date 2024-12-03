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

/*! \internal \file
 *  \brief Implements PME GPU spline calculation and charge spreading in HIP.
 *  TODO: consider always pre-sorting particles (as in DD case).
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>
#include <climits>

#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/hip_sycl_kernel_utils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"

#include "pme_gpu_calculate_splines_hip.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_types.h"

static constexpr int sc_spreadHipMaxWarpsPerBlock = 8;

template<int parallelExecutionWidth>
static constexpr int sc_spreadMaxThreadsPerBlock = sc_spreadHipMaxWarpsPerBlock * parallelExecutionWidth;

/*! \brief Per-axis sub-grid tile dimensions for the LDS accumulation path.
 *
 * When \c ThreadsPerAtom::Order is selected, each wavefront of the
 * spline-and-spread kernel accumulates its atoms' contributions into a private
 * sub-grid in shared memory, then cooperatively flushes the sub-grid to the
 * real-space grid in global memory. This replaces many global atomicAdds with
 * (cheaper) LDS atomicAdds, plus one global atomicAdd per non-zero sub-grid
 * cell.
 *
 * The wavefront bounding-box of its atoms' gridline bases must fit in a tile
 * of size \c (TileDimX, TileDimY, TileDimZ); waves that do not fit fall back
 * to the direct per-lane atomicAdd path at the cost of one warp-wide int
 * min/max reduction.
 *
 * With PME order 4 an atom's support spans \c order cells along each axis,
 * so an axis dimension \c D accepts any wave whose atoms' gridline-base span
 * along that axis is at most \c D-order. Larger D catches more waves on the
 * fast path at the cost of LDS (and thus occupancy).
 *
 * The current default (10, 10, 8) was selected empirically on MI300X. Two
 * effects drive the choice:
 *
 *  - Fit rate: tolerated atom-base span per axis is \c D-order = (6, 6, 4)
 *    cells, i.e. 50 % more X/Y slack than the original (8, 8, 8) tile.
 *  - LDS bank conflicts: with 32 LDS banks the row-major starting bank for
 *    atom k is \c ((lx*TileY + ly)*TileZ + base_lz) mod 32. A cubic tile of
 *    side that shares a factor of 32 (8, 16) zeros out the lx contribution
 *    and serializes accesses across atoms with the same ly. (10, 10, 10) and
 *    (10, 10, 8) avoid that aliasing because 10 is coprime to 32/gcd(10,32).
 *    Padding the Z stride to an odd value to break the residual aliasing
 *    against ly mod 4 was tried and gave no measurable wall-clock gain on
 *    our benchmarks, so the LDS layout stays unpadded.
 *
 * Empirically (10, 10, 8) gave significant gains on PME-bound STMV and
 * cellulose vs (10, 10, 10), and (10, 10, 10) plateaus before (12, 12, 12).
 * Smaller-system regressions (e.g. benchpep-h) are modest. The "best" shape
 * is in principle grid-shape-dependent and could later be promoted to a
 * runtime-tuned kernel template parameter; for now we ship a single default.
 *
 * Footprint (LDS = product * 4 bytes/wave; one tile per wave of the block):
 *   (10, 10, 8)  ->  3.13 KB/wave, 25.0 KB/block (current default)
 *   (10, 10, 10) ->  3.90 KB/wave, 31.2 KB/block
 *   (12, 12, 12) ->  6.75 KB/wave, 54.0 KB/block (tight, hurts occupancy)
 *   (8, 8, 8)    ->  2.00 KB/wave, 16.0 KB/block (previous default)
 *
 * Note that the row-major LDS layout is (X, Y, Z) with Z fastest, so the
 * cooperative flush strides along Z, matching the global-grid layout
 * (gridIndexGlobal = ix*pny*pnz + iy*pnz + iz) for coalesced per-lane
 * global atomic addresses.
 */
static constexpr int c_pmeHipSpreadLdsTileDimX = 10;
static constexpr int c_pmeHipSpreadLdsTileDimY = 10;
static constexpr int c_pmeHipSpreadLdsTileDimZ = 8;
static constexpr int c_pmeHipSpreadLdsTileSize =
        c_pmeHipSpreadLdsTileDimX * c_pmeHipSpreadLdsTileDimY * c_pmeHipSpreadLdsTileDimZ;

/*! \brief
 * Charge spreading onto the grid.
 * This corresponds to the CPU function spread_coefficients_bsplines_thread().
 * Optional second stage of the spline_and_spread_kernel.
 *
 * When \c threadsPerAtom is \c ThreadsPerAtom::Order, the wave first tries
 * to accumulate into a private sub-grid tile in \p sm_subGrid (sized
 * \c c_pmeHipSpreadLdsTileDim{X,Y,Z}, one tile per wave of the block) and
 * flushes the tile to \c gm_grid at the end. If the bounding box of the
 * wave's atoms' gridline bases does not fit into the tile, the wave falls
 * back to the direct per-lane atomicAdd path (identical to the original
 * implementation).
 *
 * All 64/32 lanes of the wave must enter this function; lanes whose
 * \p laneValid is \c false contribute \c INT_MAX / \c INT_MIN sentinels to
 * the bounding-box reduction and do not issue any atomicAdd. Keeping the
 * invocation wave-uniform is required for the \c __shfl_xor-based bbox
 * reduction to see well-defined values on every lane.
 *
 * \tparam     order                PME interpolation order.
 * \tparam     wrapX                Whether the grid overlap in dimension X should be wrapped.
 * \tparam     wrapY                Whether the grid overlap in dimension Y should be wrapped.
 * \tparam     gridIndex            The index of the grid to use in the kernel.
 * \tparam     threadsPerAtom       How many threads work on each atom
 * \tparam     parallelExecutionWidth  Wavefront width (32 or 64).
 *
 * \param[in]     kernelParams        Input PME HIP data in constant memory.
 * \param[in]     atomCharge          Atom charge/coefficient of atom processed by thread.
 * \param[in]     sm_gridlineIndices  Atom gridline indices in shared memory.
 * \param[in]     sm_theta            Atom spline values in shared memory.
 * \param[in,out] sm_subGrid          Per-block LDS sub-grid arena, sized at
 *                                    least \c numWaves * c_pmeHipSpreadLdsTileSize
 *                                    floats when the LDS path is enabled;
 *                                    unused otherwise (may be \c nullptr / stub).
 * \param[in]     laneValid           Whether this lane's atom index is within
 *                                    bounds and should contribute to the grid.
 */
template<int order, bool wrapX, bool wrapY, int gridIndex, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
__device__ __forceinline__ void spread_charges(const PmeGpuKernelParams kernelParams,
                                               const float*             atomCharge,
                                               const int* __restrict__ sm_gridlineIndices,
                                               const float* __restrict__ sm_theta,
                                               float* __restrict__ sm_subGrid,
                                               const bool laneValid)
{
    /* Global memory pointer to the output grid */
    float* __restrict__ gm_grid = kernelParams.grid.d_realGrid[gridIndex];

    // Number of atoms processed by a single warp in spread and gather
    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int atomsPerWarp = parallelExecutionWidth / threadsPerAtomValue;

    const int nx  = kernelParams.grid.realGridSize[XX];
    const int ny  = kernelParams.grid.realGridSize[YY];
    const int nz  = kernelParams.grid.realGridSize[ZZ];
    const int pny = kernelParams.grid.realGridSizePadded[YY];
    const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

    const int offx = 0, offy = 0, offz = 0; // unused for now

    const int atomIndexLocal = threadIdx.z;
    const int ithz           = threadIdx.x;

    /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
    const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = atomIndexLocal / atomsPerWarp;

    const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);

    // Gridline-base reads; may be uninitialised for !laneValid lanes, so we
    // only feed them to the bbox reduction via sentinels below.
    const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX] - offx;
    const int iyBase = sm_gridlineIndices[atomIndexLocal * DIM + YY] - offy;
    const int izBase = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] - offz;

    int iz = izBase + ithz;
    if (iz >= nz)
    {
        iz -= nz;
    }

    const int   splineIndexZ = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
    const float thetaZ       = sm_theta[splineIndexZ];

    // NB: the short-circuit here is load-bearing -- for !laneValid lanes the
    // dereferenced *atomCharge is a garbage value from the padded allocation,
    // and the assert inside pme_gpu_check_atom_charge would fire in debug builds.
    const bool active = laneValid && pme_gpu_check_atom_charge(*atomCharge);

    if constexpr (threadsPerAtom == ThreadsPerAtom::Order)
    {
        constexpr int tileDimX = c_pmeHipSpreadLdsTileDimX;
        constexpr int tileDimY = c_pmeHipSpreadLdsTileDimY;
        constexpr int tileDimZ = c_pmeHipSpreadLdsTileDimZ;
        constexpr int tileSize = c_pmeHipSpreadLdsTileSize;

        // Wave-wide bounding box of the atoms' gridline bases. Inactive
        // lanes seed the reduction with sentinels that are neutral for
        // min/max so they do not stretch the bbox. If *every* lane in the
        // wave is inactive the bbox degenerates (ixMinBase > ixMaxBase);
        // that's fine because bboxFits below evaluates to false and we
        // take the (no-op) fallback.
        int ixMinBase = active ? ixBase : INT_MAX;
        int ixMaxBase = active ? ixBase : INT_MIN;
        int iyMinBase = active ? iyBase : INT_MAX;
        int iyMaxBase = active ? iyBase : INT_MIN;
        int izMinBase = active ? izBase : INT_MAX;
        int izMaxBase = active ? izBase : INT_MIN;

        // Fused six-way wave reduction. Each round runs all six independent
        // min/max chains in parallel, so the dependency-chain length is
        // log2(parallelExecutionWidth) shuffle+op steps total, instead of
        // 6 * log2(parallelExecutionWidth) for six sequential warp reductions.
        // This exposes the ILP to the compiler. We avoid rocprim::warp_reduce
        // because its hard VirtualWaveSize <= physical-wave-size static_assert
        // prevents the wave64 kernel from compiling on RDNA targets in this
        // multi-arch TU.
        //
        // We tried several DPP-based variants on top of this (xor 1 / xor 2
        // via quad_perm; a full DPP path using row_shr + row_bcast on CDNA
        // and row_xmask on RDNA with a final readlane broadcast) but the
        // measured benefit on PME spread wasn't large enough to justify the
        // extra ~150 lines and the wave32/wave64 code divergence. The plain
        // __shfl_xor loop is what survived.
#pragma unroll
        for (int offset = parallelExecutionWidth / 2; offset > 0; offset /= 2)
        {
            const int ixMinOther = __shfl_xor(ixMinBase, offset, parallelExecutionWidth);
            const int ixMaxOther = __shfl_xor(ixMaxBase, offset, parallelExecutionWidth);
            const int iyMinOther = __shfl_xor(iyMinBase, offset, parallelExecutionWidth);
            const int iyMaxOther = __shfl_xor(iyMaxBase, offset, parallelExecutionWidth);
            const int izMinOther = __shfl_xor(izMinBase, offset, parallelExecutionWidth);
            const int izMaxOther = __shfl_xor(izMaxBase, offset, parallelExecutionWidth);
            ixMinBase            = (ixMinOther < ixMinBase) ? ixMinOther : ixMinBase;
            ixMaxBase            = (ixMaxOther > ixMaxBase) ? ixMaxOther : ixMaxBase;
            iyMinBase            = (iyMinOther < iyMinBase) ? iyMinOther : iyMinBase;
            iyMaxBase            = (iyMaxOther > iyMaxBase) ? iyMaxOther : iyMaxBase;
            izMinBase            = (izMinOther < izMinBase) ? izMinOther : izMinBase;
            izMaxBase            = (izMaxOther > izMaxBase) ? izMaxOther : izMaxBase;
        }

        const int bboxX = ixMaxBase - ixMinBase + order;
        const int bboxY = iyMaxBase - iyMinBase + order;
        const int bboxZ = izMaxBase - izMinBase + order;

        // bboxFits is wave-uniform (all lanes see the same reduced min/max).
        const bool bboxFits = (bboxX <= tileDimX) && (bboxY <= tileDimY) && (bboxZ <= tileDimZ);

        if (bboxFits)
        {
            float* const sm_waveTile = sm_subGrid + warpIndex * tileSize;

            // Lane id within the wave. For ThreadsPerAtom::Order the block
            // has blockDim = (order, 1, atomsPerBlock), so lane within the
            // wave is atomWarpIndex*order + ithz and ranges over
            // [0, parallelExecutionWidth).
            const int laneId = atomWarpIndex * order + ithz;

            const int bboxSize = bboxX * bboxY * bboxZ;

            // Cooperatively zero the wave's tile
#pragma unroll
            for (int i = laneId; i < tileSize; i += parallelExecutionWidth)
            {
                sm_waveTile[i] = 0.0f;
            }
            __builtin_amdgcn_wave_barrier();

            // Spread into LDS. Only active lanes contribute.
            if (active)
            {
                const int lxBase = ixBase - ixMinBase;
                const int lyBase = iyBase - iyMinBase;
                const int lz     = (izBase + ithz) - izMinBase;

#pragma unroll
                for (int ithy = 0; ithy < order; ithy++)
                {
                    const int splineIndexY =
                            getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
                    const float thetaY = sm_theta[splineIndexY];
                    const float Val    = thetaZ * thetaY * (*atomCharge);
                    GMX_DEVICE_ASSERT(isfinite(Val));
                    const int ly = lyBase + ithy;

#pragma unroll
                    for (int ithx = 0; ithx < order; ithx++)
                    {
                        const int lx = lxBase + ithx;
                        const int splineIndexX =
                                getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
                        const float thetaX = sm_theta[splineIndexX];
                        GMX_DEVICE_ASSERT(isfinite(thetaX));
                        const int ldsIdx = (lx * tileDimY + ly) * tileDimZ + lz;
                        atomicAdd(sm_waveTile + ldsIdx, thetaX * Val);
                    }
                }
            }
            __builtin_amdgcn_wave_barrier();

            // Cooperative flush: one atomicAdd per non-zero cell of the
            // wave's bbox-cropped tile. Global-wrap is applied here
            // (tile coordinates are always unwrapped).
            for (int i = laneId; i < bboxSize; i += parallelExecutionWidth)
            {
                const int lz  = i % bboxZ;
                const int lyx = i / bboxZ;
                const int ly  = lyx % bboxY;
                const int lx  = lyx / bboxY;

                const int   ldsIdx = (lx * tileDimY + ly) * tileDimZ + lz;
                const float val    = sm_waveTile[ldsIdx];
                if (val != 0.0F)
                {
                    int ix  = ixMinBase + lx;
                    int iy  = iyMinBase + ly;
                    int izG = izMinBase + lz;
                    if constexpr (wrapX)
                    {
                        if (ix >= nx)
                        {
                            ix -= nx;
                        }
                    }
                    if constexpr (wrapY)
                    {
                        if (iy >= ny)
                        {
                            iy -= ny;
                        }
                    }
                    if (izG >= nz)
                    {
                        izG -= nz;
                    }
                    const int gridIndexGlobal = ix * pny * pnz + iy * pnz + izG;
                    atomicAdd(gm_grid + gridIndexGlobal, val);
                }
            }
            return;
        }
        // bbox didn't fit -> fall through to the direct-atomic path below.
    }

    /* Direct-atomic fallback path (identical semantics to the original kernel). */
    if (active)
    {
        /* loop not used if order*order threads per atom */
        const int ithyMin = (threadsPerAtom == ThreadsPerAtom::Order) ? 0 : threadIdx.y;
        const int ithyMax = (threadsPerAtom == ThreadsPerAtom::Order) ? order : threadIdx.y + 1;
        for (int ithy = ithyMin; ithy < ithyMax; ithy++)
        {
            int iy = iyBase + ithy;
            if (wrapY & (iy >= ny))
            {
                iy -= ny;
            }

            const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
            float       thetaY = sm_theta[splineIndexY];
            const float Val    = thetaZ * thetaY * (*atomCharge);
            GMX_DEVICE_ASSERT(isfinite(Val));
            const int offset = iy * pnz + iz;

#pragma unroll
            for (int ithx = 0; (ithx < order); ithx++)
            {
                int ix = ixBase + ithx;
                if (wrapX & (ix >= nx))
                {
                    ix -= nx;
                }
                const int gridIndexGlobal = ix * pny * pnz + offset;
                const int splineIndexX =
                        getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
                const float thetaX = sm_theta[splineIndexX];
                GMX_DEVICE_ASSERT(isfinite(thetaX));
                GMX_DEVICE_ASSERT(isfinite(gm_grid[gridIndexGlobal]));
                atomicAdd(gm_grid + gridIndexGlobal, thetaX * Val);
            }
        }
    }
}

/*! \brief
 * A spline computation and charge spreading kernel function.
 *
 * Two tuning parameters can be used for additional performance. For small systems and for debugging
 * writeGlobal should be used removing the need to recalculate the theta values in the gather kernel.
 * Similarly for useOrderThreads large systems order threads per atom gives higher performance than order*order threads
 *
 * \tparam     order                PME interpolation order.
 * \tparam     computeSplines       A boolean which tells if the spline parameter and
 *                                  gridline indices' computation should be performed.
 * \tparam     spreadCharges        A boolean which tells if the charge spreading should be performed.
 * \tparam     wrapX                A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam     wrapY                A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \tparam     numGrids             The number of grids to use in the kernel. Can be 1 or 2.
 * \tparam     writeGlobal          A boolean which tells if the theta values and gridlines should be written to global memory.
 * \tparam     threadsPerAtom       How many threads work on each atom
 * param[in]  kernelParams         Input PME HIP data in constant memory.
 */
template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
LAUNCH_BOUNDS_EXACT_SINGLE(sc_spreadMaxThreadsPerBlock<parallelExecutionWidth>)
__global__ void pmeSplineAndSpreadKernel(const PmeGpuKernelParams kernelParams)
{
    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int atomsPerBlock = sc_spreadMaxThreadsPerBlock<parallelExecutionWidth> / threadsPerAtomValue;
    // Number of atoms processed by a single warp in spread and gather
    constexpr int atomsPerWarp = parallelExecutionWidth / threadsPerAtomValue;
    // Gridline indices, ivec
    __shared__ int sm_gridlineIndices[atomsPerBlock * DIM];
    // Charges
    __shared__ float sm_coefficients[atomsPerBlock];
    // Spline values
    __shared__ float sm_theta[atomsPerBlock * DIM * order];
    /* Per-wave sub-grid arena for the LDS accumulation path. Sized to one
     * tile per wave of the block when the optimization is active; otherwise
     * stubbed out at 1 float so the symbol still exists for the call site. */
    constexpr bool useLdsSubGrid = spreadCharges && threadsPerAtom == ThreadsPerAtom::Order;
    constexpr int  subGridStorage =
            useLdsSubGrid ? (sc_spreadHipMaxWarpsPerBlock * c_pmeHipSpreadLdsTileSize) : 1;
    __shared__ float sm_subGrid[subGridStorage];
    float            dtheta;

    float3 atomX;
    float  atomCharge;

    const int blockIndex      = blockIdx.y * gridDim.x + blockIdx.x;
    const int atomIndexOffset = blockIndex * atomsPerBlock + kernelParams.pipelineAtomStart;

    /* Thread index w.r.t. block */
    const int threadLocalId =
            (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / parallelExecutionWidth;

    /* Atom index w.r.t. warp */
    const int atomWarpIndex = threadIdx.z % atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;
    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    /* Early return for fully empty blocks at the end
     * (should only happen for billions of input atoms)
     */
    if (atomIndexOffset >= kernelParams.atoms.nAtoms)
    {
        return;
    }
    /* Charges, required for both spline and spread */
    if constexpr (c_useAtomDataPrefetch)
    {
        pme_gpu_stage_atom_data<float, atomsPerBlock, 1, !computeSplines>(
                sm_coefficients, &kernelParams.atoms.d_coefficients[0][kernelParams.pipelineAtomStart]);
        __syncthreads();
        atomCharge = sm_coefficients[atomIndexLocal];
    }
    else
    {
        atomCharge = kernelParams.atoms.d_coefficients[0][atomIndexGlobal];
    }

    if constexpr (computeSplines)
    {
        const float3* __restrict__ gm_coordinates = asFloat3(kernelParams.atoms.d_coordinates);
        if constexpr (c_useAtomDataPrefetch)
        {
            // Coordinates
            __shared__ float3 sm_coordinates[atomsPerBlock];

            /* Staging coordinates */
            pme_gpu_stage_atom_data<float3, atomsPerBlock, 1, !computeSplines>(
                    sm_coordinates, gm_coordinates + kernelParams.pipelineAtomStart);
            __syncthreads();
            atomX = sm_coordinates[atomIndexLocal];
        }
        else
        {
            atomX = gm_coordinates[atomIndexGlobal];
        }
        calculate_splines<order, atomsPerBlock, atomsPerWarp, false, writeGlobal, numGrids, parallelExecutionWidth>(
                kernelParams, atomIndexOffset, atomX, atomCharge, sm_theta, &dtheta, sm_gridlineIndices);
        __builtin_amdgcn_wave_barrier();
    }
    else
    {
        /* Staging the data for spread
         * (the data is assumed to be in GPU global memory with proper layout already,
         * as in after running the spline kernel)
         */
        /* Spline data - only thetas (dthetas will only be needed in gather) */
        pme_gpu_stage_atom_data<float, atomsPerBlock, DIM * order, !computeSplines>(
                sm_theta, kernelParams.atoms.d_theta);
        /* Gridline indices */
        pme_gpu_stage_atom_data<int, atomsPerBlock, DIM, !computeSplines>(
                sm_gridlineIndices, kernelParams.atoms.d_gridlineIndices);

        __syncthreads();
    }

    /* Spreading */
    if constexpr (spreadCharges)
    {
        /* The boundary check is passed to spread_charges so that every lane
         * of the wave enters the function; the LDS path relies on wave-uniform
         * participation for its __shfl_xor-based bbox reduction. */
        const bool laneValid =
                (atomIndexGlobal < kernelParams.atoms.nAtoms)
                && (!kernelParams.usePipeline || (atomIndexGlobal < kernelParams.pipelineAtomEnd));
        spread_charges<order, wrapX, wrapY, 0, threadsPerAtom, parallelExecutionWidth>(
                kernelParams, &atomCharge, sm_gridlineIndices, sm_theta, sm_subGrid, laneValid);
    }
    if constexpr (numGrids == 2)
    {
        __syncthreads();
        if constexpr (c_useAtomDataPrefetch)
        {
            pme_gpu_stage_atom_data<float, atomsPerBlock, 1, !computeSplines>(
                    sm_coefficients, &kernelParams.atoms.d_coefficients[1][kernelParams.pipelineAtomStart]);
            __syncthreads();
            atomCharge = sm_coefficients[atomIndexLocal];
        }
        else
        {
            atomCharge = kernelParams.atoms.d_coefficients[1][atomIndexGlobal];
        }
        if constexpr (spreadCharges)
        {
            const bool laneValid =
                    (atomIndexGlobal < kernelParams.atoms.nAtoms)
                    && (!kernelParams.usePipeline || (atomIndexGlobal < kernelParams.pipelineAtomEnd));
            spread_charges<order, wrapX, wrapY, 1, threadsPerAtom, parallelExecutionWidth>(
                    kernelParams, &atomCharge, sm_gridlineIndices, sm_theta, sm_subGrid, laneValid);
        }
    }
}

//! Kernel instantiations
/* Disable the "explicit template instantiation 'PmeSplineAndSpreadKernel<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * It is only explicitly instantiated in this translation unit, so we should be safe.
 */
CLANG_DIAGNOSTIC_IGNORE("-Wweak-template-vtables")

#define INSTANTIATE_3(order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth)                     \
    template __global__ void                                                                                                                   \
    pmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth>( \
            PmeGpuKernelParams kernelParams);

#define INSTANTIATE_2(order, numGrids, threadsPerAtom, parallelExecutionWidth)                 \
    INSTANTIATE_3(order, true, true, numGrids, true, threadsPerAtom, parallelExecutionWidth);  \
    INSTANTIATE_3(order, true, false, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_3(order, false, true, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_3(order, true, true, numGrids, false, threadsPerAtom, parallelExecutionWidth);

#define INSTANTIATE(order, parallelExecutionWidth)                                 \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::OrderSquared, parallelExecutionWidth); \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::OrderSquared, parallelExecutionWidth);

INSTANTIATE(4, 32);
INSTANTIATE(4, 64);

CLANG_DIAGNOSTIC_RESET
