/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  \brief Implements helper routines for PME gather and spline routines.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"

#include "pme.cuh"
#include "pme_gpu_utils.h"
#include "pme_grid.h"

//! Controls if the atom and charge data is prefeched into shared memory or loaded per thread from global
static const bool c_useAtomDataPrefetch = true;

/*! \brief
 * General purpose function for loading atom-related data from global to shared memory.
 *
 * \tparam[in] T                 Data type (float/int/...)
 * \tparam[in] atomsPerBlock     Number of atoms processed by a block - should be accounted for in the size of the shared memory array.
 * \tparam[in] dataCountPerAtom  Number of data elements per single atom (e.g. DIM for an rvec coordinates array).
 * \param[in]  kernelParams      Input PME CUDA data in constant memory.
 * \param[out] sm_destination    Shared memory array for output.
 * \param[in]  gm_source         Global memory array for input.
 */
template<typename T, const int atomsPerBlock, const int dataCountPerAtom>
__device__ __forceinline__ void pme_gpu_stage_atom_data(const PmeGpuCudaKernelParams kernelParams,
                                                        T* __restrict__ sm_destination,
                                                        const T* __restrict__ gm_source)
{
    static_assert(c_usePadding,
                  "With padding disabled, index checking should be fixed to account for spline "
                  "theta/dtheta pr-warp alignment");
    const int blockIndex       = blockIdx.y * gridDim.x + blockIdx.x;
    const int threadLocalIndex = ((threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x) + threadIdx.x;
    const int localIndex       = threadLocalIndex;
    const int globalIndexBase = blockIndex * atomsPerBlock * dataCountPerAtom;
    const int globalIndex     = globalIndexBase + localIndex;
    const int globalCheck =
            pme_gpu_check_atom_data_index(globalIndex, kernelParams.atoms.nAtoms * dataCountPerAtom);
    if ((localIndex < atomsPerBlock * dataCountPerAtom) & globalCheck)
    {
        assert(isfinite(float(gm_source[globalIndex])));
        sm_destination[localIndex] = gm_source[globalIndex];
    }
}

/*! \brief
 * PME GPU spline parameter and gridline indices calculation.
 * This corresponds to the CPU functions calc_interpolation_idx() and make_bsplines().
 * First stage of the whole kernel.
 *
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] atomsPerBlock        Number of atoms processed by a block - should be accounted for
 *                                  in the sizes of the shared memory arrays.
 * \tparam[in] atomsPerWarp         Number of atoms processed by a warp
 * \tparam[in] writeSmDtheta        Bool controling if the theta derivative should be written to shared memory. Enables calculation of dtheta if set.
 * \tparam[in] writeGlobal          A boolean which tells if the theta values and gridlines should be written to global memory. Enables calculation of dtheta if set.
 * \param[in]  kernelParams         Input PME CUDA data in constant memory.
 * \param[in]  atomIndexOffset      Starting atom index for the execution block w.r.t. global memory.
 * \param[in]  atomX                Atom coordinate of atom processed by thread.
 * \param[in]  atomCharge           Atom charge/coefficient of atom processed by thread.
 * \param[out] sm_theta             Atom spline values in the shared memory.
 * \param[out] sm_dtheta            Derivative of atom spline values in shared memory.
 * \param[out] sm_gridlineIndices   Atom gridline indices in the shared memory.
 */

template<const int order, const int atomsPerBlock, const int atomsPerWarp, const bool writeSmDtheta, const bool writeGlobal>
__device__ __forceinline__ void calculate_splines(const PmeGpuCudaKernelParams kernelParams,
                                                  const int                    atomIndexOffset,
                                                  const float3                 atomX,
                                                  const float                  atomCharge,
                                                  float* __restrict__ sm_theta,
                                                  float* __restrict__ sm_dtheta,
                                                  int* __restrict__ sm_gridlineIndices)
{
    /* Global memory pointers for output */
    float* __restrict__ gm_theta         = kernelParams.atoms.d_theta;
    float* __restrict__ gm_dtheta        = kernelParams.atoms.d_dtheta;
    int* __restrict__ gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    /* Fractional coordinates */
    __shared__ float sm_fractCoords[atomsPerBlock * DIM];

    /* Thread index w.r.t. block */
    const int threadLocalId =
            (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / warp_size;
    /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
    const int atomWarpIndex = threadIdx.z % atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;

    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;
    /* Spline contribution index in one dimension */
    const int threadLocalIdXY = (threadIdx.y * blockDim.x) + threadIdx.x;
    const int orderIndex      = threadLocalIdXY / DIM;
    /* Dimension index */
    const int dimIndex = threadLocalIdXY % DIM;

    /* Multi-purpose index of rvec/ivec atom data */
    const int sharedMemoryIndex = atomIndexLocal * DIM + dimIndex;

    float splineData[order];

    const int localCheck = (dimIndex < DIM) && (orderIndex < 1);
    const int globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);

    /* we have 4 threads per atom, but can only use 3 here for the dimensions */
    if (localCheck && globalCheck)
    {
        /* Indices interpolation */

        if (orderIndex == 0)
        {
            int   tableIndex, tInt;
            float n, t;
            assert(atomIndexLocal < DIM * atomsPerBlock);
            /* Accessing fields in fshOffset/nXYZ/recipbox/... with dimIndex offset
             * puts them into local memory(!) instead of accessing the constant memory directly.
             * That's the reason for the switch, to unroll explicitly.
             * The commented parts correspond to the 0 components of the recipbox.
             */
            switch (dimIndex)
            {
                case XX:
                    tableIndex = kernelParams.grid.tablesOffsets[XX];
                    n          = kernelParams.grid.realGridSizeFP[XX];
                    t          = atomX.x * kernelParams.current.recipBox[dimIndex][XX]
                        + atomX.y * kernelParams.current.recipBox[dimIndex][YY]
                        + atomX.z * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;

                case YY:
                    tableIndex = kernelParams.grid.tablesOffsets[YY];
                    n          = kernelParams.grid.realGridSizeFP[YY];
                    t = /*atomX.x * kernelParams.current.recipBox[dimIndex][XX] + */ atomX.y
                                * kernelParams.current.recipBox[dimIndex][YY]
                        + atomX.z * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;

                case ZZ:
                    tableIndex = kernelParams.grid.tablesOffsets[ZZ];
                    n          = kernelParams.grid.realGridSizeFP[ZZ];
                    t          = /*atomX.x * kernelParams.current.recipBox[dimIndex][XX] + atomX.y * kernelParams.current.recipBox[dimIndex][YY] + */ atomX
                                .z
                        * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;
            }
            const float shift = c_pmeMaxUnitcellShift;
            /* Fractional coordinates along box vectors, adding a positive shift to ensure t is positive for triclinic boxes */
            t    = (t + shift) * n;
            tInt = (int)t;
            assert(sharedMemoryIndex < atomsPerBlock * DIM);
            sm_fractCoords[sharedMemoryIndex] = t - tInt;
            tableIndex += tInt;
            assert(tInt >= 0);
            assert(tInt < c_pmeNeighborUnitcellCount * n);

            // TODO have shared table for both parameters to share the fetch, as index is always same?
            // TODO compare texture/LDG performance
            sm_fractCoords[sharedMemoryIndex] +=
                    fetchFromParamLookupTable(kernelParams.grid.d_fractShiftsTable,
                                              kernelParams.fractShiftsTableTexture, tableIndex);
            sm_gridlineIndices[sharedMemoryIndex] =
                    fetchFromParamLookupTable(kernelParams.grid.d_gridlineIndicesTable,
                                              kernelParams.gridlineIndicesTableTexture, tableIndex);
            if (writeGlobal)
            {
                gm_gridlineIndices[atomIndexOffset * DIM + sharedMemoryIndex] =
                        sm_gridlineIndices[sharedMemoryIndex];
            }
        }

        /* B-spline calculation */

        const int chargeCheck = pme_gpu_check_atom_charge(atomCharge);
        if (chargeCheck)
        {
            float div;
            int o = orderIndex; // This is an index that is set once for PME_GPU_PARALLEL_SPLINE == 1

            const float dr = sm_fractCoords[sharedMemoryIndex];
            assert(isfinite(dr));

            /* dr is relative offset from lower cell limit */
            splineData[order - 1] = 0.0f;
            splineData[1]         = dr;
            splineData[0]         = 1.0f - dr;

#pragma unroll
            for (int k = 3; k < order; k++)
            {
                div               = 1.0f / (k - 1.0f);
                splineData[k - 1] = div * dr * splineData[k - 2];
#pragma unroll
                for (int l = 1; l < (k - 1); l++)
                {
                    splineData[k - l - 1] =
                            div * ((dr + l) * splineData[k - l - 2] + (k - l - dr) * splineData[k - l - 1]);
                }
                splineData[0] = div * (1.0f - dr) * splineData[0];
            }

            const int thetaIndexBase =
                    getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
            const int thetaGlobalOffsetBase = atomIndexOffset * DIM * order;
            /* only calculate dtheta if we are saving it to shared or global memory */
            if (writeSmDtheta || writeGlobal)
            {
                /* Differentiation and storing the spline derivatives (dtheta) */
#pragma unroll
                for (o = 0; o < order; o++)
                {
                    const int thetaIndex =
                            getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);

                    const float dtheta = ((o > 0) ? splineData[o - 1] : 0.0f) - splineData[o];
                    assert(isfinite(dtheta));
                    assert(thetaIndex < order * DIM * atomsPerBlock);
                    if (writeSmDtheta)
                    {
                        sm_dtheta[thetaIndex] = dtheta;
                    }
                    if (writeGlobal)
                    {
                        const int thetaGlobalIndex  = thetaGlobalOffsetBase + thetaIndex;
                        gm_dtheta[thetaGlobalIndex] = dtheta;
                    }
                }
            }

            div                   = 1.0f / (order - 1.0f);
            splineData[order - 1] = div * dr * splineData[order - 2];
#pragma unroll
            for (int k = 1; k < (order - 1); k++)
            {
                splineData[order - k - 1] = div
                                            * ((dr + k) * splineData[order - k - 2]
                                               + (order - k - dr) * splineData[order - k - 1]);
            }
            splineData[0] = div * (1.0f - dr) * splineData[0];

            /* Storing the spline values (theta) */
#pragma unroll
            for (o = 0; o < order; o++)
            {
                const int thetaIndex =
                        getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);
                assert(thetaIndex < order * DIM * atomsPerBlock);
                sm_theta[thetaIndex] = splineData[o];
                assert(isfinite(sm_theta[thetaIndex]));
                if (writeGlobal)
                {
                    const int thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;
                    gm_theta[thetaGlobalIndex] = splineData[o];
                }
            }
        }
    }
}
