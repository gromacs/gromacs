/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *  SYCL implementation of coordinate layout conversion kernel and its wrapper
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/nbnxm_gpu_buffer_ops_internal.h"

#include "nbnxm_sycl_types.h"

using mode = sycl::access_mode;

namespace gmx
{

/*! \brief SYCL kernel for transforming position coordinates from rvec to nbnxm layout.
 *
 * Processes cells from multiple grids in a single kernel launch.
 * itemIdx.get(0) covers atoms within a cell, itemIdx.get(1) maps to a
 * global cell index across all grids. The grid membership of each cell
 * is determined using \p gridParams: columnsPrefix[] holds the prefix sum
 * of per-grid cell counts and binOffset[] the starting bin for each
 * grid.
 *
 * \param[out]    gm_xq                Coordinates buffer in nbnxm layout.
 * \param[in]     gm_x                 Coordinates buffer.
 * \param[in]     gm_atomIndex         Atom index mapping.
 * \param[in]     gm_numAtomsPerCell   Array of number of atoms per cell (all grids, stride numCellsMax).
 * \param[in]     gm_cellToBin         Array of bin indices per cell (all grids, stride numCellsMax).
 * \param[in]     numAtomsPerBin       Number of atoms per bin.
 * \param[in]     gridBegin            Index of first grid in gridset.
 * \param[in]     numCellsMax          Max cells per grid (stride for per-grid arrays).
 * \param[in]     gridParams           Per-grid parameters (prefix sums and bin offsets).
 */
static auto nbnxmKernelTransformXToXq(Float4* __restrict__ gm_xq,
                                      const Float3* __restrict__ gm_x,
                                      const int* __restrict__ gm_atomIndex,
                                      const int* __restrict__ gm_numAtomsPerCell,
                                      const int* __restrict__ gm_cellToBin,
                                      int                  numAtomsPerBin,
                                      int                  gridBegin,
                                      int                  numCellsMax,
                                      FusedXToXqGridParams gridParams)
{
    return [=](sycl::id<2> itemIdx)
    {
        const int globalCell = itemIdx.get(1);

        // Start assuming the cell belongs to grid 0, using constant [0] indices.
        int localCell = globalCell;
        int gridIdx   = gridBegin;
        int binOffset = gridParams.binOffset[0];

        // Last-write-wins scan over strictly increasing columnsPrefix[]: the highest matching g is
        // the correct grid. Full unrolling with constant indices is required to avoid local memory
        // spills; unused entries use INT_MAX sentinels so they never match.
#pragma unroll
        for (int g = 1; g < c_maxGridsPerKernelLaunch; g++)
        {
            if (globalCell >= gridParams.columnsPrefix[g])
            {
                localCell = globalCell - gridParams.columnsPrefix[g];
                gridIdx   = g + gridBegin;
                binOffset = gridParams.binOffset[g];
            }
        }

        // Access per-grid device arrays using numCellsMax stride
        const int cxy      = numCellsMax * gridIdx + localCell;
        const int numAtoms = gm_numAtomsPerCell[cxy];
        const int offset   = (binOffset + gm_cellToBin[cxy]) * numAtomsPerBin;

        const int threadIndex = itemIdx.get(0);

        // Perform layout conversion of each element.
        if (threadIndex < numAtoms)
        {
            const float  q              = gm_xq[threadIndex + offset][3];
            const Float3 xNew           = gm_x[gm_atomIndex[threadIndex + offset]];
            gm_xq[threadIndex + offset] = Float4(xNew[0], xNew[1], xNew[2], q);
        }
    };
}

// SYCL 1.2.1 requires providing a unique type for a kernel. Should not be needed for SYCL2020.
class NbnxmKernelTransformXToXqName;

void launchNbnxmKernelTransformXToXq(const FusedXToXqLaunchParams& launchParams,
                                     NbnxmGpu*                     nb,
                                     DeviceBuffer<Float3>          d_x,
                                     const DeviceStream&           deviceStream)
{
    const sycl::range<2> globalSize{ static_cast<size_t>(launchParams.maxNumAtomsPerCell),
                                     static_cast<size_t>(launchParams.totalNumCells) };
    sycl::queue          q = deviceStream.stream();

    auto kernelFunctionBuilder = nbnxmKernelTransformXToXq;
    syclSubmitWithoutCghOrEvent<NbnxmKernelTransformXToXqName>(q,
                                                               kernelFunctionBuilder,
                                                               globalSize,
                                                               nb->atdat->xq.get_pointer(),
                                                               d_x.get_pointer(),
                                                               nb->atomIndices.get_pointer(),
                                                               nb->numAtomsPerCell.get_pointer(),
                                                               nb->cellToBin.get_pointer(),
                                                               launchParams.numAtomsPerBin,
                                                               launchParams.gridBegin,
                                                               launchParams.numCellsMax,
                                                               launchParams.gridParams);
}

} // namespace gmx
