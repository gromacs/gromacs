/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#ifndef GMX_NBNXM_NBNXM_GPU_BUFFER_OPS_INTERNAL_H
#define GMX_NBNXM_NBNXM_GPU_BUFFER_OPS_INTERNAL_H

/*! \internal \file
 *  \brief
 *  Wrapper for the backend-specific coordinate layout conversion functionality
 *
 *  \ingroup module_nbnxm
 */
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"

class DeviceStream;

namespace gmx
{
struct NbnxmGpu;
class Grid;

/*! \brief Maximum number of grids supported in a single fused kernel launch */
static constexpr int c_maxGridsPerKernelLaunch = 8;

/*! \brief Per-grid parameters passed to the fused X-to-Xq kernel.
 *
 * This struct is passed by value as a kernel argument.
 */
struct FusedXToXqGridParams
{
    //! Prefix-sum start offsets: columnsPrefix[g] = total cells in grids before g.
    //! Indices [numGrids..c_maxGridsPerKernelLaunch-1] are INT_MAX sentinels (never matched).
    int columnsPrefix[c_maxGridsPerKernelLaunch];
    //! Bin offset for each grid.
    int binOffset[c_maxGridsPerKernelLaunch];
};

/*! \brief All parameters needed to launch the fused X-to-Xq kernel, computed
 *  once at pair-list setup and stored per interaction locality in NbnxmGpu.
 */
struct FusedXToXqLaunchParams
{
    //! Per-grid topology parameters passed to the kernel.
    FusedXToXqGridParams gridParams;
    //! Total number of cells across all grids (= grid dimension Y).
    int totalNumCells = 0;
    //! Maximum atoms-per-cell across all grids (determines grid dimension X).
    int maxNumAtomsPerCell = 0;
    //! Atoms per bin, identical for all grids.
    int numAtomsPerBin = 0;
    //! Index of the first grid in the gridset (0 = local, 1 = non-local).
    int gridBegin = 0;
    //! Maximum cells per grid; used as stride for per-grid device arrays.
    int numCellsMax = 0;
};

/*! \brief Launch the fused coordinate layout conversion kernel.
 *
 * Processes all grids in a single kernel launch by mapping cells from
 * all grids into a flat set of GPU blocks. Uses blockIdx.y for global
 * cell indexing across grids, and a prefix-sum lookup to determine
 * which grid each cell belongs to. All parameters are pre-computed at
 * pair-list setup and stored per interaction locality in
 * NbnxmGpu::xToXqLaunchParams.
 *
 * \param[in]     launchParams  Pre-computed kernel launch parameters.
 * \param[in,out] nb            Nbnxm main structure.
 * \param[in]     d_x           Source atom coordinates.
 * \param[in]     deviceStream  Device stream for kernel submission.
 */
void launchNbnxmKernelTransformXToXq(const FusedXToXqLaunchParams& launchParams,
                                     NbnxmGpu*                     nb,
                                     DeviceBuffer<Float3>          d_x,
                                     const DeviceStream&           deviceStream);

} // namespace gmx

#endif // GMX_NBNXM_NBNXM_GPU_BUFFER_OPS_INTERNAL_H
