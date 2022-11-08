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

namespace Nbnxm
{

/*! \brief SYCL kernel for transforming position coordinates from rvec to nbnxm layout.
 *
 * \param         cgh                 SYCL's command group handler.
 * \param[out]    a_xq                Coordinates buffer in nbnxm layout.
 * \param[in]     a_x                 Coordinates buffer.
 * \param[in]     a_atomIndex         Atom index mapping.
 * \param[in]     a_numAtoms          Array of number of atoms.
 * \param[in]     a_cellIndex         Array of cell indices.
 * \param[in]     cellOffset          First cell.
 * \param[in]     numAtomsPerCell     Number of atoms per cell.
 * \param[in]     columnsOffset       Index if the first column in the cell.
 */
static auto nbnxmKernelTransformXToXq(sycl::handler&                           cgh,
                                      DeviceAccessor<Float4, mode::read_write> a_xq,
                                      DeviceAccessor<Float3, mode::read>       a_x,
                                      DeviceAccessor<int, mode::read>          a_atomIndex,
                                      DeviceAccessor<int, mode::read>          a_numAtoms,
                                      DeviceAccessor<int, mode::read>          a_cellIndex,
                                      int                                      cellOffset,
                                      int                                      numAtomsPerCell,
                                      int                                      columnsOffset)
{
    a_xq.bind(cgh);
    a_x.bind(cgh);
    a_atomIndex.bind(cgh);
    a_numAtoms.bind(cgh);
    a_cellIndex.bind(cgh);

    return [=](sycl::id<2> itemIdx) {
        // Map cell-level parallelism to y component of block index.
        const int cxy = itemIdx.get(1) + columnsOffset;

        const int numAtoms = a_numAtoms[cxy];
        const int offset   = (cellOffset + a_cellIndex[cxy]) * numAtomsPerCell;

        const int threadIndex = itemIdx.get(0);

        // Perform layout conversion of each element.
        if (threadIndex < numAtoms)
        {
            const float  q             = a_xq[threadIndex + offset][3];
            const Float3 xNew          = a_x[a_atomIndex[threadIndex + offset]];
            a_xq[threadIndex + offset] = Float4(xNew[0], xNew[1], xNew[2], q);
        }
    };
}

// SYCL 1.2.1 requires providing a unique type for a kernel. Should not be needed for SYCL2020.
class NbnxmKernelTransformXToXqName;

void launchNbnxmKernelTransformXToXq(const Grid&          grid,
                                     NbnxmGpu*            nb,
                                     DeviceBuffer<Float3> d_x,
                                     const DeviceStream&  deviceStream,
                                     unsigned int         numColumnsMax,
                                     int                  gridId)
{
    const unsigned int numColumns  = grid.numColumns();
    const unsigned int numAtomsMax = grid.numCellsColumnMax() * grid.numAtomsPerCell();
    GMX_ASSERT(numColumns <= numColumnsMax, "Grid has more columns than allowed");

    const sycl::range<2> globalSize{ numAtomsMax, numColumns };
    sycl::queue          q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = nbnxmKernelTransformXToXq(cgh,
                                                nb->atdat->xq,
                                                d_x,
                                                nb->atomIndices,
                                                nb->cxy_na,
                                                nb->cxy_ind,
                                                grid.cellOffset(),
                                                grid.numAtomsPerCell(),
                                                numColumnsMax * gridId);
        cgh.parallel_for<NbnxmKernelTransformXToXqName>(globalSize, kernel);
    });
}

} // namespace Nbnxm
