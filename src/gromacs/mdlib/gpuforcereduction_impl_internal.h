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
/*! \internal \file
 *
 * \brief Declares vendor-specific function to launch force reduction kernel
 *
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_mdlib
 */

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"

class DeviceStream;

namespace gmx
{

/*! \brief Backend-specific function to launch GPU Force Reduction kernel.
 *
 * In pseudocode:
 *
 * \code{.cpp}
 * for (int i = 0; i < numAtoms; i++) {
 *     totalForce = d_nbnxmForceToAdd[d_cell[i]]
 *     if (accumulate)
 *         totalForce += d_baseForce[atomStart + i]
 *     if (addRvecForce)
 *         totalForce += d_rvecForceToAdd[atomStart + i]
 *     d_baseForce[atomStart + i] = totalForce[i]
 * }
 * \endcode
 *
 * \param numAtoms Number of atoms subject to reduction.
 * \param atomStart First atom index (for \p d_rvecForceToAdd and \p d_baseForce).
 * \param addRvecForce When \c false, \p d_rvecForceToAdd is ignored.
 * \param accumulate When \c false, the previous values of \p d_baseForce are discarded.
 * \param d_nbnxmForceToAdd Buffer containing Nbnxm forces in Nbnxm layout.
 * \param d_rvecForceToAdd Optional buffer containing arbitrary forces in linear layout.
 * \param d_baseForce Destination buffer for forces in linear layout.
 * \param d_cell Atom index to Nbnxm cell index.
 * \param deviceStream Device stream for kernel submission.
 * \param d_forcesReadyNvshmemFlags  NVSHMEM signals from PME to PP force transfer.
 * \param forcesReadyNvshmemFlagsCounter Tracks NVSHMEM signal from PME to PP force transfer.
 */
void launchForceReductionKernel(int                        numAtoms,
                                int                        atomStart,
                                bool                       addRvecForce,
                                bool                       accumulate,
                                const DeviceBuffer<Float3> d_nbnxmForceToAdd,
                                const DeviceBuffer<Float3> d_rvecForceToAdd,
                                DeviceBuffer<Float3>       d_baseForce,
                                DeviceBuffer<int>          d_cell,
                                const DeviceStream&        deviceStream,
                                DeviceBuffer<uint64_t>     d_forcesReadyNvshmemFlags,
                                const uint64_t             forcesReadyNvshmemFlagsCounter);

} // namespace gmx
