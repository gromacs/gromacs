/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief Declarations for NVSHMEM initialization/finalize class.
 * gmxNvshmemHandle takes the MPI communicator and initializes the
 * NVSHMEM over all the ranks involved in the given MPI communicator.
 * This is a collective call for all the ranks in the given MPI comm.
 * After NVSHMEM initialization all NVSHMEM APIs can be safely used.
 *
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 *
 * \ingroup module_gpu_utils
 * \inlibraryapi
 */

#ifndef GMX_NVSHMEM_UTILS_H_
#define GMX_NVSHMEM_UTILS_H_

#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"

class gmxNvshmemHandle
{

private:
    MPI_Comm nvshmem_mpi_comm_;

public:
    gmxNvshmemHandle(const gmx::MDLogger& mdlog, MPI_Comm comm);

    ~gmxNvshmemHandle();

    //! Number of signal buffers types used for PP Halo exchange
    static const int numOfPpHaloExSyncBufs = 3;
    //! Size for the each of the 3 signal buffers used for PP Halo exchange
    int ppHaloExPerSyncBufSize_ = 0;
    //! Allocation size for the signal buffers used for PP Halo exchange
    int ppHaloExSyncBufSize_ = -1;
    //! Allocation capacity for the signal buffers used for PP Halo exchange
    int ppHaloExSyncBufCapacity_ = -1;

    DeviceBuffer<uint64_t> d_ppHaloExSyncBase_;

    /*! \brief Handles NVSHEM signal buffer initialization
     *
     * Allocates and initializes the signal buffers used in NVSHMEM enabled
     * PP Halo exchange.
     */
    void allocateAndInitSignalBufs(int                  totalDimsAndPulses,
                                   const DeviceContext& deviceContext_,
                                   const DeviceStream*  localStream);
};

#endif
