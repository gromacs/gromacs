/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Declaration of GPU PME-PP Communication.
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_PME_PP_COMM_GPU_H
#define GMX_PME_PP_COMM_GPU_H

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

class DeviceContext;
class DeviceStream;
class GpuEventSynchronizer;

namespace gmx
{

class DeviceStreamManager;

/*! \libinternal

 * \brief Manages communication related to GPU buffers between this
 * PME rank and its PP rank. */
class PmePpCommGpu
{

public:
    /*! \brief Creates PME-PP GPU communication object
     * \param[in] comm              Communicator used for simulation
     * \param[in] pmeRank           Rank of PME task
     * \param[in] pmeCpuForceBuffer Buffer for PME force in CPU memory
     * \param[in] deviceContext     GPU context.
     * \param[in] deviceStream      GPU stream.
     */
    PmePpCommGpu(MPI_Comm                    comm,
                 int                         pmeRank,
                 gmx::HostVector<gmx::RVec>* pmeCpuForceBuffer,
                 const DeviceContext&        deviceContext,
                 const DeviceStream&         deviceStream);
    ~PmePpCommGpu();

    /*! \brief Perform steps required when buffer size changes
     * \param[in]  size   Number of elements in buffer
     */
    void reinit(int size);

    /*! \brief
     * Pull data from PME GPU directly using CUDA Memory copy.
     * \param[out] recvPtr  Buffer to receive PME force data
     * \param[in] recvSize Number of elements to receive
     * \param[in] recvPmeForceToGpu Whether receive is to GPU, otherwise CPU
     */
    void receiveForceFromPme(RVec* recvPtr, int recvSize, bool recvPmeForceToGpu);

    /*! \brief Push coordinates buffer directly to GPU memory on PME task
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     * \param[in] coordinatesReadyOnDeviceEvent Event recorded when coordinates are available on device
     */
    void sendCoordinatesToPmeFromGpu(DeviceBuffer<RVec>    sendPtr,
                                     int                   sendSize,
                                     GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief Push coordinates buffer from host memory directly to GPU memory on PME task
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     */
    void sendCoordinatesToPmeFromCpu(RVec* sendPtr, int sendSize);

    /*! \brief
     * Return pointer to buffer used for staging PME force on GPU
     */
    DeviceBuffer<gmx::RVec> getGpuForceStagingPtr();

    /*! \brief
     * Return pointer to event recorded when forces are ready
     */
    GpuEventSynchronizer* getForcesReadySynchronizer();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
