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
 * \brief Declaration of class which sends PME Force from GPU memory to PP task
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_PMEFORCESENDERGPU_H
#define GMX_PMEFORCESENDERGPU_H

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

class GpuEventSynchronizer;
class DeviceContext;
class DeviceStream;

/*! \libinternal
 * \brief Contains information about the PP ranks that partner this PME rank. */
struct PpRanks
{
    //! The MPI rank ID of this partner PP rank.
    int rankId = -1;
    //! The number of atoms to communicate with this partner PP rank.
    int numAtoms = -1;
};

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \libinternal
 * \brief Manages sending forces from PME-only ranks to their PP ranks. */
class PmeForceSenderGpu
{

public:
    /*! \brief Creates PME GPU Force sender object
     * \param[in] pmeForcesReady  Event synchronizer marked when PME forces are ready on the GPU
     * \param[in] comm            Communicator used for simulation
     * \param[in] deviceContext   GPU context
     * \param[in] ppRanks         List of PP ranks
     */
    PmeForceSenderGpu(GpuEventSynchronizer*  pmeForcesReady,
                      MPI_Comm               comm,
                      const DeviceContext&   deviceContext,
                      gmx::ArrayRef<PpRanks> ppRanks);
    ~PmeForceSenderGpu();

    /*! \brief
     * Sets location of force to be sent to each PP rank
     * \param[in] d_f   force buffer in GPU memory
     */
    void setForceSendBuffer(DeviceBuffer<RVec> d_f);

    /*! \brief
     * Send force to PP rank (used with Thread-MPI)
     * \param[in] ppRank                   PP rank to receive data
     * \param[in] numAtoms                 number of atoms to send
     * \param[in] sendForcesDirectToPpGpu  whether forces are transferred direct to remote GPU memory
     */
    void sendFToPpPeerToPeer(int ppRank, int numAtoms, bool sendForcesDirectToPpGpu);

    /*! \brief
     * Send force to PP rank (used with Lib-MPI)
     * \param[in] sendbuf  force buffer in GPU memory
     * \param[in] offset   starting element in buffer
     * \param[in] numBytes number of bytes to transfer
     * \param[in] ppRank   PP rank to receive data
     * \param[in] request  MPI request to track asynchronous MPI call status
     */
    void sendFToPpGpuAwareMpi(DeviceBuffer<RVec> sendbuf, int offset, int numBytes, int ppRank, MPI_Request* request);

    void waitForEvents();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
