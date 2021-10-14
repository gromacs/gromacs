/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declaration of class which receives coordinates to GPU memory on PME task
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_PMECOORDINATERECEIVERGPU_H
#define GMX_PMECOORDINATERECEIVERGPU_H

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

class DeviceStream;
class DeviceContext;

struct PpRanks;

namespace gmx
{

template<typename>
class ArrayRef;

class PmeCoordinateReceiverGpu
{

public:
    /*! \brief Creates PME GPU coordinate receiver object
     *
     * For multi-GPU runs, the PME GPU can receive coordinates from
     * multiple PP GPUs. Data from these distinct communications can
     * be handled separately in the PME spline/spread kernel, allowing
     * pipelining which overlaps computation and communication. The
     * class methods are designed to called seperately for each remote
     * PP rank, and internally a different stream is used for each
     * remote PP rank to allow overlapping.
     *
     * \param[in] comm            Communicator used for simulation
     * \param[in] deviceContext   GPU context
     * \param[in] ppRanks         List of PP ranks
     */
    PmeCoordinateReceiverGpu(MPI_Comm comm, const DeviceContext& deviceContext, gmx::ArrayRef<PpRanks> ppRanks);
    ~PmeCoordinateReceiverGpu();

    /*! \brief
     * Re-initialize: set atom ranges and, for thread-MPI case,
     * send coordinates buffer address to PP rank
     * This is required after repartitioning since atom ranges and
     * buffer allocations may have changed.
     * \param[in] d_x   coordinates buffer in GPU memory
     */
    void reinitCoordinateReceiver(DeviceBuffer<RVec> d_x);


    /*! \brief
     * Receive coordinate synchronizer pointer from the PP ranks.
     * \param[in] ppRank  PP rank to receive the synchronizer from.
     */
    void receiveCoordinatesSynchronizerFromPpCudaDirect(int ppRank);

    /*! \brief
     * Used for lib MPI, receives co-ordinates from PP ranks
     * \param[in] recvbuf   coordinates buffer in GPU memory
     * \param[in] numAtoms  starting element in buffer
     * \param[in] numBytes  number of bytes to transfer
     * \param[in] ppRank    PP rank to send data
     */
    void launchReceiveCoordinatesFromPpCudaMpi(DeviceBuffer<RVec> recvbuf, int numAtoms, int numBytes, int ppRank);

    /*! \brief
     * For lib MPI, wait for coordinates from any PP rank
     * For thread MPI, enqueue PP co-ordinate transfer event received from PP
     * rank determined from pipeline stage into given stream
     * \param[in] pipelineStage  stage of pipeline corresponding to this transfer
     * \param[in] deviceStream   stream in which to enqueue the wait event.
     * \returns                  rank of sending PP task
     */
    int synchronizeOnCoordinatesFromPpRank(int pipelineStage, const DeviceStream& deviceStream);

    /*! \brief Perform above synchronizeOnCoordinatesFromPpRanks for all PP ranks,
     * enqueueing all events to a single stream
     * \param[in] deviceStream   stream in which to enqueue the wait events.
     */
    void synchronizeOnCoordinatesFromAllPpRanks(const DeviceStream& deviceStream);

    /*! \brief
     * Return pointer to stream associated with specific PP rank sender index
     * \param[in] senderIndex    Index of sender PP rank.
     */
    DeviceStream* ppCommStream(int senderIndex);

    /*! \brief
     * Returns range of atoms involved in communication associated with specific PP rank sender
     * index \param[in] senderIndex    Index of sender PP rank.
     */
    std::tuple<int, int> ppCommAtomRange(int senderIndex);

    /*! \brief
     * Return number of PP ranks involved in PME-PP communication
     */
    int ppCommNumSenderRanks();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
