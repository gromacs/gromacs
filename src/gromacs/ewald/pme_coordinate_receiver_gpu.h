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
 * \brief Declaration of class which receives coordinates to GPU memory on PME task
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_PMECOORDINATERECEIVERGPU_H
#define GMX_PMECOORDINATERECEIVERGPU_H

#include <memory>
#include <tuple>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/vectypes.h"

class DeviceStream;
class DeviceContext;
class GpuEventSynchronizer;

struct PpRanks;

namespace gmx
{

template<typename>
class ArrayRef;

/*! \libinternal
 * \brief Manages receiving coordinates on PME-only ranks from their PP ranks.
 *
 * For multi-GPU runs, the PME GPU can receive coordinates from
 * multiple PP GPUs. Data from these distinct communications can
 * be handled separately in the PME spline/spread kernel, allowing
 * pipelining which overlaps computation and communication.
 *
 * Note that the PME rank always transfers coordinates from each PP
 * rank each step, even from empty domains. */
class PmeCoordinateReceiverGpu
{

public:
    /*! \brief Creates PME GPU coordinate receiver object
     *
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
     * Launch the receive of coordinates from a PP rank.
     *
     * With lib-MPI, a nonblocking GPU-aware receive of the coordinate data is posted.
     *
     * With thread-MPI, a nonblocking receive is posted for a pointer to the event
     * signalling the completion of the coordinate copy (recvbuf, atomOffset, and numBytes are ignored).
     * The coordinate copy itself is done by the PP rank directly into the buffer
     * whose address was sent by \c reinitCoordinateReceiver.
     *
     * \param[in] recvbuf      coordinates buffer in GPU memory
     * \param[in] atomOffset   starting element in buffer
     * \param[in] numBytes     number of bytes to transfer
     * \param[in] senderIndex  Index of PP rank within those involved in communication with this PME rank
     */
    void launchReceiveCoordinatesFromPp(DeviceBuffer<RVec> recvbuf, int atomOffset, int numBytes, int senderIndex);

    /*! \brief
     * Wait for a coordinate receive from any PP rank.
     *
     * With thread-MPI, this waits for the event pointer; the GPU copy
     * may still be in flight when this function returns.
     *
     * The returned sender index may correspond to a PP rank with an
     * empty domain, which transfers no particles but still takes part
     * in the communication.
     *
     * \returns  tuple with the index of the sending PP rank within the
     *           set of this rank's peers, and, with thread-MPI, the event
     *           that signals the completion of the D2D copy.
     */
    std::tuple<int, GpuEventSynchronizer*> waitForCoordinatesFromAnyPpRank();

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
     * Return number of PP ranks communicating with this PME rank
     */
    int ppCommNumRanks();

    /*! \brief Mark an event in the sender stream \p senderIndex
     * (which must be valid) and enqueue it into \p stream.
     */
    void insertAsDependencyIntoStream(int senderIndex, const DeviceStream& stream);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
