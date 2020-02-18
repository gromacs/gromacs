/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

class DeviceStream;
struct PpRanks;

namespace gmx
{

template<typename>
class ArrayRef;

class PmeCoordinateReceiverGpu
{

public:
    /*! \brief Creates PME GPU coordinate receiver object
     * \param[in] pmeStream       CUDA stream used for PME computations
     * \param[in] comm            Communicator used for simulation
     * \param[in] ppRanks         List of PP ranks
     */
    PmeCoordinateReceiverGpu(const DeviceStream& pmeStream, MPI_Comm comm, gmx::ArrayRef<PpRanks> ppRanks);
    ~PmeCoordinateReceiverGpu();

    /*! \brief
     * send coordinates buffer address to PP rank
     * \param[in] d_x   coordinates buffer in GPU memory
     */
    void sendCoordinateBufferAddressToPpRanks(DeviceBuffer<RVec> d_x);


    /*! \brief
     * launch receive of coordinate data from PP rank
     * \param[in] ppRank  PP rank to send data
     */
    void launchReceiveCoordinatesFromPpCudaDirect(int ppRank);

    /*! \brief
     * enqueue wait for coordinate data from PP ranks
     */
    void enqueueWaitReceiveCoordinatesFromPpCudaDirect();

private:
    class Impl;
    gmx::PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
