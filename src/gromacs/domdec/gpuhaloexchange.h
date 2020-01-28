/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declaration of GPU halo exchange.
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_GPUHALOEXCHANGE_H
#define GMX_DOMDEC_GPUHALOEXCHANGE_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_domdec_t;
class GpuEventSynchronizer;

namespace gmx
{

/*! \libinternal
 * \brief Manages GPU Halo Exchange object */
class GpuHaloExchange
{

public:
    /*! \brief Creates GPU Halo Exchange object.
     *
     * Coordinate Halo exchange will be performed in \c
     * StreamNonLocal, and the \c communicateHaloCoordinates
     * method must be called before any subsequent operations that
     * access non-local parts of the coordinate buffer (such as
     * the non-local non-bonded kernels). It also must be called
     * after the local coordinates buffer operations (where the
     * coordinates are copied to the device and hence the \c
     * coordinatesReadyOnDeviceEvent is recorded). Force Halo exchange
     * will be performed in \c streamNonLocal (also potentally
     * with buffer clearing in \c streamLocal)and the \c
     * communicateHaloForces method must be called after the
     * non-local buffer operations, after the local force buffer
     * has been copied to the GPU (if CPU forces are present), and
     * before the local buffer operations. The force halo exchange
     * does not yet support virial steps.
     *
     * \param [inout] dd                       domdec structure
     * \param [in]    mpi_comm_mysim           communicator used for simulation
     * \param [in]    streamLocal              local NB CUDA stream.
     * \param [in]    streamNonLocal           non-local NB CUDA stream.
     */
    GpuHaloExchange(gmx_domdec_t* dd, MPI_Comm mpi_comm_mysim, void* streamLocal, void* streamNonLocal);
    ~GpuHaloExchange();

    /*! \brief
     *
     * Initialization for GPU halo exchange of coordinates buffer
     * \param [in] d_coordinateBuffer   pointer to coordinates buffer in GPU memory
     * \param [in] d_forcesBuffer   pointer to coordinates buffer in GPU memory
     */
    void reinitHalo(DeviceBuffer<float> d_coordinateBuffer, DeviceBuffer<float> d_forcesBuffer);


    /*! \brief GPU halo exchange of coordinates buffer.
     *
     * Must be called after local setCoordinates (which records an
     * event when the coordinate data has been copied to the
     * device).
     * \param [in] box  Coordinate box (from which shifts will be constructed)
     * \param [in] coordinatesReadyOnDeviceEvent event recorded when coordinates have been copied to device
     */
    void communicateHaloCoordinates(const matrix box, GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief GPU halo exchange of force buffer.
     * \param[in] accumulateForces  True if forces should accumulate, otherwise they are set
     */
    void communicateHaloForces(bool accumulateForces);

    /*! \brief Get the event synchronizer for the forces ready on device.
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getForcesReadyOnDeviceEvent();

private:
    class Impl;
    gmx::PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
