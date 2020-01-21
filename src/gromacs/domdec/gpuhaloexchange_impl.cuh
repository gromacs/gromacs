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
/*! \internal \file
 *
 * \brief Declares CUDA implementation of GPU Halo Exchange.
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_GPUHALOEXCHANGE_IMPL_H
#define GMX_DOMDEC_GPUHALOEXCHANGE_IMPL_H

#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

/*! \brief switch for whether coordinates or force halo is being applied */
enum class HaloQuantity
{
    HaloCoordinates,
    HaloForces
};

/*! \internal \brief Class with interfaces and data for GPU Halo Exchange */
class GpuHaloExchange::Impl
{

public:
    /*! \brief Creates GPU Halo Exchange object.
     *
     * \param [inout] dd                       domdec structure
     * \param [in]    mpi_comm_mysim           communicator used for simulation
     * \param [in]    localStream              local NB CUDA stream
     * \param [in]    nonLocalStream           non-local NB CUDA stream
     */
    Impl(gmx_domdec_t* dd, MPI_Comm mpi_comm_mysim, void* localStream, void* nonLocalStream);
    ~Impl();

    /*! \brief
     * (Re-) Initialization for GPU halo exchange
     * \param [in] d_coordinatesBuffer  pointer to coordinates buffer in GPU memory
     * \param [in] d_forcesBuffer   pointer to forces buffer in GPU memory
     */
    void reinitHalo(float3* d_coordinatesBuffer, float3* d_forcesBuffer);


    /*! \brief
     * GPU halo exchange of coordinates buffer
     * \param [in] box  Coordinate box (from which shifts will be constructed)
     * \param [in] coordinatesReadyOnDeviceEvent event recorded when coordinates have been copied to device
     */
    void communicateHaloCoordinates(const matrix box, GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief  GPU halo exchange of force buffer
     * \param[in] accumulateForces  True if forces should accumulate, otherwise they are set
     */
    void communicateHaloForces(bool accumulateForces);

    /*! \brief Get the event synchronizer for the forces ready on device.
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getForcesReadyOnDeviceEvent();

private:
    /*! \brief Data transfer wrapper for GPU halo exchange
     * \param [inout] d_ptr      pointer to coordinates or force buffer in GPU memory
     * \param [in] haloQuantity  switch on whether X or F halo exchange is being performed
     * \param [in] coordinatesReadyOnDeviceEvent event recorded when coordinates have been copied to device
     */
    void communicateHaloData(float3*               d_ptr,
                             HaloQuantity          haloQuantity,
                             GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief Data transfer for GPU halo exchange using CUDA memcopies
     * \param [inout] sendPtr    address to send data from
     * \param [in] sendSize      number of atoms to be sent
     * \param [in] sendRank      rank to send data to
     * \param [inout] remotePtr  remote address to recv data
     * \param [in] recvRank      rank to recv data from
     */
    void communicateHaloDataWithCudaDirect(void* sendPtr, int sendSize, int sendRank, void* remotePtr, int recvRank);

    //! Domain decomposition object
    gmx_domdec_t* dd_ = nullptr;
    //! map of indices to be sent from this rank
    gmx::HostVector<int> h_indexMap_;
    //! device copy of index map
    int* d_indexMap_ = nullptr;
    //! number of elements in index map array
    int indexMapSize_ = -1;
    //! number of elements allocated in index map array
    int indexMapSizeAlloc_ = -1;
    //! device buffer for sending packed data
    float3* d_sendBuf_ = nullptr;
    //! number of atoms in sendbuf array
    int sendBufSize_ = -1;
    //! number of atoms allocated in sendbuf array
    int sendBufSizeAlloc_ = -1;
    //! device buffer for receiving packed data
    float3* d_recvBuf_ = nullptr;
    //! maximum size of packed buffer
    int maxPackedBufferSize_ = 0;
    //! number of atoms in recvbuf array
    int recvBufSize_ = -1;
    //! number of atoms allocated in recvbuf array
    int recvBufSizeAlloc_ = -1;
    //! rank to send data to for X
    int sendRankX_ = 0;
    //! rank to recv data from for X
    int recvRankX_ = 0;
    //! rank to send data to for F
    int sendRankF_ = 0;
    //! rank to recv data from for F
    int recvRankF_ = 0;
    //! send copy size from this rank for X
    int xSendSize_ = 0;
    //! recv copy size to this rank for X
    int xRecvSize_ = 0;
    //! send copy size from this rank for F
    int fSendSize_ = 0;
    //! recv copy size to this rank for F
    int fRecvSize_ = 0;
    //! number of home atoms - offset of local halo region
    int numHomeAtoms_ = 0;
    //! remote GPU coordinates buffer pointer for pushing data
    void* remoteXPtr_ = nullptr;
    //! remote GPU force buffer pointer for pushing data
    void* remoteFPtr_ = nullptr;
    //! Periodic Boundary Conditions for this rank
    bool usePBC_ = false;
    //! force shift buffer on device
    float3* d_fShift_ = nullptr;
    //! Event triggered when halo transfer has been launched with direct CUD memory copy
    GpuEventSynchronizer* haloDataTransferLaunched_ = nullptr;
    //! MPI communicator used for simulation
    MPI_Comm mpi_comm_mysim_;
    //! CUDA stream for local non-bonded calculations
    cudaStream_t localStream_ = nullptr;
    //! CUDA stream for non-local non-bonded calculations
    cudaStream_t nonLocalStream_ = nullptr;
    //! full coordinates buffer in GPU memory
    float3* d_x_ = nullptr;
    //! full forces buffer in GPU memory
    float3* d_f_ = nullptr;
    //! An event recorded once the exchanged forces are ready on the GPU
    GpuEventSynchronizer fReadyOnDevice_;
};

} // namespace gmx

#endif
