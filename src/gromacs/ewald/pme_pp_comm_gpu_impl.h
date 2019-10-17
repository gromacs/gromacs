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
 * \brief Declares CUDA implementation class for PME-PP communications
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#ifndef GMX_PME_PP_COMM_GPU_IMPL_H
#define GMX_PME_PP_COMM_GPU_IMPL_H

#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

/*! \internal \brief Class with interfaces and data for CUDA version of PME-PP Communication */
class PmePpCommGpu::Impl
{

public:
    /*! \brief Creates PME-PP GPU communication object.
     * \param[in] comm            Communicator used for simulation
     * \param[in] pmeRank         Rank of PME task
     */
    Impl(MPI_Comm comm, int pmeRank);
    ~Impl();

    /*! \brief Perform steps required when buffer size changes
     * \param[in]  size   Number of elements in buffer
     */
    void reinit(int size);

    /*! \brief Pull force buffer directly from GPU memory on PME
     * rank to either GPU or CPU memory on PP task using CUDA
     * Memory copy.
     *
     * recvPtr should be in GPU or CPU memory if recvPmeForceToGpu
     * is true or false, respectively. If receiving to GPU, this
     * method should be called before the local GPU buffer
     * operations. If receiving to CPU it should be called
     * before forces are reduced with the other force
     * contributions on the CPU. It will automatically wait for
     * remote PME force data to be ready.
     *
     * \param[out] recvPtr CPU buffer to receive PME force data
     * \param[in] recvSize Number of elements to receive
     * \param[in] receivePmeForceToGpu Whether receive is to GPU, otherwise CPU
     */
    void receiveForceFromPmeCudaDirect(void* recvPtr, int recvSize, bool receivePmeForceToGpu);


    /*! \brief Push coordinates buffer directly to GPU memory on PME
     * task, from either GPU or CPU memory on PP task using CUDA
     * Memory copy. sendPtr should be in GPU or CPU memory if
     * sendPmeCoordinatesFromGpu is true or false respectively. If
     * sending from GPU, this method should be called after the
     * local GPU coordinate buffer operations. The remote PME task will
     * automatically wait for data to be copied before commencing PME force calculations.
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     * \param[in] sendPmeCoordinatesFromGpu Whether send is from GPU, otherwise CPU
     * \param[in] coordinatesReadyOnDeviceEvent Event recorded when coordinates are available on device
     */
    void sendCoordinatesToPmeCudaDirect(void*                 sendPtr,
                                        int                   sendSize,
                                        bool                  sendPmeCoordinatesFromGpu,
                                        GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief
     * Return pointer to buffer used for staging PME force on GPU
     */
    void* getGpuForceStagingPtr();

    /*! \brief
     * Return pointer to event recorded when forces are ready
     */
    void* getForcesReadySynchronizer();

private:
    //! CUDA stream used for the communication operations in this class
    cudaStream_t pmePpCommStream_ = nullptr;
    //! Remote location of PME coordinate data buffer
    void* remotePmeXBuffer_ = nullptr;
    //! Remote location of PME force data buffer
    void* remotePmeFBuffer_ = nullptr;
    //! communicator for simulation
    MPI_Comm comm_;
    //! Rank of PME task
    int pmeRank_ = -1;
    //! Buffer for staging PME force on GPU
    rvec* d_pmeForces_ = nullptr;
    //! number of atoms in PME force staging array
    int d_pmeForcesSize_ = -1;
    //! number of atoms allocated in recvbuf array
    int d_pmeForcesSizeAlloc_ = -1;
    //! Event recorded when PME forces are ready on PME task
    GpuEventSynchronizer forcesReadySynchronizer_;
    //! Event recorded when coordinates have been transferred to PME task
    GpuEventSynchronizer pmeCoordinatesSynchronizer_;
};

} // namespace gmx

#endif
