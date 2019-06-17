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
 * \brief Declares CUDA implementation class for domdec
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_domdec
 */
#ifndef GMX_DD_GPU_IMPL_H
#define GMX_DD_GPU_IMPL_H

#include "gromacs/domdec/domdec_cuda.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/gpu_utils/hostallocator.h"

namespace gmx
{

/*! \brief switch for whether position or force halo is being applied */
enum class HaloXorF
{
    HaloX, HaloF
};

/*! \internal \brief Class with interfaces and data for CUDA version of LINCS. */
class DomdecCuda::Impl
{

    public:
        /*! \brief Creates Domdec GPU object.
         *
         * \param [inout] dd      domdec structure
         */
        Impl(gmx_domdec_t gmx_unused *dd);
        ~Impl();

        /*! \brief
         * Initialization for GPU halo exchange of position buffer
         * \param [inout] dd      domdec structure, to be populated with index and param info for exchange
         * \param [in] box        box matrix required for shift information
         * \param [in] d_x_ptr   pointer to position buffer in GPU memory
         * \param [in] stream_ptr   pointer to CUDA stream to be used for init operations
         */
        void initHaloExchange(gmx_domdec_t gmx_unused *dd, matrix gmx_unused box,
                              rvec gmx_unused *d_x_ptr, void gmx_unused *stream_ptr);


        /*! \brief
         * GPU halo exchange of position buffer
         * \param [inout] d_x_ptr   pointer to position buffer in GPU memory
         * \param [in] streamNonLocal       CUDA stream to be used for buffer packing operation
         */
        void applyXHaloExchange(rvec gmx_unused *d_x_ptr, void gmx_unused *streamNonLocal);

        /*! \brief
         * GPU halo exchange of force buffer
         * \param [inout] d_f_ptr   pointer to force buffer in GPU memory
         * \param [inout] fshift    force shift array
         * \param [in] stream       CUDA stream to be used for buffer packing operation
         */
        void applyFHaloExchange(rvec gmx_unused *d_f_ptr, rvec gmx_unused *fshift, void gmx_unused *stream);


    private:

        /*! \brief Data transfer wrapper for GPU halo exchange
         * \param [inout] d_ptr    pointer to position or force buffer in GPU memory
         * \param [in] streamNonLocal_ptr  pointer to CUDA stream to be used for buffer packing operation
         * \param [in] haloXorF    switch on whether X or F halo exchange is being performed
         */
        void haloDataTransfer(rvec gmx_unused *d_ptr, void gmx_unused *streamNonLocal_ptr, HaloXorF haloXorF);

        /*! \brief Data transfer for GPU halo exchange using CUDA-aware MPI
         * \param [inout] sendPtr    address to send data from
         * \param [in] sendSize      number of atoms to be sent
         * \param [in] sendRank      rank to send data to
         * \param [inout] recvPtr    address to recv data in
         * \param [in] recvSize      number of atoms to be received
         * \param [in] recvRank      rank to recv data from
         * \param [in] streamNonLocal_ptr    pointer to CUDA stream to be used for buffer packing operation
         */
        void haloDataTransferCudaMPI(void *sendPtr, int sendSize, int sendRank, void *recvPtr, int recvSize, int recvRank, void *streamNonLocal_ptr);

        /*! \brief Data transfer for GPU halo exchange using CUDA memcopies
         * \param [inout] sendPtr    address to send data from
         * \param [in] sendSize      number of atoms to be sent
         * \param [in] sendRank      rank to send data to
         * \param [inout] remotePtr    remote address to recv data
         * \param [in] recvRank      rank to recv data from
         * \param [in] streamNonLocal_ptr    Pointer to CUDA stream to be used for buffer packing operation
         */
        void haloDataTransferCudaDirect(void *sendPtr, int sendSize, int sendRank, void* remotePtr, int recvRank, void *streamNonLocal_ptr);


        /*! \brief Enable peer access between GPUs where supported */
        void setupCudaDevicePeerAccess();

        //! map of indices to be sent from this rank
        gmx::HostVector<int>       h_indexMap_;
        //! device copy of index map
        int                       *d_indexMap_;
        //! number of elements in index map array
        int                        indexMap_size_         = 0;
        //! number of elements allocated in index map array
        int                        indexMap_size_alloc_   = 0;
        //! host buffer for sending packed data
        gmx::HostVector<gmx::RVec> h_sendBuf_;
        //! device buffer for sending packed data
        rvec                      *d_sendBuf_;
        //! number of atoms in sendbuf array
        int                        sendBuf_size_          = 0;
        //! number of atoms allocated in sendbuf array
        int                        sendBuf_size_alloc_    = 0;
        //! host buffer for receiving packed data
        gmx::HostVector<gmx::RVec> h_recvBuf_;
        //! device buffer for receiving packed data
        rvec                      *d_recvBuf_;
        //! number of atoms in recvbuf array
        int                        recvBuf_size_           = 0;
        //! number of atoms allocated in recvbuf array
        int                        recvBuf_size_alloc_     = 0;
        //! rank to send data to for X
        int                        sendRankX_               = 0;
        //! rank to recv data from for X
        int                        recvRankX_               = 0;
        //! rank to send data to for F
        int                        sendRankF_               = 0;
        //! rank to recv data from for F
        int                        recvRankF_               = 0;
        //! send copy size from this rank for X
        int                        xSendSize_              = 0;
        //! recv copy size to this rank for X
        int                        xRecvSize_              = 0;
        //! send copy size from this rank for F
        int                        fSendSize_              = 0;
        //! recv copy size to this rank for F
        int                        fRecvSize_              = 0;
        //! offset of local halo region
        int                        localOffset_            = 0;
        //! remote GPU position buffer pointer for pushing data
        void                      *remoteXPtr_            = 0;
        //! remote GPU force buffer pointer for pushing data
        void                      *remoteFPtr_            = 0;
        //! Periodic Boundary Conditions for this rank
        bool                       usePBC_                 = false;
        //! Whether shift force update is required
        bool                       bShiftForcesNeedPbc_    = false;
        //! force shift buffer on device
        float3*                    d_fShift_;
        //! X shift values
        float3                     xShift_;
        //! Event triggered when halo transfer has been launched with direct CUD memory copy
        GpuEventSynchronizer      *haloDataTransferLaunched_ = nullptr;
};

} // namespace gmx

#endif
