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

        /*! \brief Receive address of GPU force buffer on remote PME task, to allow this rank to
         * pull data directly from GPU memory space using a CUDA memory copy.
         */
        void receiveForceBufferAddress();

        /*! \brief Pull force buffer directly from GPU memory on PME
         * task to CPU memory on PP task using CUDA Memory copy.
         *
         * recvPtr should be in CPU memory. This method should be
         * called forces are reduced with the other force
         * contributions. It will automatically wait for remote PME
         * force data to be ready.
         * \param[out] recvPtr CPU buffer to receive PME force data
         * \param[in] recvSize Number of elements to receive
         */
        void receiveForceFromPmeCudaDirect(void *recvPtr, int recvSize);

    private:
        //! CUDA stream used for the communication operations in this class
        cudaStream_t           pmePpCommStream_ = nullptr;
        //! Remote location of PME force data buffer
        void                  *remotePmeFBuffer_ = nullptr;
        //! communicator for simulation
        MPI_Comm               comm_;
        //! Rank of PME task
        int                    pmeRank_ = -1;

};

} // namespace gmx

#endif
