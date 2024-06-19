/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 *
 * \brief Definitions for NVSHMEM initialization/finalize class.
 * gmxNvshmemHandle takes the MPI communicator and initializes the
 * NVSHMEM over all the ranks involved in the given MPI communicator.
 * This is a collective call for all the ranks in the given MPI comm.
 * After NVSHMEM initialization all NVSHMEM APIs can be safely used.
 *
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 *
 * \ingroup module_gpu_utils
 * \inlibraryapi
 */

#include "nvshmem_utils.h"

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#if GMX_NVSHMEM
#    include <nvshmem.h>
#endif

gmxNvshmemHandle::gmxNvshmemHandle(MPI_Comm comm)
{
#if GMX_NVSHMEM
    // Duplicate the existing communicator for NVSHMEM usage as the communicator
    // should remain valid from nvshmem init to destruction.
    MPI_Comm_dup(comm, &nvshmem_mpi_comm_);
    nvshmemx_init_attr_t attr;
    attr.mpi_comm = (void*)&nvshmem_mpi_comm_;

    int ret = nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);

    GMX_RELEASE_ASSERT(ret == 0, "NVSHMEM init failed");
#else
    GMX_UNUSED_VALUE(nvshmem_mpi_comm_);
    GMX_UNUSED_VALUE(comm);
#endif
}

gmxNvshmemHandle::~gmxNvshmemHandle()
{
#if GMX_NVSHMEM
    // Call nvshmem_finalize before destroying the MPI Comm.
    nvshmem_finalize();
    MPI_Comm_free(&nvshmem_mpi_comm_);
#endif
}
