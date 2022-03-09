/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#include "gromacs/utility/mpiinfo.h"

#include <cstdlib>

// need to include gmxapi.h here as mpi.h needs to be included before mpi-ext.h
#include "gromacs/utility/gmxmpi.h"

#if HAVE_MPI_EXT
#    include <mpi-ext.h>
#endif

namespace gmx
{

GpuAwareMpiStatus checkMpiCudaAwareSupport()
{
#if MPI_SUPPORTS_CUDA_AWARE_DETECTION
    // With OMPI version <=4.x, this function doesn't check if UCX PML is built with CUDA-support
    // or if CUDA is disabled at runtime.
    // Expect this function to work only if OMPI uses OB1 PML
    // This is a known issue (https://github.com/open-mpi/ompi/issues/7963) and fix for this is
    // expected soon (written March 2021)
    GpuAwareMpiStatus status = (MPIX_Query_cuda_support() == 1) ? GpuAwareMpiStatus::Supported
                                                                : GpuAwareMpiStatus::NotSupported;
#else
    GpuAwareMpiStatus status = GpuAwareMpiStatus::NotKnown;
#endif

    if (status != GpuAwareMpiStatus::Supported && getenv("GMX_FORCE_GPU_AWARE_MPI") != nullptr)
    {
        status = GpuAwareMpiStatus::Forced;
    }

    return status;
}

} // namespace gmx
