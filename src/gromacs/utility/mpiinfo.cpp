/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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

#include "gromacs/utility/mpiinfo.h"

// need to include gmxapi.h here as mpi.h needs to be included before mpi-ext.h
#include "gromacs/utility/gmxmpi.h"

#if HAVE_CUDA_AWARE_MPI
#    include <mpi-ext.h>
#endif

namespace gmx
{

CudaAwareMpiStatus checkMpiCudaAwareSupport()
{
#if defined(MPIX_CUDA_AWARE_SUPPORT)
    // With OMPI version <=4.x, this function doesn't check if UCX PML is built with CUDA-support
    // or if CUDA is disabled at runtime.
    // Expect this function to work only if OMPI uses OB1 PML
    // This is a known issue (https://github.com/open-mpi/ompi/issues/7963) and fix for this is
    // expected soon (written March 2021)
    CudaAwareMpiStatus status = (MPIX_Query_cuda_support() == 1) ? CudaAwareMpiStatus::Supported
                                                                 : CudaAwareMpiStatus::NotSupported;
#else
    CudaAwareMpiStatus status = CudaAwareMpiStatus::NotKnown;
#endif

    return status;
}

} // namespace gmx
