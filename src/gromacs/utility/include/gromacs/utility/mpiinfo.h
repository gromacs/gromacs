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
#ifndef GMX_UTILITY_MPI_INFO_H
#define GMX_UTILITY_MPI_INFO_H

#include <string_view>

namespace gmx
{
/*! \brief Enum describing GPU-aware support in underlying MPI library.
 *
 * Ordinal, so that the lowest value can represent the minimal level of
 * support found across a set of devices, perhaps across nodes or ranks.
 */
enum class GpuAwareMpiStatus : int
{
    NotSupported = 0, //!< GPU-aware support NOT available or not known.
    Forced,           //!< GPU-aware support forced using env variable
    Supported,        //!< GPU-aware support available.
};

//! Return the string obtained from the MPI library via MPI_Get_library_version
std::string_view mpiLibraryVersionString();

//! Return whether GROMACS is linked against an MPI library describing itself as Intel MPI
bool usingIntelMpi();

/*! \brief
 * Wrapper on top of \c MPIX_Query_cuda_support()
 * For MPI implementations which don't support this function, it returns \c NotSupported.
 * Even when an MPI implementation does support this function, MPI library might not be
 * robust enough to detect CUDA-aware support at runtime correctly e.g. when UCX PML is used
 * or CUDA is disabled at runtime
 *
 * \returns     CUDA-aware status in MPI implementation */
GpuAwareMpiStatus checkMpiCudaAwareSupport();

/*! \brief
 * Wrapper on top of \c MPIX_Query_hip_support() or \c MPIX_Query_rocm_support().
 * For MPI implementations which don't support this function, it returns \c NotSupported.
 *
 * Currently, this function is only supported by MPICH and OpenMPI 5.0-rc, and is not very reliable.
 *
 * \returns     HIP-aware status in MPI implementation */
GpuAwareMpiStatus checkMpiHipAwareSupport();

/*! \brief
 * Wrapper on top of \c MPIX_Query_ze_support() (for MPICH) or custom
 * logic (for IntelMPI).
 *
 * For other MPI implementations which perhaps don't support the above
 * function, it returns NotSupported.
 *
 * \returns     LevelZero-aware status in MPI implementation */
GpuAwareMpiStatus checkMpiZEAwareSupport();


} // namespace gmx

#endif
