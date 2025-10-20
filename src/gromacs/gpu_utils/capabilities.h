/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_CAPABILITIES_H
#define GMX_GPU_UTILS_CAPABILITIES_H

/*! \libinternal \file
 *  \brief Declares the GPU capabilities on a configuration specific basis.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */

#include "config.h"

namespace gmx
{

//! Collection of GPU capabilities for this configuration
struct GpuConfigurationCapabilities
{
    //! Whether this configuration supports running NBNXM kernels on the device
    static constexpr bool Nonbonded = GMX_GPU;
    //! Whether this configuration supports running buffer operation kernels on the device
    static constexpr bool BufferOps = GMX_GPU && !GMX_GPU_OPENCL;
    //! Whether this configuration supports running PME kernels on the device
    static constexpr bool Pme = GMX_GPU && !GMX_GPU_HIP;
    //! Whether this configuration supports querying device stream for completion of all submitted tasks
    static constexpr bool StreamQuery = GMX_GPU_CUDA || GMX_GPU_HIP;
    //! Whether this configuration supports using the param lookup table
    static constexpr bool PmeParamLookupTable = GMX_GPU_CUDA || GMX_GPU_HIP;
    //! Whether this configuration supports PME pipelining
    static constexpr bool PmePipelining = GMX_GPU && !GMX_GPU_OPENCL;
    //! Whether this configuration supports PME ThreadsPerAtomOrder
    static constexpr bool PmeSupportsThreadsPerAtomOrder = GMX_GPU && !GMX_GPU_OPENCL;
    //! Whether this configuration supports PME max gridsize setting
    static constexpr bool PmeDynamicMaxGridSize = GMX_GPU_CUDA || GMX_GPU_HIP;
    //! Whether this configuration supports PME solve kernels with less than 4 warps
    static constexpr bool PmeSolveNeedsAtLeastFourWarps = GMX_GPU && !GMX_GPU_OPENCL;
    //! Whether this configuration supports running FFT kernels on the device
    static constexpr bool Fft =
            GMX_GPU
            && ((GMX_GPU_FFT_MKL != 0) || (GMX_GPU_FFT_ROCFFT != 0) || (GMX_GPU_FFT_HIPFFT != 0)
                || (GMX_GPU_FFT_BBFFT != 0) || (GMX_GPU_FFT_ONEMATH != 0)
                || (GMX_GPU_FFT_VKFFT != 0) || (GMX_GPU_FFT_CUFFT != 0) || (GMX_GPU_FFT_CLFFT != 0));
    //! Whether this configuration supports running bonded kernels on the device
    static constexpr bool Bonded = GMX_GPU_CUDA || GMX_GPU_SYCL;
    //! Whether this configuration supports running update+LINCS+SETTLE kernels on the device
    static constexpr bool Update = GMX_GPU_CUDA || GMX_GPU_SYCL;
    //! Whether this configuration supports running the direct GPU communication path with thread-MPI
    static constexpr bool ThreadMpiCommunication = GMX_GPU_CUDA;
    //! Whether this configuration supports running the direct GPU communication path with library-MPI
    static constexpr bool LibraryMpiCommunication = GMX_GPU_CUDA || GMX_GPU_SYCL;
    //! Whether this configuration supports running the direct GPU communication path for the current build type
    static constexpr bool MpiCommunication =
            (GMX_THREAD_MPI && ThreadMpiCommunication) || (GMX_LIB_MPI && LibraryMpiCommunication);
    //! Whether this configuration supports running gpu graphs
    static constexpr bool GpuGraph = GMX_HAVE_GPU_GRAPH_SUPPORT;
    //! Whether this configuration supports running gpu pme decomposition
    static constexpr bool PmeDecomposition = (GMX_GPU_CUDA && (GMX_USE_Heffte || GMX_USE_cuFFTMp))
                                             || (GMX_GPU_SYCL && GMX_USE_Heffte)
                                             || (GMX_GPU_HIP && GMX_USE_Heffte);
    static constexpr bool TwoDPmeDecomposition = PmeDecomposition && GMX_GPU_CUDA;
};

} // namespace gmx

#endif // GMX_GPU_UTILS_CAPABILITIES_H
