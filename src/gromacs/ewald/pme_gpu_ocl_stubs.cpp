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

/*! \internal \file
 * \brief Implements stubs of high-level PME GPU functions for OpenCL.
 *
 * \author Gaurav Garg <gaugarg@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/utility/basedefinitions.h"

#include "pme_gpu_grid.h"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types_host_impl.h"

// [[noreturn]] attributes must be added in the common headers, so it's easier to silence the warning here
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

void pmeGpuGridHaloExchange(const PmeGpu* /*pmeGpu*/, gmx_wallcycle* /*wcycle*/)
{
    GMX_THROW(gmx::NotImplementedError("PME decomposition is not implemented in OpenCL"));
}

void pmeGpuGridHaloExchangeReverse(const PmeGpu* /*pmeGpu*/, gmx_wallcycle* /*wcycle*/)
{
    GMX_THROW(gmx::NotImplementedError("PME decomposition is not implemented in OpenCL"));
}

template<bool forward>
void convertPmeGridToFftGrid(const PmeGpu* /*pmeGpu*/,
                             float* /*h_fftRealGrid*/,
                             gmx_parallel_3dfft* /*fftSetup*/,
                             const int /*gridIndex*/)
{
    GMX_THROW(gmx::NotImplementedError("PME decomposition is not implemented in OpenCL"));
}

template<bool forward>
void convertPmeGridToFftGrid(const PmeGpu* /*pmeGpu*/, DeviceBuffer<float>* /*d_fftRealGrid*/, const int /*gridIndex*/)
{
    GMX_THROW(gmx::NotImplementedError("PME decomposition is not implemented in OpenCL"));
}

template void convertPmeGridToFftGrid<true>(const PmeGpu* /*pmeGpu*/,
                                            float* /*h_fftRealGrid*/,
                                            gmx_parallel_3dfft* /*fftSetup*/,
                                            const int /*gridIndex*/);

template void convertPmeGridToFftGrid<false>(const PmeGpu* /*pmeGpu*/,
                                             float* /*h_fftRealGrid*/,
                                             gmx_parallel_3dfft* /*fftSetup*/,
                                             const int /*gridIndex*/);

template void convertPmeGridToFftGrid<true>(const PmeGpu* /*pmeGpu*/,
                                            DeviceBuffer<float>* /*d_fftRealGrid*/,
                                            const int /*gridIndex*/);

template void convertPmeGridToFftGrid<false>(const PmeGpu* /*pmeGpu*/,
                                             DeviceBuffer<float>* /*d_fftRealGrid*/,
                                             const int /*gridIndex*/);

CLANG_DIAGNOSTIC_RESET
