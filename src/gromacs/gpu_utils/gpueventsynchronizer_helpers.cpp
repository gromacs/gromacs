/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *  \brief Helper functions for a GpuEventSynchronizer class.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 * \inlibraryapi
 */

#include "gmxpre.h"

#include "gpueventsynchronizer_helpers.h"

#include "config.h"

#include "gromacs/utility/gmxassert.h"

#if GMX_GPU_CUDA
// Enable event consumption tracking in debug builds, see #3988.
// In OpenCL and SYCL builds, g_useEventConsumptionCounting is constexpr true.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
extern bool g_useEventConsumptionCounting;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
bool g_useEventConsumptionCounting = (CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
                                      || CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_RELWITHDEBINFO);
#endif

namespace gmx::internal
{
void disableGpuEventConsumptionCounting()
{
    GMX_RELEASE_ASSERT(GMX_GPU_CUDA, "Can only be called in CUDA builds");
#if GMX_GPU_CUDA
    /* With threadMPI, we can have a race between different threads setting and reading this flag.
     * However, either all ranks call this function or no one does,
     * so the expected value is the same for all threads,
     * and each thread reads the flag only after callin this function (or deciding no to),
     * so we cannot have any inconsistencies. */
    g_useEventConsumptionCounting = false;
#endif
}
} // namespace gmx::internal
