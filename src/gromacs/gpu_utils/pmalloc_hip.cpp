/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Define functions for host-side memory handling when using AMD with HIP devices.
 *
 *  \author Paul Bauerl <paul.bauer.q@gmail.com>
 *  \author Julio Maia <julio.maia@amd.com>
 */

#include "gmxpre.h"

#include <cstdlib>

#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/utility/cstringutil.h"

#include "pmalloc.h"

/*! Allocates nbytes of page-locked memory.
 *  This memory should always be freed using pfree (or with the page-locked
 *  free functions provied by the HIP library).
 */
void pmalloc(void** h_ptr, size_t nbytes, const DeviceContext* /*context*/)
{
    hipError_t stat;
    char       strbuf[STRLEN];
    int        flag = hipHostMallocDefault;

    if (nbytes == 0)
    {
        *h_ptr = nullptr;
        return;
    }

    gmx::ensureNoPendingDeviceError("Could not allocate page-locked memory.");

    stat = hipHostMalloc(h_ptr, nbytes, flag);
    sprintf(strbuf, "hipHostMalloc of size %d bytes failed", static_cast<int>(nbytes));
    gmx::checkDeviceError(stat, strbuf);
}

/*! Frees page locked memory allocated with pmalloc.
 *  This function can safely be called also with a pointer to a page-locked
 *  memory allocated directly with HIP API calls.
 */
void pfree(void* h_ptr, const DeviceContext* /*context*/)
{
    hipError_t stat;

    if (h_ptr == nullptr)
    {
        return;
    }

    gmx::ensureNoPendingDeviceError("Could not free page-locked memory.");

    stat = hipHostFree(h_ptr);
    gmx::checkDeviceError(stat, "hipHostFree failed");
}

void pmallocSetDefaultDeviceContext(const DeviceContext* /*context*/)
{
    // We don't need context for HIP's pmalloc.
}
void pmallocClearDefaultDeviceContext()
{
    // We don't need context for HIP's pmalloc, so we have nothing to clear.
}
