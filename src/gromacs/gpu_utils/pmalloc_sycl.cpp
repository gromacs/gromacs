/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Define utility routines for SYCL
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/utility/smalloc.h"

#include "pmalloc.h"

/*! \brief Allocates \p nbytes of host memory. Use \c pfree to free memory allocated with this function.
 *
 *  \todo
 *  This function was copied from OpenCL implementation, not tuned for SYCL at all.
 *  Once SYCL2020 is out, might be worthwhile to look into USM and sycl::malloc_host / sycl::aligned_alloc_host.
 *  Overall, it is better to directly use sycl::buffer instead of pinned arrays. But this function
 *  is needed to compile some PME code with SYCL enabled, even if it is never used.
 *
 * \param[in,out]    h_ptr   Pointer where to store the address of the newly allocated buffer.
 * \param[in]        nbytes  Size in bytes of the buffer to be allocated.
 */
void pmalloc(void** h_ptr, size_t nbytes)
{
    /* Need a temporary type whose size is 1 byte, so that the
     * implementation of snew_aligned can cope without issuing
     * warnings. */
    auto** temporary = reinterpret_cast<std::byte**>(h_ptr);

    /* 16-byte alignment inherited from OpenCL and does not sound unreasonable */
    snew_aligned(*temporary, nbytes, 16);
}

/*! \brief Frees memory allocated with pmalloc.
 *
 * \param[in]    h_ptr   Buffer allocated with pmalloc that needs to be freed.
 */
void pfree(void* h_ptr)
{
    if (h_ptr)
    {
        sfree_aligned(h_ptr);
    }
}
