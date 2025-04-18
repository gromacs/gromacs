
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

/*! \brief Definitions of utility functions for working with H5MD attributes.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_attribute.h"

#include <hdf5.h>

#include <cstring>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/unique_cptr.h"

#include "h5md_error.h"
#include "h5md_guard.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

void setAttribute(const hid_t container, const char* name, const char* value)
{
    throwUponH5mdError(
            H5Aexists(container, name) != 0,
            "Attribute already exists, or an error occured when checking its existance.");

    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Tcopy(H5T_C_S1));
    H5Tset_size(dataType, strlen(value));
    H5Tset_strpad(dataType, H5T_STR_NULLTERM);
    H5Tset_cset(dataType, H5T_CSET_UTF8);

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(
            H5Acreate2(container, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT));

    throwUponH5mdError(H5Awrite(attribute, dataType, value) < 0, "Cannot write attribute.");
}

std::optional<std::string> getAttribute(const hid_t container, const char* name)
{
    if (!H5Aexists(container, name))
    {
        return std::nullopt;
    }
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, name, H5P_DEFAULT));
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Aget_type(attribute));

    /* Make room for string termination as well. */
    size_t            allocationSize = H5Tget_size(dataType) + 1;
    char*             value          = static_cast<char*>(calloc(allocationSize, 1));
    unique_cptr<char> pointerGuard(value);
    throwUponH5mdError(H5Aread(attribute, dataType, pointerGuard.get()) < 0,
                       "Cannot read attribute.");

    return pointerGuard.get();
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
