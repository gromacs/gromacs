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

/*! \brief Definitions of utility functions for working with H5MD groups.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_group.h"

#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/utility/basedefinitions.h"

#include "h5md_error.h"
#include "h5md_guard.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

hid_t createGroup(const hid_t container, const char* name)
{
    // create group creation property list
    const auto [propertyList, listGuard] = makeH5mdPropertyListGuard(H5Pcreate(H5P_LINK_CREATE));
    throwUponInvalidHid(propertyList, "Cannot create propertyList when creating group.");

    // Set the option to create intermediate groups if they are missing.
    H5Pset_create_intermediate_group(propertyList, 1);

    hid_t group = H5Gcreate(container, name, propertyList, H5P_DEFAULT, H5P_DEFAULT);
    throwUponInvalidHid(group, "Cannot create group.");

    return group;
}

hid_t openGroup(const hid_t container, const char* name)
{
    const hid_t group = H5Gopen(container, name, H5P_DEFAULT);
    throwUponInvalidHid(group, "Cannot open group.");

    return group;
}

hid_t openOrCreateGroup(const hid_t container, const char* name)
{
    if (objectExists(container, name))
    {
        return openGroup(container, name);
    }
    else
    {
        return createGroup(container, name);
    }
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
