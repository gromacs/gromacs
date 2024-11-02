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

/*! \brief Helper functions for managing scope guards in the H5md module.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_guard.h"

#include "config.h"

#if GMX_USE_HDF5

namespace gmx
{

std::pair<hid_t, H5mdGuard> makeH5mdGuard(const hid_t object, H5mdCloseType unboundCloser)
{
    return std::make_pair(object, sg::make_scope_guard(H5mdCloser(object, unboundCloser)));
}

std::pair<hid_t, H5mdGuard> makeH5mdGroupGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Gclose);
}

std::pair<hid_t, H5mdGuard> makeH5mdDataSetGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Dclose);
}

std::pair<hid_t, H5mdGuard> makeH5mdDataSpaceGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Sclose);
}

std::pair<hid_t, H5mdGuard> makeH5mdPropertyListGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Pclose);
}

std::pair<hid_t, H5mdGuard> makeH5mdAttributeGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Aclose);
}

std::pair<hid_t, H5mdGuard> makeH5mdTypeGuard(const hid_t object)
{
    return makeH5mdGuard(object, H5Tclose);
}

} // namespace gmx

#endif // GMX_USE_HDF5
