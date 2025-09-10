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

/*! \brief Declarations of H5md utility functions.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_UTIL_H
#define GMX_FILEIO_H5MD_UTIL_H

#include <hdf5.h>

namespace gmx
{

/*! \brief Return whether object \p name exists in the given HDF5 \p container.
 *
 * \param[in] container Handle to container in which to check.
 * \param[in] name      Name of object to look for.
 *
 * \returns True if the object exists, otherwise false.
 */
inline bool objectExists(const hid_t container, const char* name) noexcept
{
    // Return value: <0 = error, 0 = does not exist, >0 = exists
    return H5Oexists_by_name(container, name, H5P_DEFAULT) > 0;
}

/*! \brief Return whether HDF5 \p handle is valid.
 *
 * \param[in] handle HDF5 handle to check.
 *
 * \returns True if the handle is valid, otherwise false.
 */
inline bool handleIsValid(const hid_t handle) noexcept
{
    // Return value: <0 = error, 0 = invalid, >0 = valid
    return H5Iis_valid(handle) > 0;
}

} // namespace gmx

#endif // GMX_FILEIO_H5MD_UTIL_H
