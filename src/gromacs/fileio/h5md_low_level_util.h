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

/*! \brief Declarations of low-level utility functions for working with H5MD HDF5 files.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#ifndef GMX_FILEIO_GMX_H5MD_LOW_LEVEL_UTIL_H
#define GMX_FILEIO_GMX_H5MD_LOW_LEVEL_UTIL_H

#include <optional>
#include <string>

namespace gmx
{

/*! \brief Helper function for printing debug statements.
 *
 * We'd like to embed H5md diagnostic output in an exception
 * object but it can only write it directly to a POSIX stream.
 * Calling this helper method allows some useful information to
 * be passed to the user.
 */
void printHdf5ErrorsDebug();

/* \brief Throw a gmx::FileIOError if there is an error.
 *
 * \param[in] errorExists If true: throw the exception.
 * \param[in] message The message to throw.
 */
void throwUponH5mdError(const bool errorExists, const std::string& message);

/*! \brief Create an H5MD group, and intermediate groups if they do not exist
 *
 * \param[in] container  The ID of the container where the group is located, or should be created.
 * \param[in] name       The name of the group.
 * \returns the ID of the group.
 *
 * \throws FileIOError If the group cannot be created, such as if it already exists.
 */
hid_t createGroup(const hid_t container, const char* name);

/*! \brief Open an existing HDF5 group or create it if it did not exist already.
 *
 * \param[in] container  The ID of the container where the group is located, or should be created.
 * \param[in] name       The name of the group.
 * \returns the ID of the group.
 *
 * \throws FileIOError If the group cannot be found or created.
 */
hid_t openOrCreateGroup(const hid_t container, const char* name);

/*! \brief Set a string attribute value in a group or data set.
 * \param[in] container The ID of the HDF5 container, i.e., group or data set.
 * \param[in] name The name of the attribute.
 * \param[in] value The string to set as attribute value.
 *
 * \throws FileIOError If the parameter could not be set/written or if it already existed
 */
void setAttribute(const hid_t container, const char* name, const char* value);

/*! \brief Get a string attribute value from a group or data set.
 * \param[in] container The ID of the HDF5 container, i.e., group or data set.
 * \param[in] name The name of the attribute.
 * \returns the string value of the attribute, if it was found.
 *
 * \throws FileIOError If the parameter could not be read
 */
std::optional<std::string> getAttribute(const hid_t container, const char* name);

} // namespace gmx

#endif // GMX_FILEIO_GMX_H5MD_LOW_LEVEL_UTIL_H
