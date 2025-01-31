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

/*! \brief Declarations of utility functions for working with H5MD attributes.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_ATTRIBUTE_H
#define GMX_FILEIO_H5MD_ATTRIBUTE_H

#include <optional>
#include <string>

#include "h5md.h"

namespace gmx
{

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

#endif // GMX_FILEIO_H5MD_ATTRIBUTE_H
