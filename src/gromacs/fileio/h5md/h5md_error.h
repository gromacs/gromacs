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

/*! \brief Declarations of error utility functions for the H5md module.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_ERROR_H
#define GMX_FILEIO_H5MD_ERROR_H

#include <hdf5.h>

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

/*! \brief Throw a gmx::FileIOError if the HDF5 ID is invalid (H5Iis_valid <= 0).
 *
 * \param[in] id The id to check.
 * \param[in] message The message to throw.
 */
void throwUponInvalidHid(const hid_t id, const std::string& message);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_ERROR_H
