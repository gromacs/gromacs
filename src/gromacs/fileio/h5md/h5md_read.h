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

/*! \brief Declarations of utility functions for reading data from H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_READ_H
#define GMX_FILEIO_H5MD_READ_H

#include <hdf5.h>

namespace gmx
{

/*! \brief Read the value at a given frame index in a data set.
 *
 * \tparam ValueType Type of value.
 * \param[in] dataSet Open handle to the data set.
 * \param[in] index Index of frame to read.
 * \param[out] value Single value buffer to read value into.
 */
template<typename ValueType>
void readFrame(const hid_t dataSet, const hsize_t index, ValueType& value);

extern template void readFrame<float>(const hid_t, const hsize_t, float&);

extern template void readFrame<double>(const hid_t, const hsize_t, double&);

extern template void readFrame<int32_t>(const hid_t, const hsize_t, int32_t&);

extern template void readFrame<int64_t>(const hid_t, const hsize_t, int64_t&);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_READ_H
