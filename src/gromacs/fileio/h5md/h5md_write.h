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

/*! \brief Declarations of utility functions for writing data to H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_WRITE_H
#define GMX_FILEIO_H5MD_WRITE_H

#include <hdf5.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

/*! \brief Write a value to a given index in a data set.
 *
 * \tparam ValueType Type of value.
 * \param[in] dataSet Open handle to the data set.
 * \param[in] index Index of set to write into.
 * \param[in] value Single value buffer with value to write.
 */
template<typename ValueType>
void writeFrame(const hid_t dataSet, const hsize_t index, const ValueType& value);

/*! \brief Write a list of BasicVectors to a given index in a data set.
 *
 * \tparam ValueType Type of BasicVector.
 * \param[in] dataSet Open handle to the data set.
 * \param[in] index Index of set to write into.
 * \param[in] values Array with BasicVectors to write.
 */
template<typename ValueType>
void writeFrame(const hid_t dataSet, const hsize_t index, const ArrayRef<const BasicVector<ValueType>> values);

extern template void writeFrame(const hid_t, const hsize_t, const int32_t&);

extern template void writeFrame(const hid_t, const hsize_t, const int64_t&);

extern template void writeFrame(const hid_t, const hsize_t, const float&);

extern template void writeFrame(const hid_t, const hsize_t, const double&);

extern template void writeFrame(const hid_t, const hsize_t, const ArrayRef<const BasicVector<float>>);

extern template void writeFrame(const hid_t, const hsize_t, const ArrayRef<const BasicVector<double>>);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_WRITE_H
