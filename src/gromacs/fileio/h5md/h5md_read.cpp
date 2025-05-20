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

/*! \brief Definitions of utility functions for reading data from H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "h5md_read.h"

#include <cstring>

#include "gromacs/fileio/h5md/h5md_dataset.h"
#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \brief Read data at index \p index of the \p dataSet into \p readBuffer.
 *
 * Operations common for reading all types of data are performed in this function
 * before the read. This includes checking that the data type matches as well as
 * extracting the HDF5 data space which handles the read data buffers.
 *
 * WARNING: The caller must have verified that the size of \p writeBuffer matches
 * that of a single frame for the \p dataSet.
 *
 * TODO: Once we have a proper data set class it becomes easier to verify this
 * inside this function so we should return to this assumption.
 */
template<typename ValueType, int numDims>
static void readFrameData(const hid_t dataSet, const hsize_t index, const ArrayRef<ValueType> readBuffer)
{
    const hsize_t numFrames = getNumFrames<numDims>(dataSet);
    gmx::throwUponH5mdError(index >= numFrames, "Cannot read frame with frameIndex >= numFrames");

    // When reading the data we must consider that the native data type (on the current compiler
    // target) may not match the native data type that was used when the data was written.
    // The HDF5 library provides this method to determine the correct (native) type to use.
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
    const auto [nativeDataType, nativeDataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(dataType, H5T_DIR_DEFAULT));
    throwUponH5mdError(!valueTypeIsDataType(nativeDataType, makeConstArrayRef(readBuffer)),
                       "Cannot read frame from set with non-matching data type");

    const auto [frameDataSpace, frameDataSpaceGuard] =
            makeH5mdDataSpaceGuard(getFrameDataSpace<numDims>(dataSet, index));
    const auto [memoryDataSpace, memoryDataSpaceGuard] =
            makeH5mdDataSpaceGuard(getFrameMemoryDataSpace<numDims>(dataSet));

    H5Dread(dataSet, nativeDataType, memoryDataSpace, frameDataSpace, H5P_DEFAULT, readBuffer.data());
}

template<typename ValueType>
void readFrame(const hid_t dataSet, const hsize_t index, ValueType& value)
{
    constexpr int numDims = 1;

    readFrameData<ValueType, numDims>(dataSet, index, ArrayRef(&value, &value + 1));
}

template void readFrame<int32_t>(const hid_t, const hsize_t, int32_t&);

template void readFrame<int64_t>(const hid_t, const hsize_t, int64_t&);

template void readFrame<float>(const hid_t, const hsize_t, float&);

template void readFrame<double>(const hid_t, const hsize_t, double&);

} // namespace gmx
