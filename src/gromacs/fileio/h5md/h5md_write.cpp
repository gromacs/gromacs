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

/*! \brief Definitions of utility functions for writing data to H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "h5md_write.h"

#include "gromacs/fileio/h5md/h5md_dataset.h"
#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \brief Write the given data in \p writeBuffer to the \p index of the \p dataSet.
 *
 * Operations common for writing all types of data are performed in this function
 * before the write. This includes checking that the data type matches as well as
 * extracting the HDF5 data space which handles the written data buffers.
 *
 * WARNING: The caller must have verified that the size of \p writeBuffer matches
 * that of a single frame for the \p dataSet.
 *
 * TODO: Once we have a proper data set class it becomes easier to verify this
 * inside this function so we should return to this assumption.
 */
template<typename ValueType, int numDims>
static void writeFrameData(const hid_t dataSet, const hsize_t index, const ArrayRef<const ValueType> writeBuffer)
{
    const hsize_t numFrames = getNumFrames<numDims>(dataSet);
    gmx::throwUponH5mdError(index >= numFrames, "Cannot write frame with index >= numFrames");

    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
    gmx::throwUponH5mdError(!valueTypeIsDataType(dataType, writeBuffer),
                            "Cannot write frame into set with non-matching data type");

    const auto [frameDataSpace, frameDataSpaceGuard] =
            makeH5mdDataSpaceGuard(getFrameDataSpace<numDims>(dataSet, index));
    const auto [memoryDataSpace, memoryDataSpaceGuard] =
            makeH5mdDataSpaceGuard(getFrameMemoryDataSpace<numDims>(dataSet));

    gmx::throwUponH5mdError(
            H5Dwrite(dataSet, dataType, memoryDataSpace, frameDataSpace, H5P_DEFAULT, writeBuffer.data()) < 0,
            "Error writing frame data.");
}

template<typename ValueType>
void writeFrame(const hid_t dataSet, const hsize_t index, const ValueType& value)
{
    constexpr int numDims = 1;

    writeFrameData<ValueType, numDims>(dataSet, index, ArrayRef(&value, &value + 1));
}

template void writeFrame<int32_t>(const hid_t, const hsize_t, const int32_t&);

template void writeFrame<int64_t>(const hid_t, const hsize_t, const int64_t&);

template void writeFrame<float>(const hid_t, const hsize_t, const float&);

template void writeFrame<double>(const hid_t, const hsize_t, const double&);

} // namespace gmx
