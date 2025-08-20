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

/*! \brief Definitions of H5md data set manipulation routines.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_dataset.h"

#include "gromacs/utility/basedefinitions.h"

#include "h5md_error.h"
#include "h5md_guard.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace
{

using namespace gmx;

/*! \brief Return the per-dimension offset to the data corresponding to a given frame.
 *
 * Since frames are always the major axis the first value of the returned offset is
 * the given frame index. We are only considering the case when all data for that frame
 * is selected, so the offset for all other axis is 0. For example, to reach frame 15
 * of a 3d set we use the offset [15, 0, 0] to access the memory.
 *
 * \param[in] dataSet Handle to data set.
 * \param[in] frameIndex Index of frame.
 * \returns Fixed-size array with the offset.
 */
static DataSetDims getFrameChunkOffset(const hid_t dataSet, const hsize_t frameIndex)
{
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    DataSetDims chunkOffset(H5Sget_simple_extent_ndims(dataSpace), 0);
    chunkOffset[0] = frameIndex;

    return chunkOffset;
}

/* \brief Return the dimensions of the memory chunk for a single frame.
 *
 * Since frames are always the major axis the first value off the chunk dimensions is 1.
 * We only consider the case of selecting all data corresponding to a single frame so
 * all other chunk dimensions are equal to the data set dimensions. For example, the
 * chunk which contains a single frame in a data set with dimensions [30, 150, 3] is
 * [1, 150, 3].
 *
 * \param[in] dataSetDims Dimensions of data set.
 * \returns Fixed-size array with the chunk dimensions.
 */
static DataSetDims getFrameChunkDims(const DataSetDims& dataSetDims)
{
    DataSetDims chunkDims = dataSetDims;
    chunkDims[0]          = 1;

    return chunkDims;
}

} // namespace

namespace gmx
{

template<typename ValueType>
hsize_t getNumFrames(const H5mdDataSetBase<ValueType>& dataSet)
{
    return dataSet.dims()[0];
}

template<typename ValueType>
void setNumFrames(const H5mdDataSetBase<ValueType>& dataSet, hsize_t numFrames)
{
    DataSetDims dataSetDims = dataSet.dims();
    dataSetDims[0]          = numFrames;
    throwUponH5mdError(H5Dset_extent(dataSet.id(), dataSetDims.data()) < 0,
                       "Could not change size of data set");
}

template<typename ValueType>
hid_t getFrameDataSpace(const H5mdDataSetBase<ValueType>& dataSet, hsize_t frameIndex)
{
    const DataSetDims dataSetDims = dataSet.dims();
    const DataSetDims chunkOffset = getFrameChunkOffset(dataSet.id(), frameIndex);
    const DataSetDims chunkDims   = getFrameChunkDims(dataSetDims);

    const hid_t frameDataSpace = H5Dget_space(dataSet.id());
    throwUponH5mdError(
            H5Sselect_hyperslab(
                    frameDataSpace, H5S_SELECT_SET, chunkOffset.data(), nullptr, chunkDims.data(), nullptr)
                    < 0,
            "Could not select hyperslab for given frame index.");

    return frameDataSpace;
}

template<typename ValueType>
hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<ValueType>& dataSet)
{
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    const DataSetDims dataSetDims          = dataSet.dims();
    const DataSetDims chunkDims            = getFrameChunkDims(dataSetDims);

    return H5Screate_simple(H5Sget_simple_extent_ndims(dataSpace), chunkDims.data(), nullptr);
}

template hsize_t getNumFrames(const H5mdDataSetBase<int32_t>&);

template hsize_t getNumFrames(const H5mdDataSetBase<int64_t>&);

template hsize_t getNumFrames(const H5mdDataSetBase<float>&);

template hsize_t getNumFrames(const H5mdDataSetBase<double>&);

template hsize_t getNumFrames(const H5mdDataSetBase<gmx::BasicVector<float>>&);

template hsize_t getNumFrames(const H5mdDataSetBase<gmx::BasicVector<double>>&);

template void setNumFrames(const H5mdDataSetBase<int32_t>&, hsize_t);

template void setNumFrames(const H5mdDataSetBase<int64_t>&, hsize_t);

template void setNumFrames(const H5mdDataSetBase<float>&, hsize_t);

template void setNumFrames(const H5mdDataSetBase<double>&, hsize_t);

template void setNumFrames(const H5mdDataSetBase<gmx::BasicVector<float>>&, hsize_t);

template void setNumFrames(const H5mdDataSetBase<gmx::BasicVector<double>>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<int32_t>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<int64_t>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<float>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<double>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<gmx::BasicVector<float>>&, hsize_t);

template hid_t getFrameDataSpace(const H5mdDataSetBase<gmx::BasicVector<double>>&, hsize_t);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<int32_t>&);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<int64_t>&);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<float>&);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<double>&);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<gmx::BasicVector<float>>&);

template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<gmx::BasicVector<double>>&);

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
