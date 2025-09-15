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

/*! \brief Definitions of H5md frame data set manipulation routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_framedataset.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_util.h" // required for GMX_ASSERT calls

namespace gmx
{

///@{
//! \brief Return the frame dimensions for the data set dimensions.
template<typename ValueType, std::enable_if_t<std::is_arithmetic_v<ValueType>, bool> = true>
static DataSetDims getFrameDims(const DataSetDims& dataSetDims)
{
    return DataSetDims(dataSetDims.cbegin() + 1, dataSetDims.cend());
}

template<typename ValueType,
         std::enable_if_t<std::is_same_v<ValueType, gmx::BasicVector<float>> || std::is_same_v<ValueType, gmx::BasicVector<double>>, bool> = true>
static DataSetDims getFrameDims(const DataSetDims& dataSetDims)
{
    DataSetDims frameDims(dataSetDims.cbegin() + 1, dataSetDims.cend());
    throwUponH5mdError(frameDims.empty(),
                       "Data set dimensions for BasicVector<T> must be at least 1");
    throwUponH5mdError(frameDims.back() != DIM,
                       "Innermost dimension of data set for BasicVector<T> must be 3");
    frameDims.pop_back();

    return frameDims;
}
///@}

static DataSetDims getFullDims(const DataSetDims& dataSetDims)
{
    throwUponH5mdError(dataSetDims.empty(),
                       "Cannot create frame data set for 0-dimensional data set");
    DataSetDims fullDims = dataSetDims;
    if (!fullDims.empty())
    {
        fullDims[0] = 1;
    }
    return fullDims;
}

static hsize_t getNumValues(const DataSetDims& dims)
{
    hsize_t totalProduct = 1;
    for (const hsize_t d : dims)
    {
        totalProduct *= d;
    }
    return totalProduct;
}

template<typename ValueType>
H5mdFrameDataSet<ValueType>::FrameDescription::FrameDescription(const DataSetDims& dataSetDims) :
    dims_{ getFrameDims<ValueType>(dataSetDims) },
    numValues_{ getNumValues(dims_) },
    frameDimsPrimitive_{ getFullDims(dataSetDims) },
    memoryDataSpace_{ H5Screate_simple(frameDimsPrimitive_.size(), frameDimsPrimitive_.data(), nullptr) },
    frameOffset_(dataSetDims.size(), 0)
{
    throwUponInvalidHid(memoryDataSpace_, "Could not create memory data space for frame data");
}

template<typename ValueType>
hid_t H5mdFrameDataSet<ValueType>::FrameDescription::fileDataSpaceForFrame(const hsize_t frameIndex,
                                                                           const hid_t dataSetHandle) noexcept
{
    // To read or write a frame from a data set we must select a rectangular hyperslab
    // of the data. Here we define the offset, which is only along the major axis,
    // which indexes the frames.
    frameOffset_[0] = frameIndex;

    // We now select the hyperslab inside a copy of the data space for the data set.
    // The hyperslab size is given by the dimensions of a single frame, i.e. the data
    // set dimensions with the major axis value = 1 frame.
    const hid_t  fileDataSpace = H5Dget_space(dataSetHandle);
    const herr_t ret           = H5Sselect_hyperslab(
            fileDataSpace, H5S_SELECT_SET, frameOffset_.data(), nullptr, frameDimsPrimitive_.data(), nullptr);
    GMX_ASSERT(ret >= 0, "Could not select hyperslab for given frame index within file");

    return fileDataSpace;
}

template<typename ValueType>
H5mdFrameDataSet<ValueType>::H5mdFrameDataSet(H5mdDataSetBase<ValueType>&& dataSet) :
    Base(std::move(dataSet)),
    extentDimsPrimitive_{ Base::dims() },
    frameDescription_{ extentDimsPrimitive_ },
    numFrames_{ extentDimsPrimitive_[0] } // FrameDescription would throw above for dims.empty()
{
}

template<typename ValueType>
H5mdFrameDataSet<ValueType>::H5mdFrameDataSet(const hid_t container, const char* name) :
    Base(container, name),
    extentDimsPrimitive_{ Base::dims() },
    frameDescription_{ extentDimsPrimitive_ },
    numFrames_{ extentDimsPrimitive_[0] } // FrameDescription would throw above for dims.empty()
{
}

template<typename ValueType>
H5mdFrameDataSet<ValueType>::~H5mdFrameDataSet() noexcept = default;

template<typename ValueType>
H5mdFrameDataSet<ValueType>::H5mdFrameDataSet(H5mdFrameDataSet<ValueType>&&) noexcept = default;

template<typename ValueType>
H5mdFrameDataSet<ValueType>& H5mdFrameDataSet<ValueType>::operator=(H5mdFrameDataSet<ValueType>&&) noexcept = default;

template<typename ValueType>
const DataSetDims& H5mdFrameDataSet<ValueType>::frameDims() const
{
    return frameDescription_.dims();
}

template<typename ValueType>
const DataSetDims& H5mdFrameDataSet<ValueType>::extentForNumFrames(const hsize_t numFrames)
{
    extentDimsPrimitive_[0] = numFrames;
    return extentDimsPrimitive_;
}

template<typename ValueType>
hsize_t H5mdFrameDataSet<ValueType>::numFrames() const noexcept
{
    return numFrames_;
}

template<typename ValueType>
void H5mdFrameDataSet<ValueType>::readFrame(hsize_t index, ArrayRef<ValueType> values)
{
    throwUponH5mdError(index >= numFrames_, "Cannot read frame with index >= numFrames");
    throwUponH5mdError(values.size() != frameDescription_.numValues(),
                       gmx::formatString("Cannot read frame into buffer of incorrect size: "
                                         "size of frame is %llu values but size of buffer is %lu",
                                         static_cast<unsigned long long>(frameDescription_.numValues()),
                                         values.size()));

    const auto [fileDataSpace, fileDataSpaceGuard] =
            makeH5mdDataSpaceGuard(frameDescription_.fileDataSpaceForFrame(index, Base::id()));

    throwUponH5mdError(H5Dread(Base::id(),
                               Base::nativeDataType(),
                               frameDescription_.memoryDataSpace(),
                               fileDataSpace,
                               H5P_DEFAULT,
                               values.data())
                               < 0,
                       "Error reading frame data.");
}

template<typename ValueType>
void H5mdFrameDataSet<ValueType>::writeNextFrame(ArrayRef<const ValueType> values)
{
    throwUponH5mdError(values.size() != frameDescription_.numValues(),
                       gmx::formatString("Cannot write buffer of incorrect size into frame: "
                                         "size of frame is %llu values but size of buffer is %lu",
                                         static_cast<unsigned long long>(frameDescription_.numValues()),
                                         values.size()));

    throwUponH5mdError(
            H5Dset_extent(Base::id(), extentForNumFrames(numFrames_ + 1).data()) < 0,
            gmx::formatString("Could not set the number of frames in the data set to %llu",
                              static_cast<unsigned long long>(numFrames_ + 1)));

    const auto [fileDataSpace, fileDataSpaceGuard] =
            makeH5mdDataSpaceGuard(frameDescription_.fileDataSpaceForFrame(numFrames_, Base::id()));

    if (H5Dwrite(Base::id(), Base::dataType(), frameDescription_.memoryDataSpace(), fileDataSpace, H5P_DEFAULT, values.data())
        < 0)
    {
        // If our write failed we shrink the data set back to its original number of frames before
        // throwing. Ignore any error here, as we are already handling a bigger problem.
        H5Dset_extent(Base::id(), extentForNumFrames(numFrames_).data());
        throwUponH5mdError(true, "Error writing frame data.");
    }

    // Only increment frame index if the write was successful.
    ++numFrames_;
}

template class H5mdFrameDataSet<int32_t>;

template class H5mdFrameDataSet<int64_t>;

template class H5mdFrameDataSet<float>;

template class H5mdFrameDataSet<double>;

template class H5mdFrameDataSet<gmx::BasicVector<float>>;

template class H5mdFrameDataSet<gmx::BasicVector<double>>;

} // namespace gmx
