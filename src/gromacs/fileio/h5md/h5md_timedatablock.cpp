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

/*! \brief Definitions of H5md time data block manipulation routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "h5md_timedatablock.h"

#include <optional>
#include <utility>

#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

//! \brief Name of value data set inside block group per the H5md specification.
constexpr char c_valueName[] = "value";
//! \brief Name of step data set inside block group per the H5md specification.
constexpr char c_stepName[] = "step";
//! \brief Name of time data set inside block group per the H5md specification.
constexpr char c_timeName[] = "time";
//! \brief Compression algorithm used for step and time data sets.
constexpr H5mdCompression c_stepAndTimeCompression = H5mdCompression::LosslessShuffle;

template<typename ValueType>
H5mdTimeDataBlock<ValueType>::H5mdTimeDataBlock(H5mdFrameDataSet<ValueType>&&     valueDataSet,
                                                H5mdScalarFrameDataSet<int64_t>&& stepDataSet,
                                                std::optional<H5mdScalarFrameDataSet<double>>&& timeDataSet) :
    valueDataSet_{ std::move(valueDataSet) },
    stepDataSet_{ std::move(stepDataSet) },
    timeDataSet_{ std::move(timeDataSet) }
{
    const hsize_t numValues = valueDataSet_.numFrames();
    const hsize_t numSteps  = stepDataSet_.numFrames();
    const hsize_t numTimes = timeDataSet_.has_value() ? timeDataSet_.value().numFrames() : numValues;
    throwUponH5mdError(numValues != numSteps || numValues != numTimes,
                       "Input data sets have different number of frames");
    numFrames_ = static_cast<int64_t>(numValues);
}

template<typename ValueType>
H5mdTimeDataBlock<ValueType>::H5mdTimeDataBlock(const hid_t container, const char* name) :
    H5mdTimeDataBlock(std::move(makeH5mdGroupGuard(openGroup(container, name))))
{
}

template<typename ValueType>
H5mdTimeDataBlock<ValueType>::H5mdTimeDataBlock(decltype(makeH5mdGroupGuard(H5I_INVALID_HID))&& blockGroupAndGuard) :
    // Our input pair matches to: auto [blockGroup, blockGroupGuard] = blockGroupAndGuard,
    // so use the first element (blockGroup) to open all data sets.
    H5mdTimeDataBlock(
            std::move(H5mdFrameDataSet<ValueType>(blockGroupAndGuard.first, c_valueName)),
            std::move(H5mdScalarFrameDataSet<int64_t>(blockGroupAndGuard.first, c_stepName)),
            objectExists(blockGroupAndGuard.first, c_timeName)
                    ? std::make_optional(H5mdScalarFrameDataSet<double>(blockGroupAndGuard.first, c_timeName))
                    : std::nullopt)
{
}

template<typename ValueType>
std::optional<int64_t> H5mdTimeDataBlock<ValueType>::readStepAtIndex(const int64_t frameIndex)
{
    if ((frameIndex >= 0) && (frameIndex < numFrames_))
    {
        int64_t step;
        stepDataSet_.readFrame(frameIndex, &step);
        return step;
    }
    else
    {
        return std::nullopt;
    }
}

template<typename ValueType>
std::optional<double> H5mdTimeDataBlock<ValueType>::readTimeAtIndex(const int64_t frameIndex)
{
    if (timeDataSet_.has_value() && (frameIndex >= 0) && (frameIndex < numFrames_))
    {
        double time;
        timeDataSet_.value().readFrame(frameIndex, &time);
        return time;
    }
    else
    {
        return std::nullopt;
    }
}

template<typename ValueType>
bool H5mdTimeDataBlock<ValueType>::readValueAtIndex(const int64_t frameIndex, ArrayRef<ValueType> values)
{
    if ((frameIndex >= 0) && (frameIndex < numFrames_))
    {
        valueDataSet_.readFrame(frameIndex, values);
        return true;
    }
    else
    {
        return false;
    }
}

template<typename ValueType>
bool H5mdTimeDataBlock<ValueType>::readFrame(const int64_t       frameIndex,
                                             ArrayRef<ValueType> values,
                                             int64_t*            step,
                                             double*             time)
{
    if (frameIndex < numFrames_)
    {
        valueDataSet_.readFrame(frameIndex, values);
        if (step != nullptr)
        {
            stepDataSet_.readFrame(frameIndex, step);
        }
        if (timeDataSet_.has_value() && time != nullptr)
        {
            timeDataSet_->readFrame(frameIndex, time);
        }
        return true;
    }
    else
    {
        return false;
    }
}

template<typename ValueType>
void H5mdTimeDataBlock<ValueType>::writeNextFrame(ArrayRef<const ValueType> values, const int64_t step)
{
    throwUponH5mdError(
            hasTime(),
            "Must not use no-time overload for writeNextFrame when a time data set exists");
    valueDataSet_.writeNextFrame(values);
    stepDataSet_.writeNextFrame(step);
    ++numFrames_;
}

template<typename ValueType>
void H5mdTimeDataBlock<ValueType>::writeNextFrame(ArrayRef<const ValueType> values,
                                                  const int64_t             step,
                                                  const double              time)
{
    valueDataSet_.writeNextFrame(values);
    stepDataSet_.writeNextFrame(step);
    if (timeDataSet_.has_value())
    {
        timeDataSet_.value().writeNextFrame(time);
    }
    ++numFrames_;
}

template<typename ValueType>
H5mdTimeDataBlockBuilder<ValueType>::H5mdTimeDataBlockBuilder(const hid_t        container,
                                                              const std::string& name) :
    blockGroup_{ createGroup(container, name.c_str()) }, valueDataSetBuilder_(blockGroup_, c_valueName)
{
}

template<typename ValueType>
H5mdTimeDataBlockBuilder<ValueType>& H5mdTimeDataBlockBuilder<ValueType>::withCompression(H5mdCompression compression)
{
    valueDataSetBuilder_.withCompression(compression);
    return *this;
}

template<typename ValueType>
H5mdTimeDataBlockBuilder<ValueType>&
H5mdTimeDataBlockBuilder<ValueType>::withFrameDimension(ArrayRef<const hsize_t> dims)
{
    valueDataSetBuilder_.withFrameDimension(dims);
    return *this;
}

template<typename ValueType>
H5mdTimeDataBlockBuilder<ValueType>&
H5mdTimeDataBlockBuilder<ValueType>::withFrameDimension(std::initializer_list<hsize_t> dims)
{
    return withFrameDimension(ArrayRef(dims.begin(), dims.end()));
}

template<typename ValueType>
H5mdTimeDataBlock<ValueType> H5mdTimeDataBlockBuilder<ValueType>::build()
{
    return H5mdTimeDataBlock<ValueType>(
            valueDataSetBuilder_.build(),
            H5mdScalarFrameDataSet<int64_t>{ H5mdFrameDataSetBuilder<int64_t>(blockGroup_, c_stepName)
                                                     .withCompression(c_stepAndTimeCompression)
                                                     .build() },
            H5mdScalarFrameDataSet<double>{ H5mdFrameDataSetBuilder<double>(blockGroup_, c_timeName)
                                                    .withCompression(c_stepAndTimeCompression)
                                                    .build() });
}

template class H5mdTimeDataBlock<int32_t>;
template class H5mdTimeDataBlock<int64_t>;
template class H5mdTimeDataBlock<float>;
template class H5mdTimeDataBlock<double>;
template class H5mdTimeDataBlock<BasicVector<float>>;
template class H5mdTimeDataBlock<BasicVector<double>>;

template class H5mdTimeDataBlockBuilder<int32_t>;
template class H5mdTimeDataBlockBuilder<int64_t>;
template class H5mdTimeDataBlockBuilder<float>;
template class H5mdTimeDataBlockBuilder<double>;
template class H5mdTimeDataBlockBuilder<BasicVector<float>>;
template class H5mdTimeDataBlockBuilder<BasicVector<double>>;

} // namespace gmx
