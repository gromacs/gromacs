/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements classes in analysisdata.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/analysisdata.h"

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataHandleImpl
 */

namespace internal
{

/*! \internal \brief
 * Private implementation class for AnalysisDataHandle.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataHandleImpl
{
public:
    //! Creates a handle associated with the given data object.
    explicit AnalysisDataHandleImpl(AnalysisData* data) : data_(*data), currentFrame_(nullptr) {}

    //! The data object that this handle belongs to.
    AnalysisData& data_;
    //! Current storage frame object, or NULL if no current frame.
    AnalysisDataStorageFrame* currentFrame_;
};

} // namespace internal

/********************************************************************
 * AnalysisData::Impl
 */

/*! \internal \brief
 * Private implementation class for AnalysisData.
 *
 * \ingroup module_analysisdata
 */
class AnalysisData::Impl
{
public:
    //! Smart pointer type to manage a data handle implementation.
    typedef std::unique_ptr<internal::AnalysisDataHandleImpl> HandlePointer;
    //! Shorthand for a list of data handles.
    typedef std::vector<HandlePointer> HandleList;

    //! Storage implementation.
    AnalysisDataStorage storage_;
    /*! \brief
     * List of handles for this data object.
     *
     * Note that AnalysisDataHandle objects also contain (raw) pointers
     * to these objects.
     */
    HandleList handles_;
};

/********************************************************************
 * AnalysisData
 */

AnalysisData::AnalysisData() : impl_(new Impl) {}


AnalysisData::~AnalysisData() {}


void AnalysisData::setDataSetCount(int dataSetCount)
{
    GMX_RELEASE_ASSERT(impl_->handles_.empty(),
                       "Cannot change data dimensionality after creating handles");
    AbstractAnalysisData::setDataSetCount(dataSetCount);
}


void AnalysisData::setColumnCount(int dataSet, int columnCount)
{
    GMX_RELEASE_ASSERT(impl_->handles_.empty(),
                       "Cannot change data dimensionality after creating handles");
    AbstractAnalysisData::setColumnCount(dataSet, columnCount);
}


void AnalysisData::setMultipoint(bool bMultipoint)
{
    GMX_RELEASE_ASSERT(impl_->handles_.empty(), "Cannot change data type after creating handles");
    AbstractAnalysisData::setMultipoint(bMultipoint);
}


int AnalysisData::frameCount() const
{
    return impl_->storage_.frameCount();
}


AnalysisDataHandle AnalysisData::startData(const AnalysisDataParallelOptions& opt)
{
    GMX_RELEASE_ASSERT(impl_->handles_.size() < static_cast<unsigned>(opt.parallelizationFactor()),
                       "Too many calls to startData() compared to provided options");
    if (impl_->handles_.empty())
    {
        impl_->storage_.startParallelDataStorage(this, &moduleManager(), opt);
    }

    Impl::HandlePointer handle(new internal::AnalysisDataHandleImpl(this));
    impl_->handles_.push_back(std::move(handle));
    return AnalysisDataHandle(impl_->handles_.back().get());
}


void AnalysisData::finishFrameSerial(int frameIndex)
{
    impl_->storage_.finishFrameSerial(frameIndex);
}


void AnalysisData::finishData(AnalysisDataHandle handle)
{
    Impl::HandleList::iterator i;

    for (i = impl_->handles_.begin(); i != impl_->handles_.end(); ++i)
    {
        if (i->get() == handle.impl_)
        {
            break;
        }
    }
    GMX_RELEASE_ASSERT(i != impl_->handles_.end(), "finishData() called for an unknown handle");

    impl_->handles_.erase(i);

    if (impl_->handles_.empty())
    {
        impl_->storage_.finishDataStorage();
    }
}


AnalysisDataFrameRef AnalysisData::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool AnalysisData::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataHandle
 */

AnalysisDataHandle::AnalysisDataHandle() : impl_(nullptr) {}


AnalysisDataHandle::AnalysisDataHandle(internal::AnalysisDataHandleImpl* impl) : impl_(impl) {}


void AnalysisDataHandle::startFrame(int index, real x, real dx)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ == nullptr,
                       "startFrame() called twice without calling finishFrame()");
    impl_->currentFrame_ = &impl_->data_.impl_->storage_.startFrame(index, x, dx);
}


void AnalysisDataHandle::selectDataSet(int index)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "selectDataSet() called without calling startFrame()");
    impl_->currentFrame_->selectDataSet(index);
}


void AnalysisDataHandle::setPoint(int column, real value, bool bPresent)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "setPoint() called without calling startFrame()");
    impl_->currentFrame_->setValue(column, value, bPresent);
}


void AnalysisDataHandle::setPoint(int column, real value, real error, bool bPresent)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "setPoint() called without calling startFrame()");
    impl_->currentFrame_->setValue(column, value, error, bPresent);
}


void AnalysisDataHandle::setPoints(int firstColumn, int count, const real* values, bool bPresent)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "setPoints() called without calling startFrame()");
    for (int i = 0; i < count; ++i)
    {
        impl_->currentFrame_->setValue(firstColumn + i, values[i], bPresent);
    }
}


void AnalysisDataHandle::finishPointSet()
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->data_.isMultipoint(),
                       "finishPointSet() called for non-multipoint data");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "finishPointSet() called without calling startFrame()");
    impl_->currentFrame_->finishPointSet();
}


void AnalysisDataHandle::finishFrame()
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != nullptr,
                       "finishFrame() called without calling startFrame()");
    AnalysisDataStorageFrame* frame = impl_->currentFrame_;
    impl_->currentFrame_            = nullptr;
    frame->finishFrame();
}


void AnalysisDataHandle::finishData()
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Invalid data handle used");
    // Deletes the implementation pointer.
    impl_->data_.finishData(*this);
    impl_ = nullptr;
}

} // namespace gmx
