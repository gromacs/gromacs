/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements classes in datastorage.h and paralleloptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "datastorage.h"

#include <limits>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "datastorage-impl.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataParallelOptions
 */

AnalysisDataParallelOptions::AnalysisDataParallelOptions()
    : parallelizationFactor_(1)
{
}


AnalysisDataParallelOptions::AnalysisDataParallelOptions(int parallelizationFactor)
    : parallelizationFactor_(parallelizationFactor)
{
    GMX_RELEASE_ASSERT(parallelizationFactor >= 1,
                       "Invalid parallelization factor");
}


/********************************************************************
 * AnalysisDataStorageFrame
 */

AnalysisDataStorageFrame::AnalysisDataStorageFrame(AnalysisDataStorage *storage,
                                                   int columnCount, int index)
    : storage_(*storage), header_(index, 0.0, 0.0), values_(columnCount)
{
}


AnalysisDataStorageFrame::~AnalysisDataStorageFrame()
{
}


AnalysisDataPointSetRef
AnalysisDataStorageFrame::currentPoints() const
{
    std::vector<AnalysisDataValue>::const_iterator begin = values_.begin();
    std::vector<AnalysisDataValue>::const_iterator end = values_.end();
    while (begin != end && !begin->isSet())
    {
        ++begin;
    }
    while (end != begin && !(end-1)->isSet())
    {
        --end;
    }
    int firstColumn = (begin != end) ? begin - values_.begin() : 0;
    return AnalysisDataPointSetRef(header_, firstColumn,
                AnalysisDataValuesRef(begin, end));
}


void
AnalysisDataStorageFrame::clearValues()
{
    std::vector<AnalysisDataValue>::iterator i;
    for (i = values_.begin(); i != values_.end(); ++i)
    {
        i->clear();
    }
}


void
AnalysisDataStorageFrame::finishPointSet()
{
    storage_.impl_->notifyPointSet(currentPoints());
    clearValues();
}


/********************************************************************
 * AnalysisDataStorage::Impl
 */

AnalysisDataStorage::Impl::Impl()
    : data_(NULL), bMultipoint_(false),
      storageLimit_(0), pendingLimit_(1), firstFrameLocation_(0), nextIndex_(0)
{
}


AnalysisDataStorage::Impl::~Impl()
{
}


int
AnalysisDataStorage::Impl::columnCount() const
{
    GMX_ASSERT(data_ != NULL, "columnCount() called too early");
    return data_->columnCount();
}


bool
AnalysisDataStorage::Impl::isMultipoint() const
{
    return bMultipoint_;
}


int
AnalysisDataStorage::Impl::firstStoredIndex() const
{
    return frames_[firstFrameLocation_].frame->frameIndex();
}


int
AnalysisDataStorage::Impl::computeStorageLocation(int index) const
{
    if (index < firstStoredIndex() || index >= nextIndex_)
    {
        return -1;
    }
    return index % frames_.size();
}


size_t
AnalysisDataStorage::Impl::endStorageLocation() const
{
    if (storeAll())
    {
        return frames_.size();
    }
    if (frames_[0].frame->frameIndex() == 0 || firstFrameLocation_ == 0)
    {
        return frames_.size() - 1;
    }
    return firstFrameLocation_ - 1;
}


void
AnalysisDataStorage::Impl::extendBuffer(AnalysisDataStorage *storage,
                                        size_t newSize)
{
    frames_.reserve(newSize);
    while (frames_.size() < newSize)
    {
        frames_.push_back(StoredFrame(
            new AnalysisDataStorageFrame(storage, columnCount(), nextIndex_)));
        ++nextIndex_;
    }
    // The unused frame should not be included in the count.
    if (!storeAll())
    {
        --nextIndex_;
    }
}


void
AnalysisDataStorage::Impl::rotateBuffer()
{
    GMX_ASSERT(!storeAll(),
               "No need to rotate internal buffer if everything is stored");
    size_t prevFirst = firstFrameLocation_;
    size_t nextFirst = prevFirst + 1;
    if (nextFirst == frames_.size())
    {
        nextFirst = 0;
    }
    firstFrameLocation_ = nextFirst;
    StoredFrame &prevFrame = frames_[prevFirst];
    prevFrame.status = StoredFrame::eMissing;
    prevFrame.frame->header_ = AnalysisDataFrameHeader(nextIndex_ + 1, 0.0, 0.0);
    prevFrame.frame->clearValues();
    ++nextIndex_;
}


void
AnalysisDataStorage::Impl::notifyPointSet(const AnalysisDataPointSetRef &points)
{
    data_->notifyPointsAdd(points);
}


void
AnalysisDataStorage::Impl::notifyNextFrames(size_t firstLocation)
{
    if (firstLocation != firstFrameLocation_)
    {
        // firstLocation can only be zero here if !storeAll() because
        // firstFrameLocation_ is always zero for storeAll()
        int prevIndex =
            (firstLocation == 0 ? frames_.size() - 1 : firstLocation - 1);
        if (!frames_[prevIndex].isNotified())
        {
            return;
        }
    }
    size_t i = firstLocation;
    size_t end = endStorageLocation();
    while (i != end)
    {
        Impl::StoredFrame &storedFrame = frames_[i];
        if (!storedFrame.isFinished())
        {
            break;
        }
        if (storedFrame.status == StoredFrame::eFinished)
        {
            data_->notifyFrameStart(storedFrame.frame->header());
            data_->notifyPointsAdd(storedFrame.frame->currentPoints());
            data_->notifyFrameFinish(storedFrame.frame->header());
            storedFrame.status = StoredFrame::eNotified;
            if (storedFrame.frame->frameIndex() >= storageLimit_)
            {
                rotateBuffer();
            }
        }
        ++i;
        if (!storeAll() && i >= frames_.size())
        {
            i = 0;
        }
    }
}


/********************************************************************
 * AnalysisDataStorage
 */

AnalysisDataStorage::AnalysisDataStorage()
    : impl_(new Impl())
{
}


AnalysisDataStorage::~AnalysisDataStorage()
{
}


void
AnalysisDataStorage::setMultipoint(bool bMultipoint)
{
    if (bMultipoint && impl_->storageLimit_ > 0)
    {
        GMX_THROW(APIError("Storage of multipoint data not supported"));
    }
    impl_->bMultipoint_ = bMultipoint;
}


void
AnalysisDataStorage::setParallelOptions(const AnalysisDataParallelOptions &opt)
{
    impl_->pendingLimit_ = 2 * opt.parallelizationFactor() - 1;
}


AnalysisDataFrameRef
AnalysisDataStorage::tryGetDataFrame(int index) const
{
    if (impl_->isMultipoint())
    {
        return AnalysisDataFrameRef();
    }
    int storageIndex = impl_->computeStorageLocation(index);
    if (storageIndex == -1)
    {
        return AnalysisDataFrameRef();
    }
    const Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    if (!storedFrame.isAvailable())
    {
        return AnalysisDataFrameRef();
    }
    const Impl::FramePointer &frame = storedFrame.frame;
    return AnalysisDataFrameRef(frame->header(), frame->values_);
}


bool
AnalysisDataStorage::requestStorage(int nframes)
{
    if (impl_->isMultipoint())
    {
        return false;
    }

    // Handle the case when everything needs to be stored.
    if (nframes == -1)
    {
        impl_->storageLimit_ = std::numeric_limits<int>::max();
        return true;
    }
    // Check whether an earlier call has requested more storage.
    if (nframes < impl_->storageLimit_)
    {
        return true;
    }
    impl_->storageLimit_ = nframes;
    return true;
}


void
AnalysisDataStorage::startDataStorage(AbstractAnalysisData *data)
{
    // Data needs to be set before calling extendBuffer()
    impl_->data_ = data;
    setMultipoint(data->isMultipoint());
    if (!impl_->storeAll())
    {
        impl_->extendBuffer(this, impl_->storageLimit_ + impl_->pendingLimit_ + 1);
    }
}


AnalysisDataStorageFrame &
AnalysisDataStorage::startFrame(const AnalysisDataFrameHeader &header)
{
    GMX_ASSERT(header.isValid(), "Invalid header");
    Impl::StoredFrame *storedFrame;
    if (impl_->storeAll())
    {
        size_t size = header.index() + 1;
        if (impl_->frames_.size() < size)
        {
            impl_->extendBuffer(this, size);
        }
        storedFrame = &impl_->frames_[header.index()];
    }
    else
    {
        int storageIndex = impl_->computeStorageLocation(header.index());
        if (storageIndex == -1)
        {
            GMX_THROW(APIError("Out of bounds frame index"));
        }
        storedFrame = &impl_->frames_[storageIndex];
    }
    GMX_RELEASE_ASSERT(!storedFrame->isStarted(),
                       "startFrame() called twice for the same frame");
    GMX_RELEASE_ASSERT(storedFrame->frame->frameIndex() == header.index(),
                       "Inconsistent internal frame indexing");
    storedFrame->status = Impl::StoredFrame::eStarted;
    storedFrame->frame->header_ = header;
    if (impl_->isMultipoint())
    {
        impl_->data_->notifyFrameStart(header);
    }
    return *storedFrame->frame;
}


AnalysisDataStorageFrame &
AnalysisDataStorage::startFrame(int index, real x, real dx)
{
    return startFrame(AnalysisDataFrameHeader(index, x, dx));
}


AnalysisDataStorageFrame &
AnalysisDataStorage::currentFrame(int index)
{
    int storageIndex = impl_->computeStorageLocation(index);
    GMX_RELEASE_ASSERT(storageIndex >= 0, "Out of bounds frame index");
    Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    GMX_RELEASE_ASSERT(storedFrame.isStarted(),
                       "currentFrame() called for frame before startFrame()");
    GMX_RELEASE_ASSERT(!storedFrame.isFinished(),
                       "currentFrame() called for frame after finishFrame()");
    GMX_RELEASE_ASSERT(storedFrame.frame->frameIndex() == index,
                       "Inconsistent internal frame indexing");
    return *storedFrame.frame;
}


void
AnalysisDataStorage::finishFrame(int index)
{
    int storageIndex = impl_->computeStorageLocation(index);
    GMX_RELEASE_ASSERT(storageIndex >= 0, "Out of bounds frame index");
    Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    GMX_RELEASE_ASSERT(storedFrame.isStarted(),
                       "finishFrame() called for frame before startFrame()");
    GMX_RELEASE_ASSERT(!storedFrame.isFinished(),
                       "finishFrame() called twice for the same frame");
    GMX_RELEASE_ASSERT(storedFrame.frame->frameIndex() == index,
                       "Inconsistent internal frame indexing");
    storedFrame.status = Impl::StoredFrame::eFinished;
    if (impl_->isMultipoint())
    {
        // TODO: Check that the last point set has been finished
        impl_->data_->notifyFrameFinish(storedFrame.frame->header());
        if (storedFrame.frame->frameIndex() >= impl_->storageLimit_)
        {
            impl_->rotateBuffer();
        }
    }
    else
    {
        impl_->notifyNextFrames(storageIndex);
    }
}


void
AnalysisDataStorage::finishFrame(const AnalysisDataStorageFrame &frame)
{
    finishFrame(frame.frameIndex());
}

} // namespace gmx
