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
 * Implements classes in analysisdata.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/analysisdata.h"

#include <algorithm>
#include <memory>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

#include "analysisdata-impl.h"

namespace gmx
{

/********************************************************************
 * AnalysisData::Impl
 */

AnalysisData::Impl::Impl()
{
}


AnalysisData::Impl::~Impl()
{
    HandleList::const_iterator i;
    for (i = handles_.begin(); i != handles_.end(); ++i)
    {
        delete *i;
    }
}


/********************************************************************
 * AnalysisData
 */

AnalysisData::AnalysisData()
    : impl_(new Impl)
{
}


AnalysisData::~AnalysisData()
{
    delete impl_;
}


void
AnalysisData::setColumns(int ncol, bool multipoint)
{
    GMX_RELEASE_ASSERT(ncol > 0, "Number of columns must be positive");
    GMX_RELEASE_ASSERT(impl_->handles_.empty(),
                       "Cannot change data dimensionality after creating handles");
    setColumnCount(ncol);
    setMultipoint(multipoint);
    impl_->storage_.setMultipoint(multipoint);
}


AnalysisDataHandle *
AnalysisData::startData(const AnalysisDataParallelOptions &opt)
{
    GMX_RELEASE_ASSERT(impl_->handles_.size() < static_cast<unsigned>(opt.parallelizationFactor()),
                       "Too many calls to startData() compared to provided options");
    if (impl_->handles_.empty())
    {
        notifyDataStart();
        impl_->storage_.setParallelOptions(opt);
        impl_->storage_.startDataStorage(this);
    }
    else if (isMultipoint())
    {
        GMX_THROW(NotImplementedError("Parallelism not supported for multipoint data"));
    }

    std::auto_ptr<AnalysisDataHandle> handle(new AnalysisDataHandle(this));
    impl_->handles_.push_back(handle.get());
    return handle.release();
}


void
AnalysisData::finishData(AnalysisDataHandle *handle)
{
    Impl::HandleList::iterator i;

    i = std::find(impl_->handles_.begin(), impl_->handles_.end(), handle);
    GMX_RELEASE_ASSERT(i != impl_->handles_.end(),
                       "finishData() called for an unknown handle");

    impl_->handles_.erase(i);
    delete handle;

    if (impl_->handles_.empty())
    {
        notifyDataFinish();
    }
}


AnalysisDataFrameRef
AnalysisData::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool
AnalysisData::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataHandle::Impl
 */

AnalysisDataHandle::Impl::Impl(AnalysisData *data)
    : data_(*data), currentFrame_(NULL)
{
}


/********************************************************************
 * AnalysisDataHandle
 */

AnalysisDataHandle::AnalysisDataHandle(AnalysisData *data)
    : impl_(new Impl(data))
{
}


AnalysisDataHandle::~AnalysisDataHandle()
{
    delete impl_;
}


void
AnalysisDataHandle::startFrame(int index, real x, real dx)
{
    GMX_RELEASE_ASSERT(impl_->currentFrame_ == NULL,
                       "startFrame() called twice without calling finishFrame()");
    impl_->currentFrame_ =
        &impl_->data_.impl_->storage_.startFrame(index, x, dx);
}


void
AnalysisDataHandle::setPoint(int col, real y, real dy, bool present)
{
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != NULL,
                       "setPoint() called without calling startFrame()");
    impl_->currentFrame_->setValue(col, y, dy, present);
}


void
AnalysisDataHandle::setPoints(int firstcol, int n, const real *y)
{
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != NULL,
                       "setPoints() called without calling startFrame()");
    for (int i = 0; i < n; ++i)
    {
        impl_->currentFrame_->setValue(firstcol + i, y[i]);
    }
}


void
AnalysisDataHandle::finishPointSet()
{
    GMX_RELEASE_ASSERT(impl_->data_.isMultipoint(),
                       "finishPointSet() called for non-multipoint data");
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != NULL,
                       "finishPointSet() called without calling startFrame()");
    impl_->currentFrame_->finishPointSet();
}


void
AnalysisDataHandle::finishFrame()
{
    GMX_RELEASE_ASSERT(impl_->currentFrame_ != NULL,
                       "finishFrame() called without calling startFrame()");
    int index = impl_->currentFrame_->frameIndex();
    impl_->currentFrame_ = NULL;
    impl_->data_.impl_->storage_.finishFrame(index);
}


void
AnalysisDataHandle::finishData()
{
    // Calls delete this
    impl_->data_.finishData(this);
}

} // namespace gmx
