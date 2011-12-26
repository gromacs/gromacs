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
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

#include "abstractdata-impl.h"
#include "analysisdata-impl.h"

namespace gmx
{

/********************************************************************
 * AnalysisData::Impl
 ********************************************************************/

static bool
frame_index_gtr(AnalysisDataFrame *a, AnalysisDataFrame *b)
{
    return a->_index > b->_index;
}


AnalysisData::Impl::Impl(AnalysisData *data)
    : _data(*data), _pstart(0)
{
}


AnalysisData::Impl::~Impl()
{
    HandleList::const_iterator i;
    for (i = _handles.begin(); i != _handles.end(); ++i)
    {
        delete *i;
    }

    FrameList::const_iterator j;
    for (j = _pending.begin(); j != _pending.end(); ++j)
    {
        delete *j;
    }
}


void
AnalysisData::Impl::addPendingFrame(AnalysisDataFrame *fr)
{
    GMX_ASSERT(fr->_index >= _data.frameCount(),
               "addPendingFrame() called for too old frame");
    size_t pindex = fr->_index - _data.frameCount();
    if (pindex == 0)
    {
        // Just store our frame if it is the next one.
        _data.storeNextFrame(fr->_x, fr->_dx, fr->_y, fr->_dy, fr->_present);
        incrementPStart();
    }
    else
    {
        if (pindex >= _pending.size())
        {
            // TODO: We need to wait until earlier frames are ready...
        }
        // TODO: This is not thread-safe.
        pindex += _pstart;
        if (pindex > _pending.size())
        {
            pindex -= _pending.size();
        }

        int ncol = _data.columnCount();
        _pending[pindex]->_x     = fr->_x;
        _pending[pindex]->_dx    = fr->_dx;
        for (int i = 0; i < ncol; ++i)
        {
            _pending[pindex]->_y[i]       = fr->_y[i];
            _pending[pindex]->_dy[i]      = fr->_dy[i];
            _pending[pindex]->_present[i] = fr->_present[i];
        }
        _pending[pindex]->_index = fr->_index;
    }
    processPendingFrames();
}


void
AnalysisData::Impl::processPendingFrames()
{
    while (_pending[_pstart]->_index != -1)
    {
        AnalysisDataFrame *fr = _pending[_pstart];

        _data.storeNextFrame(fr->_x, fr->_dx, fr->_y, fr->_dy, fr->_present);
        fr->_index = -1;
        incrementPStart();
    }
}


void
AnalysisData::Impl::incrementPStart()
{
    size_t val = _pstart;

    ++val;
    if (val >= _pending.size())
    {
        val -= _pending.size();
    }
    _pstart = val;
}


/********************************************************************
 * AnalysisData
 */

AnalysisData::AnalysisData()
    : _impl(new Impl(this))
{
}


AnalysisData::~AnalysisData()
{
    delete _impl;
}


void
AnalysisData::setColumns(int ncol, bool multipoint)
{
    GMX_RELEASE_ASSERT(ncol > 0, "Number of columns must be positive");
    GMX_RELEASE_ASSERT(_impl->_handles.empty(),
                       "Cannot change data dimensionality after creating handles");
    setColumnCount(ncol);
    setMultipoint(multipoint);
}


AnalysisDataHandle *
AnalysisData::startData(AnalysisDataParallelOptions opt)
{
    if (_impl->_handles.empty())
    {
        startDataStore();
    }
    else if (isMultipoint())
    {
        GMX_THROW(NotImplementedError("Parallelism not supported for multipoint data"));
    }

    std::auto_ptr<AnalysisDataHandle> handle(new AnalysisDataHandle(this));
    _impl->_handles.push_back(handle.get());

    size_t oldSize = _impl->_pending.size();
    _impl->_pending.resize(2 * _impl->_handles.size() - 1);
    Impl::FrameList::iterator i;
    for (i = _impl->_pending.begin() + oldSize; i != _impl->_pending.end(); ++i)
    {
        *i = new AnalysisDataFrame();
        (*i)->allocate(columnCount());
        (*i)->_index = -1;
    }

    return handle.release();
}


void
AnalysisData::finishData(AnalysisDataHandle *handle)
{
    Impl::HandleList::iterator i;

    i = std::find(_impl->_handles.begin(), _impl->_handles.end(), handle);
    GMX_RELEASE_ASSERT(i != _impl->_handles.end(),
                       "finishData() called for an unknown handle");

    _impl->_handles.erase(i);
    delete handle;

    if (_impl->_handles.empty())
    {
        notifyDataFinish();
    }
}


/********************************************************************
 * AnalysisDataHandle::Impl
 */

AnalysisDataHandle::Impl::Impl(AnalysisData *data)
    : _data(*data)
{
    if (!_data.isMultipoint())
    {
        _frame.reset(new AnalysisDataFrame());
        _frame->allocate(_data.columnCount());
    }
}


AnalysisDataHandle::Impl::~Impl()
{
}


/********************************************************************
 * AnalysisDataHandle
 */

AnalysisDataHandle::AnalysisDataHandle(AnalysisData *data)
    : _impl(new Impl(data))
{
}


AnalysisDataHandle::~AnalysisDataHandle()
{
    delete _impl;
}


void
AnalysisDataHandle::startFrame(int index, real x, real dx)
{
    if (_impl->_data.isMultipoint())
    {
        _impl->_data.notifyFrameStart(AnalysisDataFrameHeader(index, x, dx));
    }
    else
    {
        _impl->_frame->_index = index;
        _impl->_frame->_x  = x;
        _impl->_frame->_dx = dx;
        for (int i = 0; i < _impl->_data.columnCount(); ++i)
        {
            _impl->_frame->_y[i]  = 0.0;
            _impl->_frame->_dy[i] = 0.0;
            _impl->_frame->_present[i] = false;
        }
    }
}


void
AnalysisDataHandle::addPoint(int col, real y, real dy, bool present)
{
    if (_impl->_data.isMultipoint())
    {
        _impl->_data.notifyPointsAdd(col, 1, &y, &dy, &present);
    }
    else
    {
        GMX_ASSERT(!_impl->_frame->_present[col],
                   "Data for a column set multiple times");
        _impl->_frame->_y[col] = y;
        _impl->_frame->_dy[col] = dy;
        _impl->_frame->_present[col] = present;
    }
}


void
AnalysisDataHandle::addPoints(int firstcol, int n,
                              const real *y, const real *dy,
                              const bool *present)
{
    if (_impl->_data.isMultipoint())
    {
        _impl->_data.notifyPointsAdd(firstcol, n, y, dy, present);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            addPoint(firstcol + i, y[i], dy ? dy[i] : 0.0,
                     present ? present[i] : true);
        }
    }
}


void
AnalysisDataHandle::finishFrame()
{
    if (_impl->_data.isMultipoint())
    {
        _impl->_data.notifyFrameFinish();
    }
    else
    {
        _impl->_data._impl->addPendingFrame(_impl->_frame.get());
    }
}


void
AnalysisDataHandle::addFrame(int index, real x, const real *y, const real *dy,
                             const bool *present)
{
    addFrame(index, x, 0.0, y, dy, present);
}


void
AnalysisDataHandle::addFrame(int index, real x, real dx,
                             const real *y, const real *dy,
                             const bool *present)
{
    startFrame(index, x, dx);
    addPoints(0, _impl->_data.columnCount(), y, dy, present);
    finishFrame();
}


void
AnalysisDataHandle::finishData()
{
    // Calls delete this
    _impl->_data.finishData(this);
}

} // namespace gmx
