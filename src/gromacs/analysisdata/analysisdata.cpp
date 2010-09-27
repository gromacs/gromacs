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

#include <cassert>

#include "gromacs/fatalerror/fatalerror.h"

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


int
AnalysisData::Impl::addPendingFrame(AnalysisDataFrame *fr)
{
    assert(fr->_index >= _data.frameCount());
    size_t pindex = fr->_index - _data.frameCount();
    if (pindex == 0)
    {
        int rc;

        // Just store our frame if it is the next one.
        rc = _data.storeNextFrame(fr->_x, fr->_dx,
                                  fr->_y, fr->_dy, fr->_present);
        if (rc != 0)
        {
            return rc;
        }
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
    return processPendingFrames();
}


int
AnalysisData::Impl::processPendingFrames()
{
    while (_pending[_pstart]->_index != -1)
    {
        AnalysisDataFrame *fr = _pending[_pstart];

        int rc = _data.storeNextFrame(fr->_x, fr->_dx,
                                      fr->_y, fr->_dy, fr->_present);
        if (rc != 0)
        {
            return rc;
        }
        fr->_index = -1;
        incrementPStart();
    }
    return 0;
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
    assert(ncol > 0);
    assert(_impl->_handles.empty());
    setColumnCount(ncol);
    setMultipoint(multipoint);
}


int
AnalysisData::startData(AnalysisDataHandle **handlep,
                        AnalysisDataParallelOptions opt)
{
    if (_impl->_handles.empty())
    {
        int rc = startDataStore();
        if (rc != 0)
        {
            return rc;
        }
    }
    else if (isMultipoint())
    {
        GMX_ERROR(eeNotImplemented,
                  "Parallelism not supported for multipoint data");
    }

    AnalysisDataHandle *handle = new AnalysisDataHandle(this);
    _impl->_handles.push_back(handle);
    *handlep = handle;

    _impl->_pending.resize(2 * _impl->_handles.size() - 1);
    Impl::FrameList::iterator i;
    for (i = _impl->_pending.begin(); i != _impl->_pending.end(); ++i)
    {
        *i = new AnalysisDataFrame();
        (*i)->allocate(columnCount());
        (*i)->_index = -1;
    }

    return 0;
}


int
AnalysisData::finishData(AnalysisDataHandle *handle)
{
    Impl::HandleList::iterator i;

    i = std::find(_impl->_handles.begin(), _impl->_handles.end(), handle);
    assert(i != _impl->_handles.end());

    _impl->_handles.erase(i);
    delete handle;

    if (_impl->_handles.empty())
    {
        return notifyDataFinish();
    }
    return 0;
}


/********************************************************************
 * AnalysisDataHandle::Impl
 */

AnalysisDataHandle::Impl::Impl(AnalysisData *data)
    : _data(*data), _frame(NULL)
{
    if (!_data.isMultipoint())
    {
        _frame = new AnalysisDataFrame();
        _frame->allocate(_data.columnCount());
    }
}


AnalysisDataHandle::Impl::~Impl()
{
    delete _frame;
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


int
AnalysisDataHandle::startFrame(int index, real x, real dx)
{
    if (_impl->_data.isMultipoint())
    {
        return _impl->_data.notifyFrameStart(x, dx);
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
        return 0;
    }
}


int
AnalysisDataHandle::addPoint(int col, real y, real dy, bool present)
{
    if (_impl->_data.isMultipoint())
    {
        return _impl->_data.notifyPointsAdd(col, 1, &y, &dy, &present);
    }
    else
    {
        assert(!_impl->_frame->_present[col]);
        _impl->_frame->_y[col] = y;
        _impl->_frame->_dy[col] = dy;
        _impl->_frame->_present[col] = present;
        return 0;
    }
}


int
AnalysisDataHandle::addPoints(int firstcol, int n,
                              const real *y, const real *dy,
                              const bool *present)
{
    if (_impl->_data.isMultipoint())
    {
        return _impl->_data.notifyPointsAdd(firstcol, n, y, dy, present);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            addPoint(firstcol + i, y[i], dy ? dy[i] : 0.0,
                     present ? present[i] : true);
        }
        return 0;
    }
}


int
AnalysisDataHandle::finishFrame()
{
    if (_impl->_data.isMultipoint())
    {
        return _impl->_data.notifyFrameFinish();
    }
    else
    {
        return _impl->_data._impl->addPendingFrame(_impl->_frame);
    }
}


int
AnalysisDataHandle::addFrame(int index, real x, const real *y, const real *dy,
                             const bool *present)
{
    return addFrame(index, x, 0.0, y, dy, present);
}


int
AnalysisDataHandle::addFrame(int index, real x, real dx,
                             const real *y, const real *dy,
                             const bool *present)
{
    int rc = startFrame(index, x, dx);
    if (rc == 0)
    {
        rc = addPoints(0, _impl->_data.columnCount(), y, dy, present);
    }
    if (rc == 0)
    {
        rc = finishFrame();
    }
    return rc;
}


int
AnalysisDataHandle::finishData()
{
    // Calls delete this
    return _impl->_data.finishData(this);
}

} // namespace gmx
