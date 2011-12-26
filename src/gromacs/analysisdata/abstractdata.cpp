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
 * Implements gmx::AbstractAnalysisData and gmx::AbstractAnalysisDataStored.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/abstractdata.h"

#include <memory>

// Legacy header.
#include "smalloc.h"

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

#include "abstractdata-impl.h"
#include "dataframe.h"
#include "dataproxy.h"

namespace gmx
{

/********************************************************************
 * AbstractAnalysisData::Impl
 */

AbstractAnalysisData::Impl::Impl()
    : _bDataStart(false), _bInData(false), _bInFrame(false),
      _bAllowMissing(true), _nframes(0)
{
}

AbstractAnalysisData::Impl::~Impl()
{
    ModuleList::const_iterator i;
    for (i = _modules.begin(); i != _modules.end(); ++i)
    {
        delete *i;
    }
}


void
AbstractAnalysisData::Impl::presentData(AbstractAnalysisData *data,
                                        AnalysisDataModuleInterface *module)
{
    module->dataStarted(data);
    bool bCheckMissing = _bAllowMissing
        && !(module->flags() & AnalysisDataModuleInterface::efAllowMissing);
    int ncol = data->columnCount();
    for (int i = 0; i < data->frameCount(); ++i)
    {
        real        x, dx;
        const real *y, *dy;
        const bool *present;

        if (!data->getDataWErr(i, &x, &dx, &y, &dy, &present))
        {
            GMX_THROW(APIError("Data not available when module added"));
        }
        if (bCheckMissing && present)
        {
            for (int j = 0; j < ncol; ++j)
            {
                if (!present[j])
                {
                    GMX_THROW(APIError("Missing data not supported by a module"));
                }
            }
        }
        AnalysisDataFrameHeader header(i, x, dx);
        module->frameStarted(header);
        module->pointsAdded(
                AnalysisDataPointSetRef(header, 0, ncol, y, dy, present));
        module->frameFinished();
    }
    if (!_bInData)
    {
        module->dataFinished();
    }
}


/********************************************************************
 * AbstractAnalysisData
 */

AbstractAnalysisData::AbstractAnalysisData()
    : _impl(new Impl()), _ncol(0), _bMultiPoint(false)
{
}


AbstractAnalysisData::~AbstractAnalysisData()
{
    delete _impl;
}


int
AbstractAnalysisData::frameCount() const
{
    return _impl->_nframes;
}


bool
AbstractAnalysisData::getData(int index, real *x, const real **y,
                              const bool **missing) const
{
    return getDataWErr(index, x, 0, y, 0, missing);
}


bool
AbstractAnalysisData::getErrors(int index, real *dx, const real **dy) const
{
    return getDataWErr(index, 0, dx, 0, dy, 0);
}


void
AbstractAnalysisData::addModule(AnalysisDataModuleInterface *module)
{
    std::auto_ptr<AnalysisDataModuleInterface> module_ptr(module);
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }

    if (_impl->_bDataStart)
    {
        GMX_RELEASE_ASSERT(!_impl->_bInFrame,
                           "Cannot add data modules in mid-frame");
        _impl->presentData(this, module);
    }
    if (!(module->flags() & AnalysisDataModuleInterface::efAllowMissing))
    {
        _impl->_bAllowMissing = false;
    }
    _impl->_modules.push_back(module);
    module_ptr.release();
}


void
AbstractAnalysisData::addColumnModule(int col, int span,
                                      AnalysisDataModuleInterface *module)
{
    std::auto_ptr<AnalysisDataModuleInterface> module_ptr(module);
    GMX_RELEASE_ASSERT(col >= 0 && span >= 1 && col + span <= _ncol,
                       "Invalid columns specified for a column module");
    if (_impl->_bDataStart)
    {
        GMX_THROW(NotImplementedError("Cannot add column modules after data"));
    }

    std::auto_ptr<AnalysisDataProxy> proxy(new AnalysisDataProxy(col, span, this));
    proxy->addModule(module_ptr.release());
    addModule(proxy.release());
}


void
AbstractAnalysisData::applyModule(AnalysisDataModuleInterface *module)
{
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }
    GMX_RELEASE_ASSERT(_impl->_bDataStart && !_impl->_bInData,
                       "Data module can only be applied to ready data");

    _impl->presentData(this, module);
}


void
AbstractAnalysisData::setColumnCount(int ncol)
{
    GMX_RELEASE_ASSERT(ncol > 0, "Invalid data column count");
    GMX_RELEASE_ASSERT(_ncol == 0 || _impl->_modules.empty(),
                       "Data column count cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "Data column count cannot be changed after data has been added");
    _ncol = ncol;
}


void
AbstractAnalysisData::setMultipoint(bool multipoint)
{
    GMX_RELEASE_ASSERT(_impl->_modules.empty(),
                       "Data type cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "Data type cannot be changed after data has been added");
    _bMultiPoint = multipoint;
}


/*! \internal
 * This method is not const because the dataStarted() methods of the attached
 * modules can request storage of the data.
 */
void
AbstractAnalysisData::notifyDataStart()
{
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "notifyDataStart() called more than once");
    GMX_RELEASE_ASSERT(_ncol > 0, "Data column count is not set");
    _impl->_bDataStart = _impl->_bInData = true;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        if (_ncol > 1 && !((*i)->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        {
            GMX_THROW(APIError("Data module not compatible with data object properties"));
        }
        (*i)->dataStarted(this);
    }
}


void
AbstractAnalysisData::notifyFrameStart(const AnalysisDataFrameHeader &header) const
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(!_impl->_bInFrame,
               "notifyFrameStart() called while inside a frame");
    GMX_ASSERT(header.index() == _impl->_nframes,
               "Out of order frames");
    _impl->_bInFrame = true;
    _impl->_currHeader = header;
    ++_impl->_nframes;

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->frameStarted(header);
    }
}


void
AbstractAnalysisData::notifyPointsAdd(const AnalysisDataPointSetRef &points) const
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(_impl->_bInFrame, "notifyFrameStart() not called");
    GMX_ASSERT(points.lastColumn() < columnCount(), "Invalid columns");
    GMX_ASSERT(points.frameIndex() == _impl->_currHeader.index(),
               "Points do not correspond to current frame");
    if (!_impl->_bAllowMissing && !points.allPresent())
    {
        GMX_THROW(APIError("Missing data not supported by a module"));
    }

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->pointsAdded(points);
    }
}


void
AbstractAnalysisData::notifyPointsAdd(int firstcol, int n,
                                      const real *y, const real *dy,
                                      const bool *present) const
{
    notifyPointsAdd(AnalysisDataPointSetRef(
            _impl->_currHeader, firstcol, n, y, dy, present));
}


void
AbstractAnalysisData::notifyFrameFinish() const
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(_impl->_bInFrame, "notifyFrameStart() not called");
    _impl->_bInFrame = false;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->frameFinished();
    }
}


void
AbstractAnalysisData::notifyDataFinish() const
{
    GMX_RELEASE_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_RELEASE_ASSERT(!_impl->_bInFrame,
                       "notifyDataFinish() called while inside a frame");
    _impl->_bInData = false;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->dataFinished();
    }
}


/********************************************************************
 * AnalysisDataFrame
 */

AnalysisDataFrame::AnalysisDataFrame()
    : _y(NULL), _dy(NULL), _present(NULL)
{
}

AnalysisDataFrame::~AnalysisDataFrame()
{
    sfree(_y);
    sfree(_dy);
    sfree(_present);
}


void AnalysisDataFrame::allocate(int ncol)
{
    snew(_y, ncol);
    snew(_dy, ncol);
    snew(_present, ncol);
}


/********************************************************************
 * AbstractAnalysisDataStored::Impl
 */

AbstractAnalysisDataStored::Impl::Impl()
    : _nalloc(0), _bStoreAll(false), _nextind(-1)
{
}


AbstractAnalysisDataStored::Impl::~Impl()
{
    FrameList::const_iterator i;
    for (i = _store.begin(); i != _store.end(); ++i)
    {
        delete *i;
    }
}


int
AbstractAnalysisDataStored::Impl::getStoreIndex(int index, int nframes) const
{
    // Check that the requested index is available.
    if ((index < 0 && (-index > _nalloc || -index > nframes))
        || index >= nframes || (index >= 0 && index < nframes - _nalloc))
    {
        return -1;
    }
    // Calculate the index into the storage array.
    if (index < 0)
    {
        index = _nextind + index;
        if (index < 0)
        {
            index += _nalloc;
        }
    }
    else if (nframes > _nalloc)
    {
        index %= _nalloc;
    }
    return index;
}


/********************************************************************
 * AbstractAnalysisDataStored
 */

AbstractAnalysisDataStored::AbstractAnalysisDataStored()
    : _impl(new Impl())
{
}


AbstractAnalysisDataStored::~AbstractAnalysisDataStored()
{
    delete _impl;
}


bool
AbstractAnalysisDataStored::getDataWErr(int index, real *x, real *dx,
                                        const real **y, const real **dy,
                                        const bool **present) const
{
    index = _impl->getStoreIndex(index, frameCount());
    if (index < 0)
    {
        return false;
    }

    // Retrieve the data.
    AnalysisDataFrame *fr = _impl->_store[index];
    if (x)
    {
        *x = fr->_x;
    }
    if (dx)
    {
        *dx = fr->_dx;
    }
    if (y)
    {
        *y = fr->_y;
    }
    if (dy)
    {
        *dy = fr->_dy;
    }
    if (present)
    {
        *present = fr->_present;
    }
    return true;
}


bool
AbstractAnalysisDataStored::requestStorage(int nframes)
{
    GMX_RELEASE_ASSERT(nframes >= -1, "Invalid number of frames requested");
    if (nframes == 0)
    {
        return true;
    }
    GMX_RELEASE_ASSERT(!isMultipoint(), "Storage of multipoint data not supported");

    // Handle the case when everything needs to be stored.
    if (nframes == -1)
    {
        _impl->_bStoreAll = true;
        _impl->_nalloc = 1;
        return true;
    }
    // Check whether an earier call has requested more storage.
    if (_impl->_bStoreAll || nframes < _impl->_nalloc)
    {
        return true;
    }
    // nframes previous frames plus the current one
    _impl->_nalloc = nframes + 1;
    return true;
}


void
AbstractAnalysisDataStored::setMultipoint(bool multipoint)
{
    GMX_RELEASE_ASSERT(_impl->_nalloc == 0 || !multipoint,
                       "Storage of multipoint data not supported");
    AbstractAnalysisData::setMultipoint(multipoint);
}


void
AbstractAnalysisDataStored::startDataStore()
{
    // We first notify any attached modules, because they also might request
    // some storage.
    notifyDataStart();

    int ncol = columnCount();

    // If any storage has been requested, preallocate it.
    if (_impl->_nalloc > 0)
    {
        _impl->_store.resize(_impl->_nalloc);
        for (int i = 0; i < _impl->_nalloc; ++i)
        {
            _impl->_store[i] = new AnalysisDataFrame();
            _impl->_store[i]->allocate(ncol);
        }
        _impl->_nextind = 0;
    }
}


void
AbstractAnalysisDataStored::startNextFrame(const AnalysisDataFrameHeader &header)
{
    // Start storing the frame if needed.
    if (_impl->_nalloc > 0)
    {
        if (_impl->_nextind >= _impl->_nalloc)
        {
            if (_impl->_bStoreAll)
            {
                int ncol = columnCount();

                _impl->_nalloc = _impl->_nextind + 1;
                _impl->_store.resize(_impl->_nalloc);
                for (int i = _impl->_nextind; i < _impl->_nalloc; ++i)
                {
                    _impl->_store[i] = new AnalysisDataFrame();
                    _impl->_store[i]->allocate(ncol);
                }
            }
            else
            {
                _impl->_nextind = 0;
            }
        }

        _impl->_store[_impl->_nextind]->_index = header.index();
        _impl->_store[_impl->_nextind]->_x     = header.x();
        _impl->_store[_impl->_nextind]->_dx    = header.dx();
    }

    // Notify any modules.
    notifyFrameStart(header);
}


void
AbstractAnalysisDataStored::storeThisFrame(const real *y, const real *dy,
                                           const bool *present)
{
    int ncol = columnCount();

    // Store the values if required.
    if (_impl->_nextind >= 0)
    {
        AnalysisDataFrame *fr = _impl->_store[_impl->_nextind];

        for (int i = 0; i < ncol; ++i)
        {
            fr->_y[i] = y[i];
            if (dy)
            {
                fr->_dy[i] = dy[i];
            }
            if (present)
            {
                fr->_present[i] = present[i];
            }
        }
    }

    // Notify modules of new data.
    notifyPointsAdd(0, ncol, y, dy, present);
    // The index needs to be incremented after the notifications to allow
    // the modules to use getData() properly.
    if (_impl->_nextind >= 0)
    {
        ++_impl->_nextind;
    }
    notifyFrameFinish();
}


void
AbstractAnalysisDataStored::storeNextFrame(real x, real dx, const real *y,
                                           const real *dy, const bool *present)
{
    startNextFrame(AnalysisDataFrameHeader(frameCount(), x, dx));
    storeThisFrame(y, dy, present);
}

} // namespace gmx
